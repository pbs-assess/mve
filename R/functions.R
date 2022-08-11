#' Multi-View Embedding (MVE)
#'
#' @param data [data.frame()]
#' @param index [integer()][vector()]
#' @param response [character()]
#' @param lags [list()] whose elements are one named vector of integer lags for
#'   each explanatory variable
#' @param within_row [logical()] forecast response using explanatory values
#'   from within the same row in \code{data}. This is appropriate if the
#'   response is indexed by a generating event but occurs at a later time. For
#'   example sockeye recruitment is indexed by brood year but typically occurs
#'   over the subsequent 3-5 years, so \code{within_row = TRUE} is appropriate.
#'   Note that this excludes the response from the state space reconstruction,
#'   and consequently identifies nearest neighbours by explanatory variables
#'   and their lags, but not by the resulting recruitment.
#' @param n_best [integer()]
#' @param id_cols [character()] colnames in \code{data}
#' @param id_vals [list()] of name-value pairs
#' @param cores [integer()]
#'
#' @importFrom rlang .data
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
mve <- function (data,
                 index,
                 response,
                 lags,
                 within_row = FALSE,
                 n_best = ceiling(sqrt(2^length(unlist(lags)))),
                 id_cols = NULL,
                 id_vals = NULL,
                 cores = 1L) {

  # Check arguments ------------------------------------------------------------


  # Define id columns and values -----------------------------------------------

  # Define id columns
  id_columns <- NULL
  if (length(id_cols) > 0) {
    id_columns <- data[1L, id_cols]
  }
  # Define id values
  id_values <- NULL
  if (length(id_vals) > 0) {
    id_values <- tibble::as_tibble(id_vals)
  }

  # Create subset lags ---------------------------------------------------------

  subset_lags <- create_subset_lags(lags)

  # Apply single view embedding ------------------------------------------------

  if (.Platform$OS.type == "unix") {
    outputs <- parallel::mclapply(
      X = subset_lags,
      FUN = mve::sve,
      data = data,
      index = index,
      response = response,
      within_row = within_row,
      superset = lags,
      mc.cores = cores
    )
  } else {
    outputs <- lapply(
      X = subset_lags,
      FUN = mve::sve,
      data = data,
      index = index,
      response = response,
      within_row = within_row,
      superset = lags
    )
  }

  # Bind rows ------------------------------------------------------------------

  outputs <- outputs %>%
    dplyr::bind_rows(.id = "ssr")

  print(utils::head(outputs))

  # %>%
  #   dplyr::mutate(ssr = as.numeric(.data$ssr)) %>%
  #   dplyr::relocate(.data$ssr, .before = 1)

  # Define results -------------------------------------------------------------

  results <- outputs %>%
    dplyr::arrange(.data$index, .data$rmse) %>%
    dplyr::group_by(.data$index) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$ssr, .data$index) %>%
    dplyr::group_by(.data$ssr) %>%
    dplyr::mutate(lag_rank = dplyr::lag(.data$rank, n = 1L)) %>%
    dplyr::ungroup()

  # Define forecasts -----------------------------------------------------------

  forecasts <- results %>%
    dplyr::filter(.data$lag_rank <= n_best) %>%
    dplyr::arrange(.data$index, .data$lag_rank) %>%
    dplyr::group_by(.data$index) %>%
    dplyr::mutate(mean = mean(.data$forecast, na.omit = TRUE)) %>%
    dplyr::mutate(median = stats::median(.data$forecast, na.rm = TRUE)) %>%
    dplyr::mutate(sd = stats::sd(.data$forecast, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      .data$index,
      .data$observed,
      .data$mean,
      .data$median,
      .data$sd
    ) %>%
    dplyr::distinct() %>%
    tibble::add_row(
      index = 1,
      observed = data[1, response, drop = TRUE],
      .before = 1
    ) %>%
    dplyr::bind_cols(id_values) %>%
    dplyr::relocate(colnames(id_values), .before = 1) %>%
    dplyr::bind_cols(id_columns) %>%
    dplyr::relocate(colnames(id_columns), .before = 1)

  # Define hindsight -----------------------------------------------------------

  hindsight <- results %>%
    dplyr::filter(.data$rank == 1) %>%
    dplyr::arrange(.data$index) %>%
    dplyr::select(-.data$rmse, -.data$mre) %>%
    dplyr::bind_cols(id_values) %>%
    dplyr::relocate(colnames(id_values), .before = 1) %>%
    dplyr::bind_cols(id_columns) %>%
    dplyr::relocate(colnames(id_columns), .before = 1)

  # Define summary -------------------------------------------------------------

  summary <- results %>%
    tidyr::drop_na() %>%
    dplyr::filter(.data$lag_rank <= n_best) %>%
    dplyr::select(
      -(.data$ssr:.data$mre),
      -(.data$rank:.data$lag_rank)
    ) %>%
    dplyr::summarise(dplyr::across(.fns = sum)) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = tidyselect::everything(),
      names_to = "lag",
      values_to = "count"
    ) %>%
    dplyr::mutate(total = n_best * length(index)) %>%
    dplyr::mutate(pct = .data$count / .data$total) %>%
    dplyr::bind_cols(id_values) %>%
    dplyr::relocate(colnames(id_values), .before = 1) %>%
    dplyr::bind_cols(id_columns) %>%
    dplyr::relocate(colnames(id_columns), .before = 1)

  # Return output --------------------------------------------------------------

  return(
    structure(
      list(
        forecasts = forecasts,
        hindsight = hindsight,
        results = results,
        summary = summary
      ),
      class = "mve"
    )
  )
}

#' Single-View Embedding (SVE)
#'
#' @param data [data.frame()]
#' @param index [integer()][vector()]
#' @param response [character()]
#' @param lags [list()] whose elements are one named vector of integer lags for
#'   each explanatory variable
#' @param within_row [logical()] forecast response using explanatory values
#'   from within the same row in \code{data}. This is appropriate if the
#'   response is indexed by a generating event but occurs at a later time. For
#'   example sockeye recruitment is indexed by brood year but typically occurs
#'   over the subsequent 3-5 years, so \code{within_row = TRUE} is appropriate.
#'   Note that this excludes the response from the state space reconstruction,
#'   and consequently identifies nearest neighbours by explanatory variables
#'   and their lags, but not by the resulting recruitment.
#' @param superset [list()]
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
sve <- function (data,
                 index,
                 response,
                 lags,
                 within_row = FALSE,
                 superset = NULL) {

  # Check arguments ------------------------------------------------------------




  # Identify forecast type -----------------------------------------------------

  if (within_row) {

    # Define the state space reconstruction omitting the response --------------

    ssr <- state_space_reconstruction(data = data, response = NULL, lags = lags)

    # Compute state space distances between points -----------------------------

    # - Rows in ssr are points in the SSR
    # - Each row in distances corresponds to a focal point in the SSR
    # - Each column in distances corresponds to a potential neighbour in the SSR
    # - Elements of distances correspond to the distance to a neighbour
    # - NA elements indicate disallowed neighbours for a given focal point

    distances <- state_space_distances(ssr, index)

    # Define the observed vector -----------------------------------------------

    observed <- dplyr::pull(data, response) # Past values used in forecasts

    # Compute centred and scaled forecasts -------------------------------------

    # - Create neighbour index matrix
    # - Create neighbour matrices
    # - Project neighbour matrices
    # - Compute ssr_forecasts vector

    ssr_forecasts <- state_space_forecasts(ssr, distances, within_row, observed)

    # Compute forecasts --------------------------------------------------------

    forecasts <- untransform_forecasts(observed, ssr_forecasts)

    # Exit conditional to prepare output ---------------------------------------

  } else {

    # Define the state space reconstruction including the response -------------

    ssr <- state_space_reconstruction(data, response, lags)

    # Compute state space distances between points -----------------------------

    # - Rows in ssr are points in the SSR
    # - Each row in distances corresponds to a focal point in the SSR
    # - Each column in distances corresponds to a potential neighbour in the SSR
    # - Elements of distances correspond to the distance to a neighbour
    # - NA elements indicate disallowed neighbours for a given focal point

    distances <- state_space_distances(ssr, index)

    # Define the observed vector -----------------------------------------------

    observed <- dplyr::pull(data, response) # Past values used in forecasts

    # Compute centred and scaled forecasts -------------------------------------

    # - Create neighbour index matrix
    # - Create neighbour matrices
    # - Project neighbour matrices
    # - Compute ssr_forecasts vector

    ssr_forecasts <- state_space_forecasts(ssr, distances)

    # Compute forecasts --------------------------------------------------------

    forecasts <- untransform_forecasts(observed, ssr_forecasts)

    # Exit conditional to prepare output ---------------------------------------
  }

  # Prepare output -------------------------------------------------------------

  output <- tibble::tibble(
    index = seq_along(observed),
    observed = observed,
    forecast = forecasts,
    points = count_ssr_points(distances),
    dim = rep(ncol(ssr), nrow(data)),
    rmse = rmse(observed, forecasts, running = TRUE),
    mre = mre(observed, forecasts, running = TRUE),
    superset_columns(data, lags, superset)
  )

  # Return output --------------------------------------------------------------

  return(output)
}

#' Empirical Dynamic Modelling (EDM)
#'
#' @param data [data.frame()]
#' @param index [integer()][vector()]
#' @param response [character()]
#' @param lags [list()] whose elements are one named vector of integer lags for
#'   each explanatory variable
#' @param within_row [logical()] forecast response using explanatory values
#'   from within the same row in \code{data}. This is appropriate if the
#'   response is indexed by a generating event but occurs at a later time. For
#'   example sockeye recruitment is indexed by brood year but typically occurs
#'   over the subsequent 3-5 years, so \code{within_row = TRUE} is appropriate.
#'   Note that this excludes the response from the state space reconstruction,
#'   and consequently identifies nearest neighbours by explanatory variables
#'   and their lags, but not by the resulting recruitment.
#' @param id_cols [character()] colnames in \code{data}
#' @param id_vals [list()] of name-value pairs
#' @param cores [integer()]
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
edm <- function (data,
                 index,
                 response,
                 lags,
                 within_row = FALSE,
                 id_cols = NULL,
                 id_vals = NULL,
                 cores = 1L) {

  # Return forecasts -----------------------------------------------------------

  mve::mve(data = data,
           index = index,
           response = response,
           lags = lags,
           within_row = within_row,
           n_best = 1L,
           id_cols = id_cols,
           id_vals = id_vals,
           cores = cores)
}
