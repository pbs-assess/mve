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
#' @param cores [integer()]
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
                 cores = 1L) {

  # Check arguments ------------------------------------------------------------



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

  # Average n best forecasts ---------------------------------------------------

  output <- weight_sve_outputs_by_past(outputs, n_best)

  # Return output --------------------------------------------------------------

  return(output)
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
                 cores = 1L) {

  # Return forecasts -----------------------------------------------------------

  mve::mve(data = data,
           index = index,
           response = response,
           lags = lags,
           within_row = within_row,
           n_best = 1L,
           cores = cores)
}
