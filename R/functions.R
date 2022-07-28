#' Multi-View Embedding (MVE)
#'
#' @param index [integer()][vector()]
#' @param data [data.frame()]
#' @param response [character()]
#' @param lags [list()] whose elements are one named vector of integer lags for
#'   each explanatory variable
#' @param beyond [logical()]
#' @param n_best [integer()]
#' @param cores [integer()]
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
mve <- function (index,
                 data,
                 response,
                 lags,
                 beyond = FALSE,
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
      index = index,
      data = data,
      response = response,
      beyond = beyond,
      superset = lags,
      mc.cores = cores
    )
  } else {
    outputs <- lapply(
      X = subset_lags,
      FUN = mve::sve,
      index = index,
      data = data,
      response = response,
      beyond = beyond,
      superset = lags
    )
  }
  # Average n best forecasts ---------------------------------------------------

  output <- weight_by_past(outputs, n_best)

  # Return output --------------------------------------------------------------

  return(output)
}


#' Single-View Embedding (SVE)
#'
#' @param index [integer()][vector()]
#' @param data [data.frame()]
#' @param response [character()]
#' @param lags [list()] whose elements are one named vector of integer lags for
#'   each explanatory variable
#' @param beyond [logical()]
#' @param superset [list()]
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
sve <- function (index,
                 data,
                 response,
                 lags,
                 beyond = FALSE,
                 superset = NULL) {

  # Check arguments ------------------------------------------------------------

  # Define the state space reconstruction --------------------------------------

  ssr <- state_space_reconstruction(data, response, lags)

  # Compute state space distances between points -------------------------------
  # TODO: update index and remove buffer

  # - Rows in X are points in the SSR
  # - Each row in X_distance corresponds to a focal point in the SSR
  # - Each column in X_distance corresponds to a potential neighbour in the SSR
  # - Elements of X_distance correspond to distances to neighbours
  # - NA elements indicate disallowed neighbours for a given focal point

  distances <- state_space_distances(ssr, index, buffer)

  # Compute centred and scaled forecasts ---------------------------------------

  # - Create neighbour index matrix
  # - Create neighbour matrices
  # - Project neighbour matrices
  # - Compute X_forecast vector

  ssr_forecasts <- state_space_forecasts(ssr, distances, beyond)

  # Define observed ------------------------------------------------------------

  observed <- c(dplyr::pull(data, response), NA)[seq_along(ssr_forecasts)]

  # Compute forecast -----------------------------------------------------------

  forecast <- untransform_forecasts(observed, ssr_forecasts)

  # Return results -------------------------------------------------------------
  # TODO: update index, compute rmse

  rows <- seq_along(forecast)

  tibble::tibble(
    set = rep(0:1, c(index - 1, nrow(data) - index + 2))[rows],
    time = seq_len(nrow(data) + 1L)[rows],
    points = c(0, as.vector(apply(distances, 1, function (x) sum(!is.na(x)))))[rows],
    dim = rep(ncol(ssr), nrow(data) + 1L)[rows],
    observed = observed,
    forecast = forecast,
    forecast_metrics(observed, forecast, window, metric),
    superset_columns(data, lags, superset)
  )
}

#' Empirical Dynamic Modelling (EDM)
#'
#' @param index [integer()][vector()]
#' @param data [data.frame()]
#' @param response [character()]
#' @param lags [list()] whose elements are one named vector of integer lags for
#'   each explanatory variable
#' @param beyond [logical()]
#' @param cores [integer()]
#'
#' @author Luke A. Rogers
#'
#' @return [list()]
#' @export
#'
edm <- function (index,
                 data,
                 response,
                 lags,
                 beyond = FALSE,
                 cores = 1L) {

  # Return forecasts -----------------------------------------------------------

  mve::mve(index = index,
           data = data,
           response = response,
           lags = lags,
           beyond = beyond,
           n_best = 1L,
           cores = cores)

}
