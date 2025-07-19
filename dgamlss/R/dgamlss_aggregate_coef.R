#' Aggregate local site data into pooled coefficient update.
#'
#' @param site_list List of sumarry statistic updates from local sites.
#' @param fixed_lambda If fixed penalty smooth is desired, can provide one fixed lambda. Default NULL to allow for empirical estimation of optimal lambda.
#' @param central_lambda If empirical estimation of optimal smoothing is desired, can provide central lambda from which to begin search.
#' @param search_type Options for fine or coarse grid of lambdas to check around the central_lambda.
#' @param penalty_matrix Penalty matrix for the given parameter being fit.
#'
#' @return List of new coefficient update as well as global deviance calculated at the previous coefficient. Note that global deviance is "lagged" by one communication round.
#' @export
#'
#' @examples
#' # To update the "mu" coefficients
#'
#' \dontrun{
#' current_update <- dgamlss_send_coefs(to_update = "mu", mu.coef = c(10, 1, 2),
#' sigma.coef = c(4, 3), nu.coef = 2, tau.coef = 1)
#' site1_update <- dgamlss_RS(formula = y ~ x1 + x2, sigma.formula = ~ x,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1,
#' dgamlss.update = current_update)
#' site1_update <- dgamlss_RS(formula = y ~ x1 + x2, sigma.formula = ~ x,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1,
#' dgamlss.update = current_update)
#' pooled_update <- dgamlss_aggregate_coef(list(site1_hess, site2_hess))
#' }
dgamlss_aggregate_coef <- function(site_list,
                                   fixed_lambda = NULL,
                                   central_lambda = NULL,
                                   search_type = c("coarse", "fine"),
                                   penalty_matrix = NULL) {
  search_type <- match.arg(search_type, c("coarse", "fine"))
  n_sites <- length(site_list)
  xtwx_list <- vector("list", n_sites)
  xtwy_list <- vector("list", n_sites)
  deviance_list <- vector("list", n_sites)

  for (i in 1:n_sites) {
    xtwx_list[[i]] <- site_list[[i]]$xtwx
    xtwy_list[[i]] <- site_list[[i]]$xtwy
    deviance_list[[i]] <- site_list[[i]]$deviance
  }

  if (!is.null(penalty_matrix)) {
    if (nrow(penalty_matrix) != nrow(xtwx_list[[1]])) {
      stop(paste0("Penalty matrix given by dgamlss_bs() only covers the spline terms, but does not cover fixed effects. Use 'generated_combined_penalty_matrix()' function to generate desired penalty matrix. Penalty matrix size must be ", nrow(xtwx_list[[1]])))
    }

    if (is.null(fixed_lambda)) {
      if (is.null(central_lambda)) {
        #warning("Setting central_lambda to default of 10.")
        central_lambda <- 10
      }

      if (search_type == "fine") {
        lambda_vec <- 10^(seq(log10(central_lambda) - 3, log10(central_lambda) + 3, by = 0.05))
      } else if (search_type == "coarse") {
        lambda_vec <- 10^(seq(log10(central_lambda) - 6, log10(central_lambda) + 6, by = 0.1))
      }

      proposed_coef_list <- lapply(lambda_vec, \(lambda) solve(Reduce("+", xtwx_list) + lambda * penalty_matrix) %*% Reduce("+", xtwy_list))

      return(list(proposed_coef_list = proposed_coef_list,
                  lambda_vec = lambda_vec))
    } else {
      return(solve(Reduce("+", xtwx_list) + fixed_lambda * penalty_matrix) %*% Reduce("+", xtwy_list))
    }
  }

  return(list(coef = solve(Reduce("+", xtwx_list)) %*% Reduce("+", xtwy_list),
              deviance = Reduce("+", deviance_list)))
}
