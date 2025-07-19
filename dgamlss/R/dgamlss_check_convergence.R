#' Check convergence of GAMLSS model updates
#'
#' @param current_update List containing current model coefficients. Note that the current_update list must have an item with name to.update that specifies which coefficients to check for convergence.
#' @param previous_update List containing previous model coefficients.
#' @param coef_crit Convergence criterion in terms of proportion of change for coefficients between updates. Default is 0.05, so all previous coefficients must be within 5% of the current coefficients to pass. Higher coef_crit may allow for faster convergence, but with larger error with respect to the pooled analysis.
#'
#' @return TRUE if all coefficients have converged, FALSE otherwise.
#' @export
#'
#' @examples
#' # Simulate current and previous updates
#' current_update <- dgamlss_send_coefs(to_update = "mu", mu.coef = c(0.1, 0.2, 0.3))
#' previous_update <- dgamlss_send_coefs(to_update = "mu", mu.coef = c(0.09, 0.19, 0.29))
#'
#' # Check convergence
#' is_converged <- dgamlss_check_convergence(current_update, previous_update, coef_crit = 0.05)

dgamlss_check_convergence <- function(current_update, previous_update, coef_crit = 0.05) {
  get_percent_change <- function(current_coef, previous_coef) {
    abs((current_coef - previous_coef) / current_coef)
  }
  if (current_update$to_update == "mu") {
    is_converged <- all(get_percent_change(current_update$mu.coef, previous_update$mu.coef) < coef_crit)
  }
  if (current_update$to_update == "sigma") {
    is_converged <- all(get_percent_change(current_update$sigma.coef, previous_update$sigma.coef) < coef_crit)
  }
  if (current_update$to_update == "nu") {
    is_converged <- all(get_percent_change(current_update$nu.coef, previous_update$nu.coef) < coef_crit)
  }
  if (current_update$to_update == "tau") {
    is_converged <- all(get_percent_change(current_update$tau.coef, previous_update$tau.coef) < coef_crit)
  }
  return(is_converged)
}
