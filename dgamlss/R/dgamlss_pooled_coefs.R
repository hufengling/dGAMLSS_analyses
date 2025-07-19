#' Create object containing pooled coefficients for distributed GAMLSS model
#'
#' @param current_coefs List containing updated pooled coefficients for mu, sigma, nu, and tau. Can use output from dgamlss_send_coefs() or manually create named list with names "mu.coef", "sigma.coef", "nu.coef", "tau.coef"
#' @param n_communications Total number of communications during model fitting.
#' @param global_deviance Global deviance of the distributed model. Output from dgamlss_aggregate_sites()
#'
#' @return A list containing pooled coefficients, number of communications, and global deviance. For use with dgamlss_create_s
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming site data from dgamlss_RS() have already been collected.
#'
#' n_communications <- 15 # Count total number of communications
#' last_update <- dgamlss_aggregate_sites(list(site1, site2))
#' current_coef_list <- list(mu.coef = c(10, 1, 2), sigma.coef = c(4, 3), nu.coef = 2, tau.coef = NULL)
#' pooled_coefs <- dgamlss_pooled_coefs(current_coef_list, n_communications, last_update$deviance)
#' }
dgamlss_pooled_coefs <- function(current_coefs, n_communications, global_deviance) {
  output <- list(mu_coefs = current_coefs$mu.coef,
                 sigma_coefs = current_coefs$sigma.coef,
                 nu_coefs = current_coefs$nu.coef,
                 tau_coefs = current_coefs$tau.coef,
                 n_communications = n_communications,
                 global_deviance = global_deviance)
  class(output) <- c("dgamlss_pooled_coefs")
  return(output)
}
