#' Send pooled coefficients
#'
#' @param to_update One of "mu", "sigma", "nu", "tau" specifying which parameter the local sites should send relevant summary statitics back for.
#' @param mu.coef Vector of current pooled mu coefficients
#' @param sigma.coef Vector of current pooled sigma coefficients
#' @param nu.coef Vector of current pooled nu coefficients
#' @param tau.coef Vector of current pooled tau coefficients
#'
#' @return List of pooled coefficients that can be used with dgamlss_RS
#' @export
#'
#' @examples
#' \dontrun{
#' # To update the "mu" coefficients
#'
#' data(abdom)
#' site1_data <- abdom[1:110, ]
#' site2_data <- abdom[111:610, ]
#' current_update <- dgamlss_send_coefs(to_update = "mu", mu.coef = c(10, 1, 2),
#' sigma.coef = c(4, 3), nu.coef = 2, tau.coef = 1)
#' site1_update <- dgamlss_RS(formula = y ~ x1 + x2, sigma.formula = ~ x,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1_data,
#' dgamlss.update = current_update)
#' site2_update <- dgamlss_RS(formula = y ~ x1 + x2, sigma.formula = ~ x,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site2_data,
#' dgamlss.update = current_update)
#' pooled_update <- dgamlss_aggregate_coef(list(site1_update, site2_update))
#' }
dgamlss_send_coefs <- function(to_update, #  mu, sigma, nu, tau
                               mu.coef = NULL,
                               sigma.coef = NULL,
                               nu.coef = NULL,
                               tau.coef = NULL) {
  return(list(to_update = to_update,
              mu.coef = as.vector(mu.coef),
              sigma.coef = as.vector(sigma.coef),
              nu.coef = as.vector(nu.coef),
              tau.coef = as.vector(tau.coef)))
}

