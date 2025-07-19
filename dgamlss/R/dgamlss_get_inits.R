#' Get initial values based on running gamlss on the coordinating site data.
#'
#' @param local_gamlss gamlss output from local site
#'
#' @return List of initial coefficients pulled from the local_gamlss object. Since the local_gamlss design matrix may be singular even if the pooled design matrix would not be, NAs are substituted by 0s.
#' @export
#'
#' @examples
#' \dontrun{
#' data(abdom)
#' site1_data <- abdom[1:110, ]
#' site2_data <- abdom[111:610, ]
#' site_data <- list(site1_data, site2_data)
#' site1_gamlss <- gamlss(mu.formula = y ~ x1 + x2, sigma.formula = ~ x1,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1_data)
#' initial_values <- dgamlss_get_inits(site1_gamlss)
#' }
dgamlss_get_inits <- function(local_gamlss) {
  mu_inits <- local_gamlss$mu.coefficients
  sigma_inits <- local_gamlss$sigma.coefficients
  nu_inits <- local_gamlss$nu.coefficients
  tau_inits <- local_gamlss$tau.coefficients

  inits_list <- list(mu_inits = mu_inits,
    sigma_inits = sigma_inits,
    nu_inits = nu_inits,
    tau_inits = tau_inits)

  inits_list <- lapply(inits_list, function(init) {
    if (is.null(init)) {
      return(NULL)
    }
    init[is.na(init)] <- 0
    init
  })

  return(inits_list)
}
