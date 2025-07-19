dgamlss_recenter_spline <- function(pooled_coefs, new_intercepts = c(0, 0, 0, 0)) {
  if (new_intercepts[1] != 0) {
    pooled_coefs$mu_coefs <- c(new_intercepts[1], pooled_coefs$mu_coefs - new_intercepts[1])
  }
  if (new_intercepts[2] != 0) {
    pooled_coefs$sigma_coefs <- c(new_intercepts[2], pooled_coefs$sigma_coefs - new_intercepts[2])
  }
  if (new_intercepts[3] != 0) {
    pooled_coefs$nu_coefs <- c(new_intercepts[3], pooled_coefs$nu_coefs - new_intercepts[3])
  }
  if (new_intercepts[4] != 0) {
    pooled_coefs$tau_coefs <- c(new_intercepts[4], pooled_coefs$tau_coefs - new_intercepts[4])
  }

  pooled_coefs$new_intercepts <- new_intercepts
  return(pooled_coefs)
}
