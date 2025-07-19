#' Select optimal lambda using GAIC
#'
#' @param lambda_vec Vector of lambda values
#' @param proposed_coef_list List of fitted coefficients for each lambda value. Length must be equal to length(lambda_vec). Output from dgamlss_aggregate_coef()
#' @param site_ssr_list List of vectors of sum of squared residuals for each site for each lambda value. Length of site_ssr_list should be equal to number of sites. Length of SSR vectors within site_ssr_list should be length(lambda_vec). Output from dgamlss_RS() with get_ssr = TRUE.
#' @param site_xtx_list List of X^TX matrix for the given parameter from each site. Length of site_xtx_list should be equal to number of sites. Output from dgamlss_get_inference() with pooled_coefs = NULL.
#' @param penalty_matrix Penalty matrix for the given parameter. Length of each side of penalty matrix should be equal to number of coefficients for that parameter. Same size as X^TX. Output from generate_penalty_matrix().
#' @param k GAIC multiplier. Default is k = 2 for standard AIC. BIC is k = ln(n)
#'
#' @return List of optimal lambda for the given update along with corresponding GAIC and effective degrees of freedom (EDF)
#' @export
#'
#' @examples
#' \dontrun{
#' dgamlss_select_lambda(seq(0.1, 10, by = 0.1), proposed_coef_list = dgamlss_aggregate_coef_output,
#' site_ssr_list = list(site1 = site1_ssr_vec, site2 = site2_ssr_vec),
#' site_xtx_list = list(site1 = site1_xtx, site2 = site2_xtx),
#' penalty_matrix = penalty_matrix_list$mu,
#' k = 2)
#' }
dgamlss_select_lambda <- function(lambda_vec,
                                  proposed_coef_list,
                                  site_ssr_list,
                                  site_xtx_list,
                                  penalty_matrix,
                                  k = 2) {
  if (is.null(site_xtx_list)) {
    stop("If a penalty is being used, xtx_site_list must be supplied in order to estimate effective degrees of freedom for each given lambda. Please run dgamlss_get_inference() one time at local sites -- xtx_site_list will not change between communication rounds and therefore does not need to be updated.")
  }

  xtx <- Reduce("+", site_xtx_list)
  site_ssr_df <- do.call(rbind, site_ssr_list)

  edf_vec <- numeric(length(lambda_vec))
  gaic_vec <- numeric(length(lambda_vec))
  for (i in 1:length(lambda_vec)) {
    edf_vec[i] <- sum(diag(solve(xtx + lambda_vec[[i]] * penalty_matrix) %*% xtx))
    gaic_vec[i] <- sum(site_ssr_df[, i]) + k * edf_vec[i]
  }

  print(gaic_vec)
  optimal_index <- which.min(gaic_vec)

  return(list(lambda_opt = lambda_vec[optimal_index],
              gaic_opt = gaic_vec[optimal_index],
              edf_opt = edf_vec[optimal_index],
              index = optimal_index))
}
