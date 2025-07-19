#' Aggregate local Hessians into pooled Hessian
#'
#' @param hessian_list List of local Hessians produced by dgamlss_get_hessian()
#'
#' @return Pooled hessian
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # Assuming site1_gamlss, site2_gamlss, and pooled_coefs are already defined
#'
#' site1_hess <- dgamlss_get_hessian(site1_gamlss, pooled_coefs)
#' site2_hess <- dgamlss_get_hessian(site2_gamlss, pooled_coefs)
#' pooled_hess <- dgamlss_aggregate_hessian(list(site1_hess, site2_hess))
#' }
dgamlss_aggregate_hessians <- function(hessian_list) {
    return(Reduce("+", hessian_list))
}
