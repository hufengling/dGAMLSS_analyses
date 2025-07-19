#' Title
#'
#' @param local_object Local fit of GAMLSS object. Used to get formatted matrices and coefficients to obtain X^TX matrix and formatted Hessian.
#' @param penalized_parameters String vector containing names of parameters that are to be penalized. X^TX matrices from relevant parameters will be generated.
#' @param pooled_coefs Named list of final pooled coefficients. Names must be "mu", "sigma", "nu", "tau".
#'
#' @return List containing local Hessian, local X^TX matrix for relevant penalized parameters, and local sample size.
#' @export
#'
#' @examples
#' \dontrun{
#' local_object = gamlss(y ~ x + bs_1 + bs_2 + bs_3 + bs_4, ~ x, ~ 1, data = df)
#' dgamlss_get_inference(local_object, penalized_parameters = "mu",
#' pooled_coefs = NULL) # Only gives X^TX matrix.
#' dgamlss_get_inference(local_object, penalized_parameters = "mu",
#' pooled_coefs = pooled_coefs)
#' }
dgamlss_get_inference <- function(local_object,
                                  penalized_parameters = NULL,
                                  pooled_coefs = NULL) {
  if (is.null(pooled_coefs) & is.null(penalized_parameters)) {
    warning("Both pooled_coefs and penalized_parameters are NULL. Only n = number of subjects will be returned.")
  }

  if (!is.null(pooled_coefs)) {
    # family <- as.gamlss.family(local_object$family[1])
    hessian <- dgamlss_get_hessian(local_object, pooled_coefs)
  } else {
    hessian <- NULL
  }

  if (!is.null(penalized_parameters)) {
    xtx_list <- list()
    # rss_list <- list()
    if ("mu" %in% penalized_parameters & "mu" %in% local_object$parameters) {
      xtx_list$mu <- t(local_object$mu.x) %*% local_object$mu.x
      # rss_list$mu <- sum((family$mu.linkfun(local_object$y) -
      #                       local_object$mu.x %*% pooled_coefs$mu)^2)
    }
    if ("sigma" %in% penalized_parameters & "sigma" %in% local_object$parameters) {
      xtx_list$sigma <- t(local_object$sigma.x) %*% local_object$sigma.x
      # rss_list$sigma <- sum((family$sigma.linkfun(local_object$y) -
      #                          local_object$sigma.x %*% pooled_coefs$sigma)^2)
    }
    if ("nu" %in% penalized_parameters & "nu" %in% local_object$parameters) {
      xtx_list$nu <- t(local_object$nu.x) %*% local_object$nu.x
      # rss_list$nu <- sum((family$nu.linkfun(local_object$y) -
      #                       local_object$nu.x %*% pooled_coefs$nu)^2)
    }
    if ("tau" %in% penalized_parameters & "tau" %in% local_object$parameters) {
      xtx_list$tau <- t(local_object$tau.x) %*% local_object$tau.x
      # rss_list$tau <- sum((family$tau.linkfun(local_object$y) -
      #                        local_object$tau.x %*% pooled_coefs$tau)^2)
    }
  } else {
    xtx_list <- NULL
  }

  list(hessian = hessian,
       xtx_list = xtx_list,
       n = length(local_object$y))
}
