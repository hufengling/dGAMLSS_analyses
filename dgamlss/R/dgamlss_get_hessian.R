#' Get Hessian matrix for a GAMLSS object
#'
#' @param local_object A GAMLSS object with the same function call as the global model, but created using local site data. Used to pull likelihood information created by the GAMLSS object.
#' @param pooled_coefs A list containing converged pooled coefficients for model parameters. Output from dgamlss_pooled_coefs()
#'
#' @import gamlss
#' @import stats
#' @return A Hessian matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming site1_gamlss, site2_gamlss, and pooled_coefs are already defined
#'
#' site1_hess <- dgamlss_get_hessian(site1_gamlss, pooled_coefs)
#' site2_hess <- dgamlss_get_hessian(site2_gamlss, pooled_coefs)
#' pooled_hess <- dgamlss_aggregate_hessian(list(site1_hess, site2_hess))
#' }
dgamlss_get_hessian <- function(local_object, pooled_coefs) {
  replace_local_coefs <- function(local_gamlss, pooled_coefs) {
    names(pooled_coefs$mu) <- names(local_gamlss$mu.coefficients)
    local_gamlss$mu.coefficients <- pooled_coefs$mu

    if (!is.null(pooled_coefs$sigma)) {
      names(pooled_coefs$sigma) <- names(local_gamlss$sigma.coefficients)
      local_gamlss$sigma.coefficients <- pooled_coefs$sigma
    }

    if (!is.null(pooled_coefs$nu)) {
      names(pooled_coefs$nu) <- names(local_gamlss$nu.coefficients)
      local_gamlss$nu.coefficients <- pooled_coefs$nu
    }

    if (!is.null(pooled_coefs$tau)) {
      names(pooled_coefs$tau) <- names(local_gamlss$tau.coefficients)
      local_gamlss$tau.coefficients <- pooled_coefs$tau
    }

    return(local_gamlss)
  }

  ## end of local ----------------------------------------------------------------
  local_object <- replace_local_coefs(local_object, pooled_coefs)

  if (!is.gamlss(local_object)) {
    stop(paste("This is not an gamlss object", "\n", ""))
  }

  coefBeta <- list()
  for (i in local_object$par)
  {
    if (length(eval(parse(text = paste(paste("local_object$", i, sep = ""), ".fix==TRUE", sep = "")))) != 0) {
      ff <- eval(parse(text = paste(paste(paste(local_object$family[1], "()$", sep = ""), i, sep = ""), ".linkfun", sep = "")))
      fixvalue <- ff(fitted(local_object, i)[1])
      names(fixvalue) <- paste("fixed", i, sep = " ")
      coefBeta <- c(coefBeta, fixvalue)
    } else {
      if (!is.null(unlist(attr(terms(formula(local_object, i), specials = .gamlss.sm.list), "specials")))) {
        warning(paste("Additive terms exists in the ", i, "formula. \n ",
                      "Standard errors for the linear terms maybe are not appropriate"))
      }
      nonNAcoef <- !is.na(coef(local_object, i))
      coefBeta <- c(coefBeta, coef(local_object, i)[nonNAcoef])
    }
  }

  betaCoef <- unlist(coefBeta)
  like.fun <- dgamlss_gen_likelihood(local_object)
  hess <- optimHess(betaCoef, like.fun)
  return(hess)
}
