#' Create summary object for distributed GAMLSS model
#'
#' @param local_gamlss A GAMLSS object with the same function call as the global model, but created using local site data. Used to pull function call information created by the GAMLSS object.
#' @param pooled_coefs A list containing pooled coefficients for model parameters.
#' @param global_deviance Global deviance equal to sum of site-wise deviances. Site-wise deviances are returned from dgamlss_RS().
#' @param pooled_inference Output from dgamlss_aggregate_inference().
#' @param spline_prefix Named list containing spline prefices. Names must be "mu", "sigma", "nu", and "tau", depending on which parameters include splines.
#'
#' @return A summary object of the distributed GAMLSS model which can be run through dgamlss_summary() and dgamlss_fitted_plot()
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming site1_gamlss and pooled_n are already defined.
#' # Assuming site1_hess and site2_hess has already been sent.
#'
#' last_update <- dgamlss_aggregate_coefs(list(site1_update, site2_update))
#' pooled_hess <- dgamlss_aggregate_hessian(list(site1_hess, site2_hess))
#'
#' final_coef_list <- list(mu.coef = c(10, 1, 2), sigma.coef = c(4, 3),
#' nu.coef = 2, tau.coef = 3) # Declared based on series of last_updates
#' pooled_coefs <- dgamlss_pooled_coefs(current_coef_list, n_communications,
#' last_update$deviance)
#' summary_object <- dgamlss_create_summary(site1_gamlss, pooled_coefs,
#' pooled_hess, pooled_n)
#' dgamlss_summary(summary_object)
#' }
dgamlss_create_summary <- function(local_gamlss,
                                   pooled_coefs,
                                   global_deviance,
                                   pooled_inference,
                                   spline_prefix = NULL) {
  dgamlss_calculate_vcov <- function(pooled_coefs, pooled_hessian) {
    coefs <- pooled_coefs$mu_coefs

    if (!is.null(pooled_coefs$sigma_coefs)) {
      coefs <- c(coefs, pooled_coefs$sigma_coefs)
    }
    if (!is.null(pooled_coefs$nu_coefs)) {
      coefs <- c(coefs, pooled_coefs$nu_coefs)
    }
    if (!is.null(pooled_coefs$tau_coefs)) {
      coefs <- c(coefs, pooled_coefs$tau_coefs)
    }

    varCov <- try(solve(pooled_hessian), silent = TRUE)
    if (class(varCov)[1] != "try-error") {
      rownames(varCov) <- colnames(varCov) <- rownames(pooled_hessian)
      se <- sqrt(diag(varCov))
      corr <- cov2cor(varCov)

      return(list(coef = coefs,
                  hessian = pooled_hessian,
                  vcov = varCov,
                  se = se,
                  cor = corr))
    }
    return(list(coef = coefs))
  }

  summary_out <- list()
  class(summary_out) <- c("gamlss", "gam", "glm", "lm")

  summary_out$mu.formula <- local_gamlss$mu.formula
  summary_out$sigma.formula <- local_gamlss$sigma.formula
  summary_out$nu.formula <- local_gamlss$nu.formula
  summary_out$tau.formula <- local_gamlss$tau.formula

  # summary_out$mu.fv <- local_gamlss$mu.fv
  # summary_out$sigma.fv <- local_gamlss$sigma.fv
  # summary_out$nu.fv <- local_gamlss$nu.fv
  # summary_out$tau.fv <- local_gamlss$tau.fv

  summary_out$parameters <- c("mu")
  names(pooled_coefs$mu_coefs) <- names(local_gamlss$mu.coefficients)
  summary_out$mu.coefficients <- pooled_coefs$mu_coefs

  if ("sigma" %in% (local_gamlss$parameters)) {
    names(pooled_coefs$sigma_coefs) <- names(local_gamlss$sigma.coefficients)
    summary_out$sigma.coefficients <- pooled_coefs$sigma_coefs
    summary_out$parameters <- c(summary_out$parameters, "sigma")
  }

  if ("nu" %in% (local_gamlss$parameters)) {
    names(pooled_coefs$nu_coefs) <- names(local_gamlss$nu.coefficients)
    summary_out$nu.coefficients <- pooled_coefs$nu_coefs
    summary_out$parameters <- c(summary_out$parameters, "nu")
  }

  if ("tau" %in% (local_gamlss$parameters)) {
    names(pooled_coefs$tau_coefs) <- names(local_gamlss$tau.coefficients)
    summary_out$tau.coefficients <- pooled_coefs$tau_coefs
    summary_out$parameters <- c(summary_out$parameters, "tau")
  }

  summary_out$family <- local_gamlss$family
  summary_out$call <- local_gamlss$call
  summary_out$pooled_coefs <- pooled_coefs
  summary_out$vcov <- dgamlss_calculate_vcov(pooled_coefs, pooled_inference$pooled_hessian)
  summary_out$method <- "distributed RS()"

  summary_out$mu.df <- min(length(pooled_coefs$mu_coefs), pooled_inference$edf_list$mu)
  summary_out$sigma.df <- min(length(pooled_coefs$sigma_coefs), pooled_inference$edf_list$sigma)
  summary_out$nu.df <- min(length(pooled_coefs$nu_coefs), pooled_inference$edf_list$nu)
  summary_out$tau.df <- min(length(pooled_coefs$tau_coefs), pooled_inference$edf_list$tau)
  summary_out$df.fit <- summary_out$mu.df + summary_out$sigma.df + summary_out$nu.df + summary_out$tau.df
  summary_out$noObs <- pooled_inference$pooled_n
  summary_out$N <- pooled_inference$pooled_n
  summary_out$df.residual <- summary_out$noObs - summary_out$df.fit

  summary_out$G.deviance <- global_deviance
  summary_out$aic <- summary_out$G.deviance + 2 * summary_out$df.fit
  summary_out$sbc <- summary_out$G.deviance + log(summary_out$noObs) * summary_out$df.fit

  if (!is.null(spline_prefix)) {
    if (!is.list(spline_prefix)) {
      stop("spline_prefix must be a named list with names 'mu', 'sigma', 'nu', or 'tau', depending on which parameters have splines.")
    }
    if (is.null(names(spline_prefix))) {
      if (!(names(spline_prefix) %in% c("mu", "sigma", "nu", "tau"))) {
        stop("spline_prefix must be a named list with names 'mu', 'sigma', 'nu', or 'tau', depending on which parameters have splines.")
      }
    }
    spline_prefix <- lapply(spline_prefix, function(spline_name) {
      if (substr(spline_name, nchar(spline_name), nchar(spline_name)) != "_")
        spline_name <- paste0(spline_name, "_")
      spline_name
    })
    summary_out$spline_prefix <- spline_prefix
  }

  return(summary_out)
}
