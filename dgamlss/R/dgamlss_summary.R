#' Produce summary of dgamlss_create_summary() output, similar to summary.gamlss()
#'
#' @param object Output from dgamlss_create_summary()
#' @param save to save. Default FALSE
#' @param digits Default 3
#'
#' @return Summary printout
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
dgamlss_summary <- function(object,
                            save = FALSE,
                            digits = max(3, getOption("digits") - 3)) {
  ifWarning <- rep(FALSE, length(object$parameters)) # to create warnings

  pooled_coefs <- object$pooled_coefs
  covmat <- object$vcov
  coef <- covmat$coef
  se <- covmat$se
  tvalue <- coef / se
  pvalue <- 2 * pt(-abs(tvalue), object$df.residual)
  coef.table <- cbind(coef, se, tvalue, pvalue)
  dimnames(coef.table) <- list(names(c(pooled_coefs$mu_coefs,
                                       pooled_coefs$sigma_coefs,
                                       pooled_coefs$nu_coefs,
                                       pooled_coefs$tau_coefs)),
                               c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

  cat("******************************************************************")
  cat("\nFamily: ", deparse(object$family), "\n")
  cat("\nCall: ", deparse(object$call, width.cutoff = 50), "\n", fill = TRUE)
  cat("Fitting method: ", object$method, "\n\n")

  # ================ mu ESTIMATES ========================
  if ("mu" %in% object$parameters) {
    pm <- length(pooled_coefs$mu_coefs)
    p1 <- 1:pm
    cat("------------------------------------------------------------------\n")
    cat("Mu link function: ", object$mu.link)
    cat("\n")
    cat("Mu Coefficients:")
    cat("\n")
    if (is.character(co <- object$contrasts)) {
      cat("  [contrasts: ", apply(cbind(names(co), co), 1, paste, collapse = "="), "]")
    }
    if ("mu" %in% names(object$spline_prefix)) {
      spline_rows <- grep(object$spline_prefix[["mu"]], rownames(coef.table)[p1])
      p1_spline_rows <- p1[spline_rows]
      spline_coefs <- matrix(coef[p1_spline_rows])
      spline_edf <- object$mu.df - (pm - length(spline_coefs))
      spline_vcov <- covmat$vcov[p1_spline_rows, p1_spline_rows]
      F_stat <- (t(spline_coefs) %*% solve(spline_vcov) %*% spline_coefs) / spline_edf #* (object$df.residual - length(spline_coefs)) / (object$df.residual - 1)
      F_p <- pf(F_stat, spline_edf, object$df.residual, lower.tail = F)
      smooth.table <- cbind(spline_edf, F_stat, F_p)
      dimnames(smooth.table) <- list(paste0("bs(", object$spline_prefix[["mu"]], "terms)"),
                                     c("Smooth df", "F value", "Pr(>F)"))

      printCoefmat(coef.table[p1[-spline_rows], , drop = FALSE],
                   digits = digits, signif.stars = TRUE)
      cat("\n")
      cat("Significance of fixed-effect smooth terms:")
      cat("\n")
      printCoefmat(smooth.table, digits = digits, signif.stars = TRUE,
                  P.values = TRUE, has.Pvalue = TRUE)
    } else {
      printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
    }
    cat("\n")
  }
  if ("sigma" %in% object$parameters) {
    ps <- length(pooled_coefs$sigma_coefs)
    p1 <- (pm + 1):(pm + ps)
    # cat("\n")
    cat("------------------------------------------------------------------\n")
    cat("Sigma link function: ", object$sigma.link)
    cat("\n")
    cat("Sigma Coefficients:")
    cat("\n")
    if ("sigma" %in% names(object$spline_prefix)) {
      spline_rows <- grep(object$spline_prefix[["sigma"]], rownames(coef.table)[p1])
      p1_spline_rows <- p1[spline_rows]
      spline_coefs <- matrix(coef[p1_spline_rows])
      spline_edf <- object$sigma.df - (ps - length(spline_coefs))
      spline_vcov <- covmat$vcov[p1_spline_rows, p1_spline_rows]
      F_stat <- (t(spline_coefs) %*% solve(spline_vcov) %*% spline_coefs) / spline_edf #* (object$df.residual - length(spline_coefs)) / (object$df.residual - 1)
      F_p <- pf(F_stat, spline_edf, object$df.residual, lower.tail = F)
      smooth.table <- cbind(spline_edf, F_stat, F_p)
      dimnames(smooth.table) <- list(paste0("bs(", object$spline_prefix[["sigma"]], "terms)"),
                                     c("Smooth df", "F value", "Pr(>F)"))

      printCoefmat(coef.table[p1[-spline_rows], , drop = FALSE],
                   digits = digits, signif.stars = TRUE)
      cat("\n")
      cat("Significance of fixed-effect smooth terms:")
      cat("\n")
      printCoefmat(smooth.table, digits = digits, signif.stars = TRUE,
                   P.values = TRUE, has.Pvalue = TRUE)
    } else {
      printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
    }
    cat("\n")
  }
  if ("nu" %in% object$parameters) {
    pn <- length(pooled_coefs$nu_coefs)
    p1 <- (pm + ps + 1):(pm + ps + pn)
    cat("------------------------------------------------------------------\n")
    cat("Nu link function: ", object$nu.link, "\n")
    cat("Nu Coefficients:")
    cat("\n")
    if ("nu" %in% names(object$spline_prefix)) {
      spline_rows <- grep(object$spline_prefix[["nu"]], rownames(coef.table)[p1])
      p1_spline_rows <- p1[spline_rows]
      spline_coefs <- matrix(coef[p1_spline_rows])
      spline_edf <- object$nu.df - (pn - length(spline_coefs))
      spline_vcov <- covmat$vcov[p1_spline_rows, p1_spline_rows]
      F_stat <- (t(spline_coefs) %*% solve(spline_vcov) %*% spline_coefs) / spline_edf #* (object$df.residual - length(spline_coefs)) / (object$df.residual - 1)
      F_p <- pf(F_stat, spline_edf, object$df.residual, lower.tail = F)
      smooth.table <- cbind(spline_edf, F_stat, F_p)
      dimnames(smooth.table) <- list(paste0("bs(", object$spline_prefix[["nu"]], "terms)"),
                                     c("Smooth df", "F value", "Pr(>F)"))

      printCoefmat(coef.table[p1[-spline_rows], , drop = FALSE],
                   digits = digits, signif.stars = TRUE)
      cat("\n")
      cat("Significance of fixed-effect smooth terms:")
      cat("\n")
      printCoefmat(smooth.table, digits = digits, signif.stars = TRUE,
                   P.values = TRUE, has.Pvalue = TRUE)
    } else {
      printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
    }
    cat("\n")
  }
  if ("tau" %in% object$parameters) {
    pt <- length(pooled_coefs$tau_coefs)
    p1 <- (pm + ps + pn + 1):(pm + ps + pn + pt)
    cat("------------------------------------------------------------------\n")
    cat("Tau link function: ", object$tau.link, "\n")
    cat("Tau Coefficients:")
    cat("\n")
    if ("tau" %in% names(object$spline_prefix)) {
      spline_rows <- grep(object$spline_prefix[["tau"]], rownames(coef.table)[p1])
      p1_spline_rows <- p1[spline_rows]
      spline_coefs <- matrix(coef[p1_spline_rows])
      spline_edf <- object$tau.df - (pt - length(spline_coefs))
      spline_vcov <- covmat$vcov[p1_spline_rows, p1_spline_rows]
      F_stat <- (t(spline_coefs) %*% solve(spline_vcov) %*% spline_coefs) / spline_edf #* (object$df.residual - length(spline_coefs)) / (object$df.residual - 1)
      F_p <- pf(F_stat, spline_edf, object$df.residual, lower.tail = F)
      smooth.table <- cbind(spline_edf, F_stat, F_p)
      dimnames(smooth.table) <- list(paste0("bs(", object$spline_prefix[["tau"]], "terms)"),
                                     c("Smooth df", "F value", "Pr(>F)"))

      printCoefmat(coef.table[p1[-spline_rows], , drop = FALSE], digits = digits, signif.stars = TRUE)
      cat("\n")
      cat("Significance of fixed-effect smooth terms:")
      cat("\n")
      printCoefmat(smooth.table, digits = digits, signif.stars = TRUE,
                   P.values = TRUE, has.Pvalue = TRUE)
    } else {
      printCoefmat(coef.table[p1, , drop = FALSE], digits = digits, signif.stars = TRUE)
    }
    cat("\n")
  }

  cat("------------------------------------------------------------------\n")
  cat("No. of observations in the fit: ", object$noObs, "\n")
  cat("Degrees of Freedom for the fit: ", object$df.fit)
  cat("\n")
  cat("      Residual Deg. of Freedom: ", object$df.residual, "\n")
  cat("     with total communications: ", pooled_coefs$n_communications, "\n \n")
  cat(
    "Global Deviance:    ", object$G.deviance,
    "\n            AIC:    ", object$aic,
    "\n            SBC:    ", object$sbc, "\n"
  ) # format(signif(object$sbc, digits)), "\n")
  cat("******************************************************************")
  cat("\n")

  if (save == TRUE) {
    out <- as.list(environment())
    return(out)
  }
  invisible(coef.table)
}

# ==============================================================================
