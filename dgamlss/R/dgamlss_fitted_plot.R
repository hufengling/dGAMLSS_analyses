#' Produce fitted plot of dgamlss_create_summary() output, similar to fittedPlot()
#'
#' @param object Output from dgamlss_create_summary()
#' @param x_vals See fittedPlot() function.
#' @param design_matrix_list See fittedPlot() function.
#' @param color See fittedPlot() function.
#' @param line.type See fittedPlot() function.
#' @param xlab See fittedPlot() function.
#'
#' @importFrom gamlss is.gamlss
#' @importFrom graphics par
#' @return Plot of fitted curves for one variable
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
#' fina_coef_list <- list(mu.coef = c(10, 1, 2), sigma.coef = c(4, 3),
#' nu.coef = 2, tau.coef = 3) # Declared based on series of last_updates
#' pooled_coefs <- dgamlss_pooled_coefs(current_coef_list, n_communications,
#' last_update$deviance)
#' summary_object <- dgamlss_create_summary(site1_gamlss, pooled_coefs,
#' pooled_hess, pooled_n)
#' dgamlss_fitted_plot(summary_object)
#' }
dgamlss_fitted_plot <- function(object,
                                x_vals,
                                design_matrix_list = NULL,
                                color = TRUE, line.type = FALSE, xlab = NULL) {
  family <- as.gamlss.family(object$family[1])
  if (!is.gamlss(object)) {
    stop(paste(
      "This is not an gamlss object", "\n",
      ""
    ))
  }

  if (is.null(design_matrix_list)) {
    stop(paste(
      "The design matrix list is not specified",
      "\n", ""
    ))
  } else if (!is.list((design_matrix_list))) {
    design_matrix_list <- list(design_matrix_list)
  }

  design_matrix_list <- lapply(design_matrix_list, function(item) as.matrix(item))

  nopar <- length(object$parameters)

  if (length(design_matrix_list) != nopar) {
    stop("Must provide design matrix for all parameters")
  }

  if (!all(unlist(lapply(design_matrix_list, function(item) (ncol(item) == 1) | (is.numeric(item)))))) {
    warning("All design matrices must be given in the same row order as the x_val vector.")
  }

  if (length(unique(unlist(lapply(design_matrix_list, function(item) nrow(item))))) != 1) {
    stop("All design matrices must have the same number of rows")
  }
  param <- object$parameters
  xvar <- if (is.null(xlab)) "x" else xlab
  x.o <- x_vals[order(x_vals)]
  let <- c("(a)", "(b)", "(c)", "(d)")
  if ("mu" %in% object$parameters) {
    vvv <- c(1, 1)
    mu.o <- family$mu.linkinv((design_matrix_list[[1]] %*% object$mu.coefficients))[order(x_vals)]
    if (all(abs(mu.o - mu.o[1]) < 1e-04)) {
      mu.o <- rep(mu.o[1], length(mu.o))
    }
    mat <- cbind(mu.o)
  }
  if ("sigma" %in% object$parameters) {
    vvv <- c(2, 1)
    sigma.o <- family$sigma.linkinv((design_matrix_list[[2]] %*% object$sigma.coefficients))[order(x_vals)]
    if (all(abs(sigma.o - sigma.o[1]) < 1e-04)) {
      sigma.o <- rep(sigma.o[1], length(sigma.o))
    }
    mat <- cbind(mat, sigma.o)
  }
  if ("nu" %in% object$parameters) {
    vvv <- c(3, 1)
    nu.o <- family$nu.linkinv((design_matrix_list[[3]] %*% object$nu.coefficients))[order(x_vals)]
    if (all(abs(nu.o - nu.o[1]) < 1e-04)) {
      nu.o <- rep(nu.o[1], length(nu.o))
    }
    mat <- cbind(mat, nu.o)
  }
  if ("tau" %in% object$parameters) {
    vvv <- c(2, 2)
    tau.o <- family$tau.linkinv((design_matrix_list[[4]] %*% object$tau.coefficients))[order(x_vals)]
    if (all(abs(tau.o - tau.o[1]) < 1e-04)) {
      tau.o <- rep(tau.o[1], length(tau.o))
    }
    mat <- cbind(mat, tau.o)
  }
  op <- par(
    mfrow = vvv, mar = par("mar") + c(
      0, 1, 0,
      0
    ), col.axis = "blue4", col.main = "blue4", col.lab = "blue4",
    cex = 0.5, cex.lab = 1.3, cex.axis = 1, cex.main = 1.3
  )
  for (ii in 1:nopar) {
    pp <- as.character(param[ii])
    if (color == TRUE) {
      col <- "darkgreen"
    } else {
      col <- "black"
    }
    plot(x.o, mat[, ii],
         xlab = xvar, ylab = pp, main = let[ii],
         col = col, col.axis = "mediumblue", type = "l",
         frame.plot = TRUE
    )
  }
  par(op)
}
