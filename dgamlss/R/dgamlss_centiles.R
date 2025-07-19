#' Plot GAMLSS centiles
#'
#' @param obj See GAMLSS centiles() function.
#' @param x_val See GAMLSS centiles() function. If centile estimates are jagged, this may be due to inclusion of categorical variables in the GAMLSS models. Jaggedness can be corrected by replacing obj$mu.fv, obj$sigma.fv, obj$nu.fv, obj$tau.fv values with predictions under a single category.
#' @param design_matrix_list See GAMLSS centiles() function.
#' @param cent See GAMLSS centiles() function.
#' @param legend See GAMLSS centiles() function.
#' @param ylab See GAMLSS centiles() function.
#' @param xlab See GAMLSS centiles() function.
#' @param main See GAMLSS centiles() function.
#' @param main.gsub See GAMLSS centiles() function.
#' @param xleg See GAMLSS centiles() function.
#' @param yleg See GAMLSS centiles() function.
#' @param xlim See GAMLSS centiles() function.
#' @param ylim See GAMLSS centiles() function.
#' @param save See GAMLSS centiles() function.
#' @param plot See GAMLSS centiles() function.
#' @param points See GAMLSS centiles() function.
#' @param pch See GAMLSS centiles() function.
#' @param cex See GAMLSS centiles() function.
#' @param col See GAMLSS centiles() function.
#' @param col.centiles See GAMLSS centiles() function.
#' @param lty.centiles See GAMLSS centiles() function.
#' @param lwd.centiles See GAMLSS centiles() function.
#' @param ... Additional parameters.
#'
#' @return Creates centile plot.
#' @export
#'
#' @importFrom grDevices gray
#' @importFrom graphics lines title
#'
#' @examples
#' \dontrun{
#' dgamlss_centiles(dgamlss_summary_object, x_val = df$age)
#' }
dgamlss_centiles <- function(obj,
                             x_val,
                             design_matrix_list,
                             cent = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6),
                             legend = TRUE,
                             ylab = "y",
                             xlab = "x",
                             main = NULL,
                             main.gsub = "@",
                             xleg = min(x_val), # Note main and a substitution marker (an experiment!)
                             yleg = NULL,
                             xlim = range(x_val),
                             ylim = NULL,
                             save = FALSE,
                             plot = TRUE,
                             points = TRUE, # this is new Saturday, August 14, 2010 MS
                             pch = 15,
                             cex = 0.5,
                             col = gray(0.7),
                             col.centiles = 1:length(cent) + 2,
                             lty.centiles = 1,
                             lwd.centiles = 1, # Handling for line appearance
                             ...) {
  if (!is.gamlss(obj)) stop(paste("This is not an gamlss object", "\n", ""))

  fname <- obj$family[1]
  family <- as.gamlss.family(fname)
  qfun <- paste("q", fname, sep = "")
  Title <- paste("Centile curves using", fname, sep = " ") # Note change to default title handling
  main <- if (is.null(main)) {
    paste("Centile curves using", fname, sep = " ")
  } else {
    gsub(main.gsub, Title, main)
  }

  yvar <- design_matrix_list[[1]] %*% obj$mu.coefficients

  if (is.null(ylim)) {
    ylim <- range(yvar) * c(0.8, 1.25)
  }

  if (is.null(yleg)) {
    yleg <- ylim[2]
  }

  if (plot) {
    lty.centiles <- rep(lty.centiles, length(cent))
    lwd.centiles <- rep(lwd.centiles, length(cent))
    col.centiles <- rep(col.centiles, length(cent))
    if (points == TRUE) {
      plot(x_val, yvar,
           type = "p", col = col, pch = pch, cex = cex,
           xlab = xlab, ylab = ylab, xlim = xlim, ylim, ...
      )
    } else {
      plot(x_val, yvar,
           type = "n", col = col, pch = pch,
           xlab = xlab, ylab = ylab, xlim = xlim, ylim, ...
      )
    }
    title(main)
  }
  col <- 3 # set this to 1 if you do not want colour
  lpar <- length(obj$parameters)
  ii <- 0
  per <- rep(0, length(cent))
  # per <- cent/100
  for (var in cent)
  {
    if (lpar == 1) {
      newcall <- call(qfun, var / 100,
                      mu = family$mu.linkinv(design_matrix_list[[1]] %*% obj$mu.coefficients)[order(x_val)]
      )
    } else if (lpar == 2) {
      newcall <- call(qfun, var / 100,
                      mu = family$mu.linkinv(design_matrix_list[[1]] %*% obj$mu.coefficients)[order(x_val)],
                      sigma = family$sigma.linkinv(design_matrix_list[[2]] %*% obj$sigma.coefficients)[order(x_val)]
      )
    } else if (lpar == 3) {
      newcall <- call(qfun, var / 100,
                      mu = family$mu.linkinv(design_matrix_list[[1]] %*% obj$mu.coefficients)[order(x_val)],
                      sigma = family$sigma.linkinv(design_matrix_list[[2]] %*% obj$sigma.coefficients)[order(x_val)],
                      nu = family$nu.linkinv(design_matrix_list[[3]] %*% obj$nu.coefficients)[order(x_val)]
      )
    } else {
      newcall <- call(qfun, var / 100,
                      mu = family$mu.linkinv(design_matrix_list[[1]] %*% obj$mu.coefficients)[order(x_val)],
                      sigma = family$sigma.linkinv(design_matrix_list[[2]] %*% obj$sigma.coefficients)[order(x_val)],
                      nu = family$nu.linkinv(design_matrix_list[[3]] %*% obj$nu.coefficients)[order(x_val)],
                      tau = family$tau.linkinv(design_matrix_list[[4]] %*% obj$tau.coefficients)[order(x_val)]
      )
    }
    ii <- ii + 1
    ll <- eval(newcall)
    if (plot) {
      lines(x_val, ll, col = col.centiles[ii], lty = lty.centiles[ii], lwd = lwd.centiles[ii], ...)
    }
    #per[ii] <- (1 - sum(yvar > ll) / length(yvar)) * 100
    if (!save) cat("% of cases below ", var, "centile is ", per[ii], "\n")
  }
  if (plot) { # Legend moved outside plotting loop
    if (legend == TRUE) {
      legend(list(x = xleg, y = yleg),
             legend = cent,
             col = col.centiles, lty = lty.centiles, lwd = lwd.centiles,
             ncol = 1, ...
      )
    }
  }
  if (save) {
    return(cbind(cent, per))
  }
}
