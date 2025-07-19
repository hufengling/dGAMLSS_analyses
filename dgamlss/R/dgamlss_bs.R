#' Generate B-spline basis matrix
#'
#' @param x The vector of input values from the local site.
#' @param x_min The pooled minimum value of x across all sites.
#' @param x_max The pooled maximum value of x across all sites.
#' @param spline_prefix Optional name prefix for the spline basis columns. Default "spline_basis_"
#' @param knot_positions Vector of knot positions, agreed upon across all sites. We recommend choosing knot positions using some approximation of pooled quantiles.
#' @param n_knots Alternative to specifying knot positions. If n_knots is provided, knots will be placed at identical intervals.
#' @param degree Degree of the B-spline basis functions. Default 3 for cubic B-splines.
#' @param order Order of B-spline penalty. Default 2 for cubic B-splines.
#' @param orthogonalize Whether to get spline representation where intercept is orthogonal to the non-linear representations. Can only be used with fixed effect splines or fixed penalty splines. Cannot be used with empirically penalized splines.
#'
#' @importFrom splines splineDesign
#' @return A design matrix containing the B-spline basis functions evaluated at the input values.
#' @export
#'
#' @examples
#' # Generate example data
#' x <- 1:100
#'
#' # Generate B-spline basis matrix
#' bs_matrix <- dgamlss_bs(x, x_min = 0, x_max = 100, spline_prefix = "my_spline",
#' n_knots = 10, degree = 3)
dgamlss_bs <- function(x, x_min, x_max, spline_prefix = NULL,
                       knot_positions = NULL, n_knots = NULL,
                       degree = 3, order = 2,
                       orthogonalize = TRUE) {
  bs_start <- x_min - 0.01 * (x_max - x_min)
  bs_end <- x_max + 0.01 * (x_max - x_min)

  if ((!is.null(knot_positions) & !is.null(n_knots)) | (is.null(knot_positions) & is.null(n_knots))){
    stop("Must choose one of knot_position or n_knots. Other must be set to NULL")
  }

  if (!is.null(knot_positions)) {
    n_lower <- sum(knot_positions < x_min)
    n_upper <- sum(knot_positions > x_max)
    if (n_lower < degree | n_upper < degree) {
      stop(paste("Must have at least", degree,
                 "knot_positions extending past x_min and",
                 degree, "knot_positions extending past x_max."))
    }
  }

  if (!is.null(n_knots)) {
    n_knots <- n_knots - 1 # Will place "lost" inner knot at the start
    interval <- (bs_end - bs_start) / n_knots
    knot_positions <- seq(bs_start - degree * interval, bs_end + degree * interval, by = interval)
  }

  B <- splineDesign(knot_positions, x = x, outer.ok = TRUE, ord = degree + 1)
  attr(B, "true_knots") <- knot_positions[-c(1:(degree - 1), (length(knot_positions) - (degree - 2)):length(knot_positions))]
  attr(B, "all_knots") <- knot_positions

  r <- ncol(B)
  D <- if (order == 0) diag(r) else diff(diag(r), diff = order)

  if (orthogonalize) {
    x_sim <- seq(x_min, x_max, length.out = 1000) # Generate x data from min to max
    B_sim <- splineDesign(knot_positions, x = x_sim, outer.ok = TRUE, ord = degree + 1) # Generate spline design from simulated x (need x along full range of spline basis)
    B_sim[, 1] <- rep(1, length(x_sim)) # Replace first column with intercept -- rest of basis should be orthogonal to intercept
    R <- qr.R(qr(B_sim)) # Get R from QR decomposition

    B[, 1] <- rep(1, nrow(B))
    B <- B %*% solve(R)
    #browser()
    B[, 1] <- rep(1, nrow(B))
    #B <- B[, -1]

    D <- D %*% solve(R)
    #D <- D[, -1]
    attr(B, "R") <- R
  }

  attr(B, "D") <- D

  attr(B, "P") <- t(D) %*% D

  if (is.null(spline_prefix)) {
    colnames(B) <- paste0("spline_basis_", 1:ncol(B))
  } else {
    if (substr(spline_prefix, nchar(spline_prefix), nchar(spline_prefix)) != "_") {
      spline_prefix <- paste0(spline_prefix, "_")
    }
    colnames(B) <- paste0(spline_prefix, 1:ncol(B))
  }

  #attr(B, "class") <- "smooth"

  return(B)
}
