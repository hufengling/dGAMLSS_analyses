#' Aggregate inference across sites
#'
#' @param site_list List of joined outputs from dgamlss_aggregate_inference
#' @param pooled_coefs Named list of final pooled coefficients. Names must be "mu", "sigma", "nu", "tau".
#' @param lambda_list Named list of fixed lambdas, if fixed penalty smooth terms are desired. Names must be "mu", "sigma", "nu", and "tau".
#' @param penalty_matrix_list Named list of penalty matrices, output from dgamlss_bs(). Names must be "mu", "sigma", "nu", and "tau".
#'
#' @return List containing pooled Hessian, effective degrees of freedom, and overall number of samples.
#' @export
#'
#' @examples
#' \dontrun{
#' dgamlss_aggregate_inference(site_list, pooled_coef, lambda_list, penalty_matrix_list)
#' }
dgamlss_aggregate_inference <- function(site_list,
                                        pooled_coefs,
                                        lambda_list = NULL,
                                        penalty_matrix_list = NULL) {
  n_sites <- length(site_list)
  hessian_list <- vector("list", n_sites)
  xtx_list <- list(mu = vector("list", n_sites),
                   sigma = vector("list", n_sites),
                   nu = vector("list", n_sites),
                   tau = vector("list", n_sites))
  # rss_list <- list(mu = vector("list", n_sites),
  #                  sigma = vector("list", n_sites),
  #                  nu = vector("list", n_sites),
  #                  tau = vector("list", n_sites))
  n_vector <- rep(0, n_sites)
  edf_list <- list()
  # gcv_list <- list()

  for (i in 1:n_sites) {
    hessian_list[[i]] <- site_list[[i]]$hessian

    xtx_list$mu[[i]] <- site_list[[i]]$xtx$mu
    xtx_list$sigma[[i]] <- site_list[[i]]$xtx$sigma
    xtx_list$nu[[i]] <- site_list[[i]]$xtx$nu
    xtx_list$tau[[i]] <- site_list[[i]]$xtx$tau

    # rss_list$mu[[i]] <- site_list[[i]]$rss$mu
    # rss_list$sigma[[i]] <- site_list[[i]]$rss$sigma
    # rss_list$nu[[i]] <- site_list[[i]]$rss$nu
    # rss_list$tau[[i]] <- site_list[[i]]$rss$tau

    n_vector[i] <- site_list[[i]]$n
  }

  pooled_hessian <- Reduce("+", hessian_list)
  pooled_n <- sum(n_vector)

  if (!is.null(penalty_matrix_list)) {
    penalized_parameters <- names(penalty_matrix_list)
  } else {
    penalized_parameters <- NULL
  }
  if ("mu" %in% penalized_parameters) {
    xtxinvxtx <- solve(Reduce("+", xtx_list$mu) + lambda_list$mu * penalty_matrix_list$mu) %*% Reduce("+", xtx_list$mu)
    edf_list$mu <- sum(diag(xtxinvxtx))#2 * sum(diag(xtxinvxtx)) - sum(diag(xtxinvxtx %*% t(xtxinvxtx)))

    # gcv_list$mu <- Reduce("+", rss_list$mu) * pooled_n / (pooled_n - sum(diag(xtxinvxtx)))^2
  } else {
    edf_list$mu <- length(pooled_coefs$mu)
    # gcv_list$mu <- NULL
  }

  if ("sigma" %in% penalized_parameters) {
    xtxinvxtx <- solve(Reduce("+", xtx_list$sigma) + lambda_list$sigma * penalty_matrix_list$sigma) %*% Reduce("+", xtx_list$sigma)
    edf_list$sigma <- 2 * sum(diag(xtxinvxtx)) - sum(diag(xtxinvxtx %*% t(xtxinvxtx)))

    # gcv_list$sigma <- Reduce("+", rss_list$sigma) * pooled_n / (pooled_n - sum(diag(xtxinvxtx)))^2
  } else {
    edf_list$sigma <- length(pooled_coefs$sigma)
    # gcv_list$sigma <- NULL
  }

  if ("nu" %in% penalized_parameters) {
    xtxinvxtx <- solve(Reduce("+", xtx_list$nu) + lambda_list$nu * penalty_matrix_list$nu) %*% Reduce("+", xtx_list$nu)
    edf_list$nu <- 2 * sum(diag(xtxinvxtx)) - sum(diag(xtxinvxtx %*% t(xtxinvxtx)))

    # gcv_list$nu <- Reduce("+", rss_list$nu) * pooled_n / (pooled_n - sum(diag(xtxinvxtx)))^2
  } else {
    edf_list$nu <- length(pooled_coefs$nu)
    # gcv_list$nu <- NULL
  }

  if ("tau" %in% penalized_parameters) {
    xtxinvxtx <- solve(Reduce("+", xtx_list$tau) + lambda_list$tau * penalty_matrix_list$tau) %*% Reduce("+", xtx_list$tau)
    edf_list$tau <- 2 * sum(diag(xtxinvxtx)) - sum(diag(xtxinvxtx %*% t(xtxinvxtx)))

    # gcv_list$tau <- Reduce("+", rss_list$tau) * pooled_n / (pooled_n - sum(diag(xtxinvxtx)))^2
  } else {
    edf_list$tau <- length(pooled_coefs$tau)
    # gcv_list$tau <- NULL
  }

  return(list(pooled_hessian = pooled_hessian,
              edf_list = edf_list,
              # gcv_list = gcv_list,
              pooled_n = pooled_n))
}
