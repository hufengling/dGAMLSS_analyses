#' Simulation coordinating function for distributed GAMLSS model fitting
#'
#' @param mu.formula Formula for the mu component of the GAMLSS model.
#' @param sigma.formula Formula for the sigma component of the GAMLSS model.
#' @param nu.formula Formula for the nu component of the GAMLSS model.
#' @param tau.formula Formula for the tau component of the GAMLSS model.
#' @param family GAMLSS family object specifying the distribution and link functions.
#' @param weights Optional weights for likelihood analysis.
#' @param contrasts Optional contrasts for model parameters.
#' @param site_data List of data frames containing local datasets.
#' @param all_inits Optional named list containing initial values for all coefficients.
#' @param coef_crit Convergence criterion for coefficient changes.
#' @param local_site_data Data from coordinating site. Used to set up initial values for some GAMLSS distributions.
#' @param lambda_list Named list of fixed lambdas, if fixed penalty smooth terms are desired. Names must be "mu", "sigma", "nu", and "tau".
#' @param penalty_matrix_list Named list of penalty matrices, output from dgamlss_bs(). Names must be "mu", "sigma", "nu", and "tau".
#' @param k GAIC multiplier for automated penalty selection. AIC uses k = 2. BIC uses k = log(n). Higher k corresponds to stronger penalty.
#' @param is_orthogonal TRUE/FALSE indicator of whether the dgamlss_bs() output was orthogonalized or not.
#' @param verbose TRUE/FALSE for verbosity.
#' @param ... Additional arguments to be passed to the underlying GAMLSS fitting function.
#'
#' Simulates dgamlss fitting between a coordinating site and all local sites. Only summary statistics are passed from local sites. The coordinating site then aggregates local summary statistics and sends back pooled updates.
#'
#' @import gamlss
#' @return A list containing pooled coefficients, number of communications, number of reduced updates, and global deviance.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data(abdom)
#' site1_data <- abdom[1:110, ]
#' site2_data <- abdom[111:610, ]
#' site_data <- list(site1_data, site2_data)
#' site1_gamlss <- gamlss(mu.formula = y ~ x1 + x2, sigma.formula = ~ x1,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1_data)
#' result <- dgamlss_coordinating(mu.formula = y ~ x1 + x2,
#' sigma.formula = ~ x1, nu.formula = ~ 1, tau.formula = ~ 1,
#' family = BCT, site_data, all_inits = dgamlss_get_inits(site1_gamlss))
#' }
dgamlss_coordinating_penalized <- function(mu.formula,
                                 sigma.formula = "~ 1",
                                 nu.formula = "~ 1",
                                 tau.formula = "~ 1",
                                 family = NO(),
                                 weights = NULL, # for weighted likelihood analysis
                                 contrasts = NULL, # one type of contrasts for all parameters
                                 local_site_data,
                                 site_data,
                                 all_inits = NULL,
                                 coef_crit = 0.01,
                                 lambda_list = NULL,
                                 penalty_matrix_list = NULL,
                                 k = 2,
                                 is_orthogonal = FALSE,
                                 verbose = FALSE,
                                 ...) {
  # Return one site ==========
  return_one_site <- function(site, update,
                              get_ssr = FALSE,
                              proposed_coef_list = NULL) {
    dgamlss_RS(formula = formula(mu.formula),
               sigma.formula = formula(sigma.formula),
               nu.formula = formula(nu.formula),
               tau.formula = formula(tau.formula),
               family = family,
               weights = weights,
               data = site,
               dgamlss.update = update,
               get_ssr = get_ssr,
               proposed_coef_list = proposed_coef_list,
               ...)
  }

  outer_iter <- 0
  n_communications <- 0
  n_reduced <- -1

  if (verbose) {
    print("Running local GAMLSS to initialize coefficients")
  }

  local_gamlss <- gamlss(formula = formula(mu.formula),
                         sigma.formula = formula(sigma.formula),
                         nu.formula = formula(nu.formula),
                         tau.formula = formula(tau.formula),
                         family = family,
                         weights = weights,
                         data = local_site_data,
                         control = gamlss.control(trace = verbose),
                         ...)
  # Initialize coefficients ====================================================
  if (!is.null(all_inits)) {
    mu_coef_init <- all_inits$mu
    sigma_coef_init <- all_inits$sigma
    nu_coef_init <- all_inits$nu
    tau_coef_init <- all_inits$tau
  } else {
    if (!is_orthogonal) {
      repeat_init <- lapply(penalty_matrix_list$smooth_index_list, sum) # If basis is not orthogonal, every basis element contributes to Intercept term
    } else {
      repeat_init <- lapply(penalty_matrix_list$smooth_index_list, \(x) {1}) # If basis is orthogonal, initial value only pertains to Intercept
    }

    gamlss_family <- as.gamlss.family(family)
    mu_coef_init <- c(rep(gamlss_family$mu.linkfun(mean(local_gamlss$y)), repeat_init$mu),
                      rep(0, length(local_gamlss$mu.coefficients) - repeat_init$mu))
    y <- local_gamlss$y
    if ("sigma" %in% names(gamlss_family$parameters)) {
      sigma_coef_init <- c(rep(gamlss_family$sigma.linkfun(eval(gamlss_family$sigma.initial))[1], repeat_init$sigma),
                           rep(0, length(local_gamlss$sigma.coefficients) - repeat_init$sigma))
    } else {
      sigma_coef_init <- NULL
    }
    if ("nu" %in% names(gamlss_family$parameters)) {
      nu_coef_init <- c(rep(gamlss_family$nu.linkfun(eval(gamlss_family$nu.initial))[1], repeat_init$nu),
                        rep(0, length(local_gamlss$nu.coefficients) - repeat_init$nu))
    } else {
      nu_coef_init <- NULL
    }
    if ("tau" %in% names(gamlss_family$parameters)) {
      tau_coef_init <- c(rep(gamlss_family$tau.linkfun(eval(gamlss_family$tau.initial))[1], repeat_init$tau),
                         rep(0, length(local_gamlss$tau.coefficients) - repeat_init$tau))
    } else {
      tau_coef_init <- NULL
    }
  }

  mu_inner <- mu_coef_init
  sigma_inner <- sigma_coef_init
  nu_inner <- nu_coef_init
  tau_inner <- tau_coef_init

  new_update <- dgamlss_send_coefs("mu",
                                   mu_inner,
                                   sigma_inner,
                                   nu_inner,
                                   tau_inner)
  if (verbose) {
    print("Initial coefficients below:")
    print(new_update)
  }

  old_update <- new_update
  old_outer_update <- new_update

  # Get XTX matrix =================================================
  n_communications <- n_communications + 1
  site_xtx_list <- lapply(site_data, function(site) {
    dgamlss_get_inference(gamlss_mock_fit(formula = formula(mu.formula),
                                          sigma.formula = formula(sigma.formula),
                                          nu.formula = formula(nu.formula),
                                          tau.formula = formula(tau.formula),
                                          family = BCPE(),
                                          data = site),
                          penalized_parameters = names(penalty_matrix_list[-length(penalty_matrix_list)]))
  })

  xtx_list <- list(xtx_mu = lapply(site_xtx_list, function(site) site$xtx_list$mu),
                   xtx_sigma = lapply(site_xtx_list, function(site) site$xtx_list$sigma),
                   xtx_nu = lapply(site_xtx_list, function(site) site$xtx_list$nu),
                   xtx_tau = lapply(site_xtx_list, function(site) site$xtx_list$tau))
  n <- sum(sapply(site_xtx_list, function(site) site$n))

  # Outer iteration ===========================================================
  while(TRUE) {
    outer_iter <- outer_iter + 1
    if (verbose) {
      print(paste("Outer", outer_iter))
    }
    # Update mu ================================================================================
    n_reduced <- n_reduced + 1
    is_updated <- FALSE
    new_update$to_update <- "mu"
    while(!dgamlss_check_convergence(new_update,
                                     old_update,
                                     coef_crit = coef_crit) | !is_updated) {
      if (verbose) {
        print("Updating mu")
      }

      n_communications <- n_communications + 1
      site_info <- lapply(site_data, return_one_site, update = new_update)
      if ("mu" %in% names(penalty_matrix_list)) {
        proposed_mu <- dgamlss_aggregate_coef(site_info,
                                         penalty_matrix = penalty_matrix_list$mu)
        n_communications <- n_communications + 1
        site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                get_ssr = TRUE, proposed_coef_list = proposed_mu$proposed_coef_list)
        new_mu <- dgamlss_select_lambda(proposed_mu$lambda_vec, proposed_mu$proposed_coef_list,
                                        site_ssr_list, xtx_list$xtx_mu, penalty_matrix_list$mu,
                                        k = k)
        mu_inner <- proposed_mu$proposed_coef_list[[new_mu$index]]
      } else {
        new_mu <- dgamlss_aggregate_coef(site_info)
        mu_inner <- new_mu$coef
      }

      if (verbose) {
        print(mu_inner)
        print(new_mu$deviance)
      }

      old_update <- new_update
      new_update <- dgamlss_send_coefs("mu",
                                       mu_inner,
                                       sigma_inner,
                                       nu_inner,
                                       tau_inner)
      is_updated <- TRUE
    }
    global_deviance <- new_mu$deviance

    if (dgamlss_check_convergence(new_update,
                                  old_outer_update,
                                  coef_crit = coef_crit) & (outer_iter > 1)) {
      break
    }

    # Update sigma ============================================================================================
    if (!is.null(sigma_coef_init)) {
      n_reduced <- n_reduced + 1
      is_updated <- FALSE
      new_update$to_update <- "sigma"
      while(!dgamlss_check_convergence(new_update,
                                       old_update,
                                       coef_crit = coef_crit) | !is_updated) {
        if (verbose) {
          print("Updating sigma")
        }

        n_communications <- n_communications + 1
        site_info <- lapply(site_data, return_one_site, update = new_update)
        if ("sigma" %in% names(penalty_matrix_list)) {
          proposed_sigma <- dgamlss_aggregate_coef(site_info,
                                                penalty_matrix = penalty_matrix_list$sigma)
          n_communications <- n_communications + 1
          site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                  get_ssr = TRUE, proposed_coef_list = proposed_sigma$proposed_coef_list)
          new_sigma <- dgamlss_select_lambda(proposed_sigma$lambda_vec, proposed_sigma$proposed_coef_list,
                                          site_ssr_list, xtx_list$xtx_sigma, penalty_matrix_list$sigma,
                                          k = k)
          sigma_inner <- proposed_sigma$proposed_coef_list[[new_sigma$index]]
        } else {
          new_sigma <- dgamlss_aggregate_coef(site_info)
          sigma_inner <- new_sigma$coef
        }
        if (verbose) {
          print(sigma_inner)
          print(new_sigma$deviance)
        }

        old_update <- new_update
        new_update <- dgamlss_send_coefs("sigma",
                                         mu_inner,
                                         sigma_inner,
                                         nu_inner,
                                         tau_inner)
        is_updated <- TRUE
      }
      global_deviance <- new_sigma$deviance
    }

    # Update nu ==============================================================================================================
    if (!is.null(nu_coef_init)) {
      n_reduced <- n_reduced + 1
      is_updated <- FALSE
      new_update$to_update <- "nu"
      while(!dgamlss_check_convergence(new_update,
                                       old_update,
                                       coef_crit = coef_crit) | !is_updated) {
        if (verbose) {
          print("Updating nu")
        }

        n_communications <- n_communications + 1
        site_info <- lapply(site_data, return_one_site, update = new_update)
        if ("nu" %in% names(penalty_matrix_list)) {
          proposed_nu <- dgamlss_aggregate_coef(site_info,
                                                penalty_matrix = penalty_matrix_list$nu)
          n_communications <- n_communications + 1
          site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                  get_ssr = TRUE, proposed_coef_list = proposed_nu$proposed_coef_list)
          new_nu <- dgamlss_select_lambda(proposed_nu$lambda_vec, proposed_nu$proposed_coef_list,
                                          site_ssr_list, xtx_list$xtx_nu, penalty_matrix_list$nu,
                                          k = k)
          nu_inner <- proposed_nu$proposed_coef_list[[new_nu$index]]
        } else {
          new_nu <- dgamlss_aggregate_coef(site_info)
          nu_inner <- new_nu$coef
        }
        if (verbose) {
          print(nu_inner)
          print(new_nu$deviance)
        }

        old_update <- new_update
        new_update <- dgamlss_send_coefs("nu",
                                         mu_inner,
                                         sigma_inner,
                                         nu_inner,
                                         tau_inner)
        is_updated <- TRUE
      }
      global_deviance <- new_nu$deviance
    }

    # Update tau ====================================================================================================
    if (!is.null(tau_coef_init)) {
      n_reduced <- n_reduced + 1
      is_updated <- FALSE
      new_update$to_update <- "tau"
      while(!dgamlss_check_convergence(new_update,
                                       old_update,
                                       coef_crit = coef_crit) | !is_updated) {
        if (verbose) {
          print("Updating tau")
        }

        n_communications <- n_communications + 1
        site_info <- lapply(site_data, return_one_site, update = new_update)
        if ("tau" %in% names(penalty_matrix_list)) {
          proposed_tau <- dgamlss_aggregate_coef(site_info,
                                                penalty_matrix = penalty_matrix_list$tau)
          n_communications <- n_communications + 1
          site_ssr_list <- lapply(site_data, return_one_site, update = new_update,
                                  get_ssr = TRUE, proposed_coef_list = proposed_tau$proposed_coef_list)
          new_tau <- dgamlss_select_lambda(proposed_tau$lambda_vec, proposed_tau$proposed_coef_list,
                                          site_ssr_list, xtx_list$xtx_tau, penalty_matrix_list$tau,
                                          k = k)
          tau_inner <- proposed_tau$proposed_coef_list[[new_tau$index]]
        } else {
          new_tau <- dgamlss_aggregate_coef(site_info)
          tau_inner <- new_tau$coef
        }
        if (verbose) {
          print(tau_inner)
          print(new_tau$deviance)
        }

        old_update <- new_update
        new_update <- dgamlss_send_coefs("tau",
                                         mu_inner,
                                         sigma_inner,
                                         nu_inner,
                                         tau_inner)
        is_updated <- TRUE
      }
      global_deviance <- new_tau$deviance
    }

    old_outer_update <- new_update
  }

  return(list(mu_coefs = new_update$mu.coef,
              sigma_coefs = new_update$sigma.coef,
              nu_coefs = new_update$nu.coef,
              tau_coefs = new_update$tau.coef,
              n_communications = n_communications,
              n_reduced = n_reduced,
              global_deviance = global_deviance,
              edf_vec = c(new_mu$edf_opt, new_sigma$edf_opt,
                          new_nu$edf_opt, new_tau$edf_opt)))
}
