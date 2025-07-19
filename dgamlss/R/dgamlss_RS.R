#' dgamlss analog of the gamlss() function
#'
#' @param formula mu formula
#' @param sigma.formula sigma formula
#' @param nu.formula nu formula
#' @param tau.formula tau formula
#' @param family family from gamlss.dist
#' @param data local data
#' @param weights weights
#' @param contrasts See gamlss() function.
#' @param dgamlss.update Current update from dgamlss_send_coefs()
#' @param control See gamlss() function.
#' @param i.control See gamlss() function.
#' @param get_ssr Get sum of squared residuals for empirical penalty selection. Only useful for penalized smooth term fitting, after getting a set of proposed coefficients under various penalties. Default FALSE.
#' @param proposed_coef_list List of proposed coefficients across a grid of lambda penalties, sent from coordinating site. Use with get_ssr = TRUE.
#'
#' @import gamlss
#' @import gamlss.dist
#' @import gamlss.data
#' @import survival
#' @importFrom methods is
#' @return Summary statistics for the site that should be sent to the coordinating site. Coordinating site can aggregate dgamlss_RS outputs using dgamlss_aggregate_coef()
#' @export
#'
#' @examples
#' \dontrun{
#' # To update the "mu" coefficients
#' data(abdom)
#' site1_data <- abdom[1:110, ]
#' site2_data <- abdom[111:610, ]
#' current_update <- dgamlss_send_coefs(to_update = "mu", mu.coef = c(10, 1, 2),
#' sigma.coef = c(4, 3), nu.coef = 2, tau.coef = 1)
#' site1_update <- dgamlss_RS(formula = y ~ x1 + x2, sigma.formula = ~ x,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1_data,
#' dgamlss.update = current_update)
#' site2_update <- dgamlss_RS(formula = y ~ x1 + x2, sigma.formula = ~ x,
#' nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site2_data,
#' dgamlss.update = current_update)
#' pooled_update <- dgamlss_aggregate_coef(list(site1_update, site2_update))
#' }
dgamlss_RS <- function(formula = formula(data),
                       sigma.formula = ~1,
                       nu.formula = ~1,
                       tau.formula = ~1,
                       family = NO(),
                       data,
                       weights = NULL, # for weighted likelihood analysis
                       contrasts = NULL, # one type of contrasts for all parameters
                       dgamlss.update = dgamlss.update(),
                       get_ssr = FALSE,
                       proposed_coef_list = NULL,
                       control = gamlss.control(),
                       i.control = glim.control()) {# the inner circle control (GLIM)
  ## =========================================================================================
  ## -----------------------------------------------------------------------------------------
  ## this function is used to extract the formula for the parameters other than mu
  ## -----------------------------------------------------------------------------------------
  ## =========================================================================================
  other.formula <- function(form) {
    dform <- formula(form)
    if (length(dform) == 2) {
      dform[3] <- dform[2] # taking 1 in position [3]
      dform[2] <- if (is(formula, "terms")) formula[[2]] else formula[2] # ms 31-12-08   # put y in position 2
    }
    dform
  }

  ## =========================================================================================
  ## -----------------------------------------------------------------------------------------
  ## this function is getting the smoothers at each parameter
  ## -----------------------------------------------------------------------------------------
  ## =========================================================================================
  get.smoothers <- function(term) {
    a <- attributes(term) #
    smoothers <- a$specials # S convert variable pointers to term pointers
    if (length(smoothers) > 0) {
      smoothers <- smoothers[sapply(smoothers, length) > 0]
      # smoothersR <-smoothers
      for (i in seq(along = smoothers))
      {
        tt <- smoothers[[i]]
        ff <- apply(a$factors[tt, , drop = FALSE], 2, any)
        smoothers[[i]] <- if (any(ff)) {
          seq(along = ff)[a$order == 1 & ff]
        } else {
          NULL
        }
      }
    }
    smoothers
  }

  ## get.object function =========================================================================================
  get.object <- function(what) {
    link <- eval(parse(text = (paste("family$", what, ".link", sep = ""))))
    linkfun <- eval(parse(text = (paste("family$", what, ".linkfun", sep = ""))))
    linkinv <- eval(parse(text = (paste("family$", what, ".linkinv", sep = ""))))
    dr <- eval(parse(text = (paste("family$", what, ".dr", sep = ""))))
    dldp <- switch(what,
                   "mu" = family$dldm,
                   "sigma" = family$dldd,
                   "nu" = family$dldv,
                   "tau" = family$dldt
    )
    d2ldp2 <- switch(what,
                     "mu" = family$d2ldm2,
                     "sigma" = family$d2ldd2,
                     "nu" = family$d2ldv2,
                     "tau" = family$d2ldt2
    )
    G.di <- family$G.dev.incr
    valid <- eval(parse(text = (paste("family$", what, ".valid", sep = ""))))
    object <- list(
      link = link, linkfun = linkfun, linkinv = linkinv, dr = dr,
      dldp = dldp, d2ldp2 = d2ldp2, G.di = G.di, valid = valid
    )
    if (length(eval(parse(text = (paste(what, ".smoothers", sep = ""))))) > 0) { # only if smoothing
      parAttrTermlevels <- eval(parse(text = (paste(what, ".a$term.labels", sep = ""))))
      boo <- unlist(eval(parse(text = (paste(what, ".smoothers", sep = "")))))
      who <- parAttrTermlevels[boo[order(boo)]]
      smooth.frame <- eval(parse(text = (paste(what, ".frame", sep = ""))))
      s <- matrix(0, N, length(who))
      dimnames(s) <- list(names(y), who)
      object$smooth <- s
      object$who <- who
      object$smooth.frame <- smooth.frame
    }
    object
  }

  RS_local <- function(no.warn = TRUE) {
    glim.fit <- function(family, X, y,
                         weights, fitted, offset,
                         step = 1, control = glim.control(),
                         auto, gd.tol) {
      eta <- family$linkfun(fitted)

      dldp <- family$dldp(fitted) # u score
      d2ldp2 <- family$d2ldp2(fitted) # second derivative of log-Likelihood
      d2ldp2 <- ifelse(d2ldp2 < -1e-15, d2ldp2, -1e-15)
      dr <- 1 / family$dr(eta) # deta/dmu

      w_matrix <- -(d2ldp2 / (dr * dr))
      w_matrix <- ifelse(w_matrix > 1e+10, 1e+10, w_matrix) # we need to stop the weights to go to Infty
      w_matrix <- ifelse(w_matrix < 1e-10, 1e-10, w_matrix)

      working_values <- (eta - offset) + dldp / (dr * w_matrix) # eta

      if (any(is.na(w_matrix)) || any(is.na(working_values))) {
        stop("NA's in the working vector or weights for parameter ", names(formals(family$valid)), "\n")
      }
      if (any(!is.finite(w_matrix)) || any(!is.finite(working_values))) {
        stop("Inf values in the working vector or weights for parameter ", names(formals(family$valid)), "\n")
      }

      if (get_ssr) {
        ssr_vec <- unlist(lapply(proposed_coef_list, \(coefs) {
          new_fitted <- family$linkinv(X %*% coefs + offset)
          # return(sum(diag(w_matrix) * weights * family$G.di(new_fitted)))
          return(sum(diag(w_matrix) * weights * (working_values - family$linkfun(new_fitted))^2))
        }))

        return(ssr_vec)
      }

      individual_deviances <- family$G.di(fitted) # deviance increment
      new_overall_deviance <- sum(weights * individual_deviances) # the global deviance

      site_data <- list(
        xtwx = t(X) %*% diag(as.vector(w_matrix * weights)) %*% X,
        xtwy = t(X) %*% diag(as.vector(w_matrix * weights)) %*% working_values,
        deviance = new_overall_deviance
      )

      # if (get_xtx) {
      #   site_data <- list(
      #     xtwx = t(X) %*% diag(as.vector(w_matrix * weights)) %*% X,
      #     xtwy = t(X) %*% diag(as.vector(w_matrix * weights)) %*% working_values,
      #     xtx = t(X) %*% X,
      #     deviance = new_overall_deviance
      #   )
      # }
      return(site_data)
    }

    ## the outer iteration starts here ---------------------------------------------------------------------------
    # the mean submodel
    if (dgamlss.update$to_update == "mu") {
      mu.fit <<- glim.fit(
        family = mu.object,
        X = mu.X, y = y,
        weights = w, fitted = mu, offset = mu.offset,
        step = 1, control = i.control,
        gd.tol = NULL, auto = FALSE
      )
      return(mu.fit)
    }

    # the scale-dispersion submodel
    if (dgamlss.update$to_update == "sigma") {
      sigma.fit <<- glim.fit(
        family = sigma.object,
        X = sigma.X, y = y,
        weights = w, fitted = sigma, offset = sigma.offset,
        step = 1, control = i.control,
        gd.tol = NULL, auto = FALSE
      )
      return(sigma.fit)
    }

    # the nu submodel
    if (dgamlss.update$to_update == "nu") {
      nu.fit <<- glim.fit(
        family = nu.object,
        X = nu.X, y = y,
        weights = w, fitted = nu, offset = nu.offset,
        step = 1, control = i.control,
        gd.tol = NULL, auto = FALSE
      )
      return(nu.fit)
    }

    # the tau submodel
    if (dgamlss.update$to_update == "tau") {
      tau.fit <<- glim.fit(
        family = tau.object,
        X = tau.X, y = y,
        weights = w, fitted = tau, offset = tau.offset,
        step = 1, control = i.control,
        gd.tol = NULL, auto = FALSE
      )
      return(tau.fit)
    }
  }


  ## Set up model matrix =========================================================================================
  ##       Save call for future reference
  gamlsscall <- match.call() #   the function call
  ## checking for NA in the data
  if (!missing(data) & any(is.na(data))) {
    stop("The data contains NA's, use data = na.omit(mydata)")
  }
  ##       Evaluate the model frame
  mnames <- c("", "formula", "data", "weights") #  "subset"  "na.action"
  cnames <- names(gamlsscall) # get the names of the arguments of the call
  cnames <- cnames[match(mnames, cnames, 0)] # keep only the ones that match with mnames
  mcall <- gamlsscall[cnames] # get in mcall all the relevant information but remember
  # that the first elenent will be NULL
  mcall[[1]] <- as.name("model.frame") # replace NULL with model.frame
  ##        Specials for smoothing
  mcall$formula <- if (missing(data)) {
    terms(formula, specials = .gamlss.sm.list)
  } else {
    terms(formula, specials = .gamlss.sm.list, data = data)
  }
  mu.frame <- eval(mcall, sys.parent()) # evalute the data.frame at the model.frame

  ## Declare family and deviance expression-----------------------------------------------------------------------------------------
  family <- as.gamlss.family(family) # bring first the gamlss family
  G.dev.expr <- body(family$G.dev.inc) # MS Thursday, April 11, 2002 at 10:34

  ## Extract model components -----------------------------------------------------------------------------------------
  ## This part deals with the response variable
  Y <- model.extract(mu.frame, "response") # extracting the y variable from the formula
  N <- max(length(Y), nrow(Y))
  ## extracting now the y and the binomial denominator in case we use BI or BB
  if (any(family$family %in% .gamlss.bi.list)) {
    if (NCOL(Y) == 1) {
      y <- if (is.factor(Y)) Y != levels(Y)[1] else Y
      bd <- rep(1, N)
      if (any(y < 0 | y > 1)) {
        stop("y values must be 0 <= y <= 1")
      }
    } else if (NCOL(Y) == 2) {
      if (any(abs(Y - round(Y)) > 0.001)) {
        warning("non-integer counts in a binomial GAMLSS!")
      }
      bd <- Y[, 1] + Y[, 2]
      y <- Y[, 1]
      if (any(y < 0 | y > bd)) stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005
    } else {
      stop(paste("For the binomial family, Y must be",
                 "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes",
                 "and col 2 is no. failures"))
    }
  } else if (any(family$family %in% .gamlss.multin.list)) {   # multinomial checking
    y <- if (is.factor(Y)) unclass(Y) else Y
  } else if (is.Surv(Y)) {   ## For censoring
    ## checking that the family is censored
    if (length(grep("censored", family$family[[2]])) == 0) {
      stop(paste("the family in not a censored distribution, use cens()"))
    }
    ## checking compatability of Surv object and censored distribution
    if (length(grep(attr(Y, "type"), family$family[[2]])) == 0) {
      stop(paste("the Surv object and the censored distribution are not of the same type"))
    }
    y <- Y
  } else {
    y <- Y
  }
  ## -----------------------------------------------------------------------------------------
  ## checking the permissible y values
  if (!family$y.valid(y)) { # MS Thursday, June 20, 2002 at 16:30
    stop("response variable out of range")
  }
  ## -----------------------------------------------------------------------------------------

  ## Extract the weights -----------------------------------------------------------------------------------------
  w <- model.extract(mu.frame, weights) # weights for the likelihood
  if (is.null(w)) {
    w <- rep(1, N)
  } else if (any(w < 0)) stop("negative weights not allowed")

  ##  Set up location-mean submodel: ====================================
  ##             mu.X   design matrix
  ##        mu.offset   offset in linear predictor
  ##         mu.start   starting values for mu (optional)
  mu.fit <- list() # MS Thursday, January 23, 2003 at 14:46
  mu.formula <- formula # ms Wednesday, December 29, 2004
  mu.terms <- attr(mu.frame, "terms") #   it peeks up the terms attribute
  mu.smoothers <- get.smoothers(mu.terms)
  mu.a <- attributes(mu.terms) #  from the model.frame
  mu.X <- model.matrix(mu.terms, mu.frame, contrasts) # the mean model matrix
  mu.offset <- model.extract(mu.frame, offset) # the mean-location offset
  if (is.null(mu.offset)) mu.offset <- rep(0, N)
  mu.object <- get.object("mu")
  formals(mu.object$dldp, envir = new.env()) <- alist(mu = fv) # this is to get the right GLIM arguments
  formals(mu.object$d2ldp2, envir = new.env()) <- alist(mu = fv) #
  formals(mu.object$G.di, envir = new.env()) <- alist(mu = fv) #
  formals(mu.object$valid, envir = new.env()) <- alist(mu = fv) #
  mu <- family$mu.linkinv(mu.X %*% dgamlss.update$mu.coef + mu.offset)

  ##
  ##  Set up dispersion-scale submodel: ---------------------------------------------------------------------------------------
  ##           sigma.X   design matrix
  ##           sigma.offset   offset in linear predictor
  ##       sigma.start   starting values for sigma (optional)
  if ("sigma" %in% names(family$parameters)) {
    orig.Envir <- attr(mcall$formula, ".Environment") # DS fix for Willem Thursday, March 18, 2010
    sigma.fit <- list() # MS Thursday, January 23, 2003 at 14:48
    form.sigma <- other.formula(form = sigma.formula)

    sigma.terms <- if (missing(data)) {
      terms(form.sigma, specials = .gamlss.sm.list)
    } else {
      terms(form.sigma, specials = .gamlss.sm.list, data = data)
    }
    mcall$formula <- sigma.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    sigma.frame <- eval(mcall, sys.parent())
    sigma.terms <- attr(sigma.frame, "terms")
    sigma.smoothers <- get.smoothers(sigma.terms)
    sigma.a <- attributes(sigma.terms)
    sigma.X <- model.matrix(sigma.terms, sigma.frame, contrasts)
    sigma.offset <- model.extract(sigma.frame, offset)
    if (is.null(sigma.offset)) sigma.offset <- rep(0, N)
    sigma.object <- get.object("sigma")
    formals(sigma.object$dldp, envir = new.env()) <- alist(sigma = fv) #
    formals(sigma.object$d2ldp2, envir = new.env()) <- alist(sigma = fv) #
    formals(sigma.object$G.di, envir = new.env()) <- alist(sigma = fv) #
    formals(sigma.object$valid, envir = new.env()) <- alist(sigma = fv) #
    formals(family$d2ldmdd, envir = new.env()) <- alist(sigma = sigma) #  ?? I do not think is needed
    sigma <- family$sigma.linkinv(sigma.X %*% dgamlss.update$sigma.coef + sigma.offset)
  }

  ##  Set up for the 3rd parameter submodel: ---------------------------------------------------------------------------------------
  ##            nu.X   design matrix
  ##       nu.offset   offset in linear predictor
  ##        nu.start   starting values for nu (optional)
  if ("nu" %in% names(family$parameters)) {
    nu.fit <- list() # MS Thursday, January 23, 2003 at 14:48
    form.nu <- other.formula(form = nu.formula)
    nu.terms <- if (missing(data)) {
      terms(form.nu, specials = .gamlss.sm.list)
    } else {
      terms(form.nu, specials = .gamlss.sm.list, data = data)
    }
    mcall$formula <- nu.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    nu.frame <- eval(mcall, sys.parent()) # ms Saturday, April 6, 2002 at 10:23
    nu.terms <- attr(nu.frame, "terms")
    nu.a <- attributes(nu.terms)
    nu.smoothers <- get.smoothers(nu.terms)
    nu.X <- model.matrix(form.nu, nu.frame, contrasts)
    nu.offset <- model.extract(nu.frame, offset)
    if (is.null(nu.offset)) nu.offset <- rep(0, N)
    nu.object <- get.object("nu")
    formals(nu.object$dldp, envir = new.env()) <- alist(nu = fv) # this is to get the right GLIM argument
    formals(nu.object$d2ldp2, envir = new.env()) <- alist(nu = fv) #
    formals(nu.object$G.di, envir = new.env()) <- alist(nu = fv) #
    formals(nu.object$valid, envir = new.env()) <- alist(nu = fv) #
    formals(family$d2ldmdv, envir = new.env()) <- alist(nu = nu)
    formals(family$d2ldddv, envir = new.env()) <- alist(nu = nu)
    nu <- family$nu.linkinv(nu.X %*% dgamlss.update$nu.coef + nu.offset)
  }

  ##  Set up for the 4rd parameter submodel:---------------------------------------------------------------------------------------
  ##            tau.X   design matrix
  ##       tau.offset   offset in linear predictor
  ##        tau.start   starting values for tau (optional)
  if ("tau" %in% names(family$parameters)) {
    tau.fit <- list() # MS Thursday, January 23, 2003 at 14:48
    form.tau <- other.formula(form = tau.formula)
    tau.terms <- if (missing(data)) {
      terms(form.tau, specials = .gamlss.sm.list)
    } else {
      terms(form.tau, specials = .gamlss.sm.list, data = data)
    }
    mcall$formula <- tau.terms
    attr(mcall$formula, ".Environment") <- orig.Envir # DS fix for Willem Thursday, March 18, 2010
    tau.frame <- eval(mcall, sys.parent())
    tau.terms <- attr(tau.frame, "terms") #
    tau.a <- attributes(tau.terms) #
    tau.smoothers <- get.smoothers(tau.terms)
    tau.X <- model.matrix(form.tau, tau.frame, contrasts)
    tau.offset <- model.extract(tau.frame, offset)
    if (is.null(tau.offset)) tau.offset <- rep(0, N)
    tau.object <- get.object("tau")
    formals(tau.object$dldp, envir = new.env()) <- alist(tau = fv) #
    formals(tau.object$d2ldp2, envir = new.env()) <- alist(tau = fv) #
    formals(tau.object$G.di, envir = new.env()) <- alist(tau = fv) #
    formals(tau.object$valid, envir = new.env()) <- alist(tau = fv) #
    formals(family$d2ldmdt, envir = new.env()) <- alist(tau = tau)
    formals(family$d2ldddt, envir = new.env()) <- alist(tau = tau)
    formals(family$d2ldvdt, envir = new.env()) <- alist(tau = tau)
    tau <- family$tau.linkinv(tau.X %*% dgamlss.update$tau.coef + tau.offset)
  }

  ## Check conditions for get_ssr
  if (get_ssr) {
    if (is.null(proposed_coef_list)) {
      stop("Must provide proposed_coef_list in order to get sums of squared residuals")
    }
  } else {
    if (!is.null(proposed_coef_list)) {
      warning("proposed_coef_list is provided, but get_ssr = FALSE. Ignoring proposed_coef_list.")
    }
  }
  ##  Start RS algorithm ================================
  environment(RS_local) <- environment()
  site_data <- RS_local()

  return(site_data)
}
