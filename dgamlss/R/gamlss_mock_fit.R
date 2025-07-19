gamlss_mock_fit <- function(formula = formula(data),
                            sigma.formula = ~1,
                            nu.formula = ~1,
                            tau.formula = ~1,
                            family = NO(),
                            data,
                            weights = NULL, # for weighted likelihood analysis
                            # (not the same as in GLM's)
                            contrasts = NULL, # one type of contrasts for all  parameters
                            method = RS(), # default algorithm
                            start.from = NULL, # starting from previous gamlss object
                            mu.start = NULL, # starting from given values
                            sigma.start = NULL,
                            nu.start = NULL,
                            tau.start = NULL,
                            mu.fix = FALSE, # whether the parameter is fixed
                            sigma.fix = FALSE,
                            nu.fix = FALSE,
                            tau.fix = FALSE,
                            control = gamlss.control(...),
                            i.control = glim.control(...), # the inner circle control (GLIM)
                            ...) {
  ## na.action = na.fail(), # both na.action and subset have been removed
  ##                    because while there is only one data set there are usually
  ##   subset = NULL,   four different model frames created therefore it is easier to apply
  ##                    sub-setting and na.action to the whole data set not to the
  ##                    frames
  ## ---------------------------------------------------------------------------------------
  ## require(stats) Thursday, June 10, 2004 at 09:58 MS
  # require(splines) # this will be removed with namespaces
  #----------------------------------------------------------------------------------------
  # gamlss.rc.list<-c("EX.rc","Exponential.rc") # the right censoring distribution list
  # gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial") # binomial denominators
  # .gamlss.multin.list<-c("MULTIN", "MN3", "MN4", "MN5")
  # ---------------------------------------------------------------------------------------
  # this is to replicate rqres within gamlss enviroment DS Friday, March 31, 2006 at 10:30
  rqres <- function(pfun = "pNO",
                    type = c("Continuous", "Discrete", "Mixed"),
                    censored = NULL,
                    ymin = NULL,
                    mass.p = NULL,
                    prob.mp = NULL,
                    y = y,
                    ...) { }
  body(rqres) <- eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  ## ---------------------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------------------
  ## first the definition of the three algorithms
  ## ---------------------------------------------------------------------------------------
  ## the mixing algorithm
  ## ---------------------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------------------

  ## =======================================================================================
  ## ---------------------------------------------------------------------------------------
  ## this function is used in for outputing the parameters
  ## ---------------------------------------------------------------------------------------
  parameterOut <- function(what = "mu", save) {
    out <- list()
    if (save == TRUE) {
      if (family$parameter[[what]] == TRUE && eval(parse(text = paste(what, ".fix", sep = ""))) == FALSE) {
        out$fv <- eval(parse(text = what))
        out$lp <- eval(parse(text = (paste(what, ".fit$eta", sep = ""))))
        out$wv <- eval(parse(text = (paste(what, ".fit$wv", sep = ""))))
        out$wt <- eval(parse(text = (paste(what, ".fit$wt", sep = ""))))
        out$link <- eval(parse(text = (paste(what, ".object$link", sep = ""))))
        out$terms <- eval(parse(text = (paste(what, ".terms", sep = ""))))
        out$x <- eval(parse(text = (paste(what, ".X", sep = ""))))
        out$qr <- eval(parse(text = (paste(what, ".fit$qr", sep = ""))))
        out$coefficients <- eval(parse(text = (paste(what, ".fit$coefficients", sep = ""))))
        out$offset <- eval(parse(text = (paste(what, ".offset", sep = ""))))
        out$xlevels <- .getXlevels(
          eval(parse(text = paste(what, ".terms", sep = ""))),
          eval(parse(text = paste(what, ".frame", sep = "")))
        )
        # ms Sunday, June 13 2004
        out$formula <- eval(parse(text = paste(what, ".formula", sep = "")))
        if (length(eval(parse(text = paste(what, ".smoothers", sep = "")))) > 0) {
          out$df <- eval(parse(text = paste(what, ".fit$nl.df", sep = ""))) +
            eval(parse(text = paste(what, ".fit$rank", sep = "")))
          out$nl.df <- eval(parse(text = paste(what, ".fit$nl.df", sep = "")))
          out$s <- eval(parse(text = paste(what, ".fit$smooth", sep = "")))
          out$var <- eval(parse(text = paste(what, ".fit$var", sep = "")))
          out$coefSmo <- eval(parse(text = paste(what, ".fit$coefSmo", sep = "")))
          out$lambda <- eval(parse(text = paste(what, ".fit$lambda", sep = "")))
          out$pen <- eval(parse(text = paste(what, ".fit$pen", sep = "")))
        } else {
          out$df <- eval(parse(text = paste(what, ".fit$rank", sep = "")))
          out$nl.df <- 0
          out$pen <- 0 # ms May 13, 2004
        }
      } else {
        out$fix <- eval(parse(text = paste(what, ".fix", sep = "")))
        out$df <- 0
        out$fv <- eval(parse(text = what))
      }
    } # if(save== not TRUE)
    else {
      if (family$parameter[[what]] == TRUE && eval(parse(text = paste(what, ".fix", sep = ""))) == FALSE) {
        if (length(eval(parse(text = paste(what, ".smoothers", sep = "")))) > 0) {
          out$df <- eval(parse(text = paste(what, ".fit$nl.df", sep = ""))) +
            eval(parse(text = paste(what, ".fit$rank", sep = "")))
          out$nl.df <- eval(parse(text = paste(what, ".fit$nl.df", sep = "")))
          out$terms <- eval(parse(text = (paste(what, ".terms", sep = ""))))
          out$formula <- eval(parse(text = paste(what, ".formula", sep = "")))
        } else {
          out$df <- eval(parse(text = paste(what, ".fit$rank", sep = "")))
          out$nl.df <- 0
          out$terms <- eval(parse(text = (paste(what, ".terms", sep = ""))))
          out$formula <- eval(parse(text = paste(what, ".formula", sep = "")))
        }
      } else {
        out$df <- 0
      }
    }
    out
  }
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
  ## =========================================================================================
  ## -----------------------------------------------------------------------------------------
  ## this function creates the parameter objects
  ## -----------------------------------------------------------------------------------------
  ## =========================================================================================
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
  # ==========================================================================================
  ## -----------------------------------------------------------------------------------------
  ## -----------------------------------------------------------------------------------------
  ## =========================================================================================
  ## here is where the proper gamlss function starts
  ## =========================================================================================
  ## -----------------------------------------------------------------------------------------
  ##       Save call for future reference
  gamlsscall <- match.call() #   the function call
  ## checking for NA in the data
  if (!missing(data)) {
    if (any(is.na(data))) {
      stop("The data contains NA's, use data = na.omit(mydata)")
    }
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

  ## -----------------------------------------------------------------------------------------
  ## This part deals with the family
  family <- as.gamlss.family(family) # bring first the gamlss family
  G.dev.expr <- body(family$G.dev.inc) # MS Thursday, April 11, 2002 at 10:34
  #  nopar <- family$nopar # the number of parameters for the family
  ## -----------------------------------------------------------------------------------------
  ## Now extract the model components using model.extra and model.matrix
  ## This part deals with the response variable
  Y <- model.extract(mu.frame, "response") # extracting the y variable from the formula
  if (is.null(dim(Y))) { # if y not matrix
    N <- length(Y)
  } else {
    N <- dim(Y)[1]
  } # calculate the dimension for y
  # .gamlss.bi.list <-  if (exists("gamlss.bi.list",envir=.GlobalEnv))
  #                         get("gamlss.bi.list", envir=.GlobalEnv) else .gamlss.bi.list
  ## extracting now the y and the binomial denominator in case we use BI or BB
  if (any(family$family %in% .gamlss.bi.list)) {
    if (NCOL(Y) == 1) {
      y <- if (is.factor(Y)) Y != levels(Y)[1] else Y
      bd <- rep(1, N)
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
    } else if (NCOL(Y) == 2) {
      if (any(abs(Y - round(Y)) > 0.001)) {
        warning("non-integer counts in a binomial GAMLSS!")
      }
      bd <- Y[, 1] + Y[, 2]
      y <- Y[, 1]
      if (any(y < 0 | y > bd)) stop("y values must be 0 <= y <= N") # MS Monday, October 17, 2005
    } else {
      stop(paste(
        "For the binomial family, Y must be",
        "a vector of 0 and 1's or a 2 column", "matrix where col 1 is no. successes",
        "and col 2 is no. failures"
      ))
    }
  }
  # multinomial checking
  else if (any(family$family %in% .gamlss.multin.list)) {
    y <- if (is.factor(Y)) {
      unclass(Y)
    } else {
      Y
    }
  }
  ## For censoring
  else if (is.Surv(Y)) {
    ## checking that the family is censored
    if (length(grep("censored", family$family[[2]])) == 0) {
      stop(paste("the family in not a censored distribution, use cens()"))
    }
    ## checking compatability of Surv object and censored distribution
    if (length(grep(attr(Y, "type"), family$family[[2]])) == 0) {
      stop(paste("the Surv object and the censored distribution are not of the same type"))
    }
    y <- Y
    # if (NCOL(Y) == 2)
    #   {
    #      #.event <- Y[,2]
    #      #    y  <- Y[,1]
    #      y <- Y
    #   }
    # else if (NCOL(Y) == 2)
    # stop("interval censored data are not implemented in gamlss yet")
  } else {
    y <- Y
  }
  ## -----------------------------------------------------------------------------------------
  ## checking the permissible y values
  if (!family$y.valid(y)) { # MS Thursday, June 20, 2002 at 16:30
    stop("response variable out of range")
  }
  ## -----------------------------------------------------------------------------------------
  ## this part is used if start.from is used as argument
  ## ------------start.from fitted model--------
  if (!is.null(start.from)) {
    if (!is.gamlss(start.from)) {
      stop(paste("The object in start.from is not a gamlss object", "\n", ""))
    }
    mu.start <- NULL
    sigma.start <- NULL
    nu.start <- NULL
    tau.start <- NULL
    ##               location model
    if ("mu" %in% start.from$parameters) {
      mu.start <- start.from$mu.fv
    }
    ##               scale-dispersion submodel
    if ("sigma" %in% start.from$parameters) {
      sigma.start <- start.from$sigma.fv
    }
    ##               nu submodel
    if ("nu" %in% start.from$parameters) {
      nu.start <- start.from$nu.fv
    }
    ##               tau submodel
    if ("tau" %in% start.from$parameters) {
      tau.start <- start.from$tau.fv
    }
  }
  ## -----------------------------------------------------------------------------------------
  ## checking the parameter.fix
  if (!is.logical(mu.fix)) stop("mu.fix should be logical TRUE or FALSE")
  if (!is.logical(sigma.fix)) stop("sigma.fix should be logical TRUE or FALSE")
  if (!is.logical(nu.fix)) stop("nu.fix should be logical TRUE or FALSE")
  if (!is.logical(tau.fix)) stop("tau.fix should be logical TRUE or FALSE")
  ## -----------------------------------------------------------------------------------------
  ## extract the weights
  w <- model.extract(mu.frame, weights) # weights for the likelihood
  if (is.null(w)) {
    w <- rep(1, N)
  } else if (any(w < 0)) stop("negative weights not allowed") #
  #   else if (!all(trunc(w)==w)) warning("weights should be integer values \n",
  #         " indicating number of observations with identical values \n") #
  ## =========================================================================================
  ## -----------------------------------------------------------------------------------------
  ##  Set up location-mean submodel:
  ##             mu.X   design matrix
  ##        mu.offset   offset in linear predictor
  ##         mu.start   starting values for mu (optional)
  ## -----------------------------------------------------------------------------------------
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
  ## initial values for mu
  if (!is.null(mu.start)) {
    mu <- if (length(mu.start) > 1) mu.start else rep(mu.start, N)
  } else {
    (eval(family$mu.initial))
  } # MS: Friday, March 29, 2002 at 11:27
  ## ---------------------------------------------------------------------------------------
  ##  Set up dispersion-scale submodel:
  ##           sigma.X   design matrix
  ##           sigma.offset   offset in linear predictor
  ##       sigma.start   starting values for sigma (optional)
  ## ---------------------------------------------------------------------------------------
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
    ## initial values for sigma
    if (!is.null(sigma.start)) {
      sigma <- if (length(sigma.start) > 1) sigma.start else rep(sigma.start, N)
    } else {
      eval(family$sigma.initial)
    } #
  }
  ## -----------------------------------------------------------------------------------------
  ##  Set up for the 3rd parameter submodel:
  ##            nu.X   design matrix
  ##       nu.offset   offset in linear predictor
  ##        nu.start   starting values for nu (optional)
  ## -----------------------------------------------------------------------------------------
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
    ## initial values for nu
    if (!is.null(nu.start)) {
      nu <- if (length(nu.start) > 1) nu.start else rep(nu.start, N)
    } else {
      eval(family$nu.initial)
    } #
  }
  ## -----------------------------------------------------------------------------------------
  ##  Set up for the 4rd parameter submodel:
  ##            tau.X   design matrix
  ##       tau.offset   offset in linear predictor
  ##        tau.start   starting values for tau (optional)
  ## -----------------------------------------------------------------------------------------
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
    ## initial values for tau
    if (!is.null(tau.start)) {
      tau <- if (length(tau.start) > 1) tau.start else rep(tau.start, N)
    } else {
      eval(family$tau.initial)
    } #
  }

  ## -----------------------------------------------------------------------------------------
  ## -----------------------------------------------------------------------------------------
  ## =========================================================================================
  ##  Checking whether proper algorithm  (RS, CG or mixed)
  ## =========================================================================================
  ## -----------------------------------------------------------------------------------------
  name.method <- substitute(method)
  name.method <- deparse(name.method[1])
  list.methods <- c("RS()", "CG()", "mixed()")
  i.method <- pmatch(name.method, list.methods, nomatch = 0)
  if (!i.method) stop("Method must be RS(), CG() or mixed()")
  ## -----------------------------------------------------------------------------------------
  ## fitting the model
  fiter <- 0
  conv <- TRUE
  method <- substitute(method)
  ## -----------------------------------------------------------------------------------------
  ## -----------------------------------------------------------------------------------------
  ##  Getting the GAMLSS object out
  ## ----------------------------------------------------------------------------------------
  ## first the general output
  ## calculate the Global deviance again
  G.dev.incr <- eval(G.dev.expr)
  G.dev <- sum(w * G.dev.incr)
  out <- list(
    family = family$family, parameters = names(family$parameters),
    call = gamlsscall, y = y, control = control, weights = w,
    G.deviance = G.dev, N = N, rqres = family$rqres, iter = fiter,
    type = family$type, method = method, contrasts = contrasts
  )
  # , na.action=na.act
  out$converged <- conv
  out$residuals <- eval(family$rqres)
  noObs <- if (all(trunc(w) == w)) sum(w) else N
  out$noObs <- noObs
  ## binomial denominator
  if (any(family$family %in% .gamlss.bi.list)) out$bd <- bd
  ## -----------------------------------------------------------------------------------------
  saveParam <- control$save
  ##  Output for mean model: ----------------------------------------------------------------
  if ("mu" %in% names(family$parameters)) {
    out <- c(out, mu = parameterOut(what = "mu", save = saveParam))
  } else {
    out$mu.df <- 0
  }
  ## define now the degrees of freedom for the fit and residuals
  out$df.fit <- out$mu.df
  out$pen <- out$mu.pen
  out$df.residual <- noObs - out$mu.df
  ## Output for dispersion model: ----------------------------------------------------------
  if ("sigma" %in% names(family$parameters)) {
    out <- c(out, sigma = parameterOut(what = "sigma", save = saveParam))
    out$df.fit <- out$mu.df + out$sigma.df
    out$pen <- out$mu.pen + out$sigma.pen
    out$df.residual <- noObs - out$mu.df - out$sigma.df
  }
  ##  output for nu ------------------------------------------------------------------------
  if ("nu" %in% names(family$parameters)) {
    out <- c(out, nu = parameterOut(what = "nu", save = saveParam))
    out$df.fit <- out$mu.df + out$sigma.df + out$nu.df
    out$df.residual <- noObs - out$mu.df - out$sigma.df - out$nu.df
    out$pen <- out$mu.pen + out$sigma.pen + out$nu.pen
  }
  ##  output for tau -----------------------------------------------------------------------
  if ("tau" %in% names(family$parameters)) {
    out <- c(out, tau = parameterOut(what = "tau", save = saveParam))
    out$df.fit <- out$mu.df + out$sigma.df + out$nu.df + out$tau.df
    out$pen <- out$mu.pen + out$sigma.pen + out$nu.pen + out$tau.pen
    out$df.residual <- noObs - out$mu.df - out$sigma.df - out$nu.df - out$tau.df
  }
  ## =======================================================================================
  out$P.deviance <- out$G.deviance + out$pen # ms Thursday, May 13, 2004
  out$aic <- G.dev + 2 * out$df.fit
  out$sbc <- G.dev + log(noObs) * out$df.fit
  # MS Thursday, April 22, 2004 at 11:26
  #  if ((ls(1,pattern="fiter")=="fiter")) rm(fiter, envir = as.environment(1))
  # MS Thursday, January 8, 2004 at 17:52
  class(out) <- c("gamlss", "gam", "glm", "lm")
  out
}
