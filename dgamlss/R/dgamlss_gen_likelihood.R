#  Mikis Stasinopoulos 23-11-12
#  this function generates the likelihood of a gamlss object
#  It takes an gamlss object and generates a function with arguments the parameters
#  of the model.
#  It is created so we can replace the current vcov() function which fails in a lot of occations
#  at the moment  parametric model are allowed where the additve are treated as fixed
#  I would be nice if the lambda are also used so we can get standard errors for them
#  In order to do that the X W and G matrices have to be saved after bp()
# TO DO
# i)   check what happents with offset's?? (OK 2-2-13)
# ii)  create a vcov function to replace the (OK done 12/12/12)
# iii)  maybe the function has to have as default the current fitted beta values (OK)
# iv)  what happents with cencored/ truncated/ log logit distributions???
# v)   what we do is the inverse fails?
# vi) fixed distribution parameters fails
# vii) what happends if x+pb(x) or
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
dgamlss_gen_likelihood <- function(object) {
  if (!is.gamlss(object)) stop("needs a gamlss object")
  #    fam <- as.gamlss.family(object$family) this is changed to get the right link functions
  fam <- if (is.null(object$call$family)) {
    as.gamlss.family(NO)
  } else {
    as.gamlss.family(object$call$family)
  }
  fname <- object$family[1]
  dfun <- paste("d", fname, sep = "")
  # pfun <- paste("p",fname,sep="")
  pdf <- eval(parse(text = dfun))
  # cdf <- eval(parse(text=pfun)) # may be we need this for censored data
  nopar <- length(object$par)
  y <- object$y
  w <- object$weights
  X <- list()
  links <- list()
  coefs <- list()
  Smo <- list()
  offSet <- list()
  binomialIndex <- FALSE
  if (fname %in% .gamlss.bi.list) {
    bd <- object$bd # binomial denominator
    binomialIndex <- TRUE
  }
  # # establish whether data is called
  # if (any(grepl("data", names(object$call)))) {
  #   exitData <- TRUE
  #   ##    DS The idea here is if na.omit(data) is used to do the correct thing
  #   DaTa <- if (startsWith(as.character(object$call["data"]), "na.omit")) {
  #     eval(parse(text = as.character(object$call["data"])))
  #   } else {
  #     eval(parse(as.character(object$call["data"])))
  #   }
  # } else {
  #   exitData <- FALSE
  # }
  for (i in object$par)
  {
    if (length(eval(parse(text = paste(paste("object$", i, sep = ""), ".fix==TRUE", sep = "")))) != 0) {
      ff <- eval(parse(text = paste(paste(paste(object$family[1], "()$", sep = ""), i, sep = ""), ".linkfun", sep = "")))
      fixvalue <- ff(fitted(object, i)[1])
      names(fixvalue) <- paste("fixed", i, sep = " ")
      coefs[[i]] <- fixvalue
      X[[i]] <- matrix(rep(1, length(object$y)), ncol = 1)

      links[i] <- eval(parse(text = paste(paste(paste(object$family[1], "()$", sep = ""), i, sep = ""), ".link", sep = "")))
      Smo[[i]] <- matrix(rep(0, length(object$y)), ncol = 1)
      offSet[[i]] <- rep(0, length(object$y))
    } else {
      coefs[[i]] <- eval(parse(text = paste(paste("object$", i, sep = ""), ".coefficients", sep = "")))
      notNAcoef <- !is.na(coefs[[i]])
      X[[i]] <- eval(parse(text = paste(paste("object$", i, sep = ""), ".x", sep = "")))
      if (any(is.na(coefs[[i]]))) { # this is to ensure that are non NA
        coefs[[i]] <- coefs[[i]][notNAcoef]
        X[[i]] <- as.matrix(X[[i]][, notNAcoef])
      }
      links[i] <- eval(parse(text = paste(paste("object$", i, sep = ""), ".link", sep = "")))
      Smo[[i]] <- eval(parse(text = paste(paste("object$", i, sep = ""), ".s", sep = "")))
      offSet[[i]] <- eval(parse(text = paste(paste("object$", i, sep = ""), ".offset", sep = "")))
    }
  }
  switch(nopar,
         { #  1 parameter
           lik.fun <- function(par) {
             lmu <- length(coefs[["mu"]])
             if (length(par) != lmu) stop("par is not the right length")
             eta.mu <- if (is.null(Smo[["mu"]])) {
               X[["mu"]] %*% par[1:lmu] + offSet[["mu"]]
             } else {
               X[["mu"]] %*% par[1:lmu] + rowSums(Smo[["mu"]]) + offSet[["mu"]]
             }
             mu <- fam$mu.linkinv(eta.mu)
             if (binomialIndex) {
               -sum(w * pdf(y, mu = mu, bd = bd, log = TRUE))
             } else {
               -sum(w * pdf(y, mu = mu, log = TRUE))
             }
           }
           thebetas <- coefs[["mu"]]
           formals(lik.fun) <- alist(par = thebetas)
         },
         { # 2 parameter
           lik.fun <- function(par) {
             lmu <- length(coefs[["mu"]])
             lsigma <- length(coefs[["sigma"]])
             tl <- lmu + lsigma
             if (length(par) != tl) stop("par is not the right length")
             eta.mu <- if (is.null(Smo[["mu"]])) {
               X[["mu"]] %*% par[1:lmu] + offSet[["mu"]]
             } else {
               X[["mu"]] %*% par[1:lmu] + rowSums(Smo[["mu"]]) + offSet[["mu"]]
             }
             mu <- fam$mu.linkinv(eta.mu)
             eta.sigma <- if (is.null(Smo[["sigma"]])) {
               X[["sigma"]] %*% par[(lmu + 1):(lmu + lsigma)] + offSet[["sigma"]]
             } else {
               X[["sigma"]] %*% par[(lmu + 1):(lmu + lsigma)] + rowSums(Smo[["sigma"]]) + offSet[["sigma"]]
             }
             sigma <- fam$sigma.linkinv(eta.sigma)
             if (binomialIndex) {
               -sum(w * pdf(y, mu = mu, sigma = sigma, bd = bd, log = TRUE))
             } else {
               -sum(w * pdf(y, mu = mu, sigma = sigma, log = TRUE))
             }
           }
           thebetas <- c(coefs[["mu"]], coefs[["sigma"]])
           formals(lik.fun) <- alist(par = thebetas)
         },
         { # 3 parameter
           lik.fun <- function(par) {
             lmu <- length(coefs[["mu"]])
             lsigma <- length(coefs[["sigma"]])
             lnu <- length(coefs[["nu"]])
             tl <- lmu + lsigma + lnu
             if (length(par) != tl) stop("par is not the right length")
             eta.mu <- if (is.null(Smo[["mu"]])) {
               X[["mu"]] %*% par[1:lmu] + offSet[["mu"]]
             } else {
               X[["mu"]] %*% par[1:lmu] + rowSums(Smo[["mu"]]) + offSet[["mu"]]
             }
             mu <- fam$mu.linkinv(eta.mu)
             eta.sigma <- if (is.null(Smo[["sigma"]])) {
               X[["sigma"]] %*% par[(lmu + 1):(lmu + lsigma)] + offSet[["sigma"]]
             } else {
               X[["sigma"]] %*% par[(lmu + 1):(lmu + lsigma)] + rowSums(Smo[["sigma"]]) + offSet[["sigma"]]
             }
             sigma <- fam$sigma.linkinv(eta.sigma)
             eta.nu <- if (is.null(Smo[["nu"]])) {
               X[["nu"]] %*% par[(lmu + lsigma + 1):(lmu + lsigma + lnu)] + offSet[["nu"]]
             } else {
               X[["nu"]] %*% par[(lmu + lsigma + 1):(lmu + lsigma + lnu)] + rowSums(Smo[["nu"]]) + offSet[["nu"]]
             }
             nu <- fam$nu.linkinv(eta.nu)
             if (binomialIndex) {
               -sum(w * pdf(y, mu = mu, sigma = sigma, nu = nu, bd = bd, log = TRUE))
             } else {
               -sum(w * pdf(y, mu = mu, sigma = sigma, nu = nu, log = TRUE))
             }
           }
           thebetas <- c(coefs[["mu"]], coefs[["sigma"]], coefs[["nu"]])
           formals(lik.fun) <- alist(par = thebetas)
         },
         { # 4 parameter
           lik.fun <- function(par) {
             lmu <- length(coefs[["mu"]])
             lsigma <- length(coefs[["sigma"]])
             lnu <- length(coefs[["nu"]])
             ltau <- length(coefs[["tau"]])
             tl <- lmu + lsigma + lnu + ltau
             if (length(par) != tl) stop("par is not the right length")
             eta.mu <- if (is.null(Smo[["mu"]])) {
               X[["mu"]] %*% par[1:lmu] + offSet[["mu"]]
             } else {
               X[["mu"]] %*% par[1:lmu] + rowSums(Smo[["mu"]]) + offSet[["mu"]]
             }
             mu <- fam$mu.linkinv(eta.mu)
             eta.sigma <- if (is.null(Smo[["sigma"]])) {
               X[["sigma"]] %*% par[(lmu + 1):(lmu + lsigma)] + offSet[["sigma"]]
             } else {
               X[["sigma"]] %*% par[(lmu + 1):(lmu + lsigma)] + rowSums(Smo[["sigma"]]) + offSet[["sigma"]]
             }
             sigma <- fam$sigma.linkinv(eta.sigma)
             eta.nu <- if (is.null(Smo[["nu"]])) {
               X[["nu"]] %*% par[(lmu + lsigma + 1):(lmu + lsigma + lnu)] + offSet[["nu"]]
             } else {
               X[["nu"]] %*% par[(lmu + lsigma + 1):(lmu + lsigma + lnu)] + rowSums(Smo[["nu"]]) + offSet[["nu"]]
             }
             nu <- fam$nu.linkinv(eta.nu)
             eta.tau <- if (is.null(Smo[["tau"]])) {
               X[["tau"]] %*% par[(lmu + lsigma + lnu + 1):(lmu + lsigma + lnu + ltau)] + offSet[["tau"]]
             } else {
               X[["tau"]] %*% par[(lmu + lsigma + lnu + 1):(lmu + lsigma + lnu + ltau)] + rowSums(Smo[["tau"]]) + offSet[["tau"]]
             }
             tau <- fam$tau.linkinv(eta.tau)
             -sum(w * pdf(y, mu = mu, sigma = sigma, nu = nu, tau = tau, log = TRUE))
           }
           thebetas <- c(coefs[["mu"]], coefs[["sigma"]], coefs[["nu"]], coefs[["tau"]])
           formals(lik.fun) <- alist(par = thebetas)
         }
  )
  lik.fun
}
#-------------------------------------------------------------------------------
