dgamlss_predict <- function(dgamlss_object, local_object) {
  gamlss_family <- as.gamlss.family(dgamlss_object$family[1])
  dgamlss_object$y <- local_object$y
  dgamlss_object$mu.x <- local_object$mu.x
  dgamlss_object$sigma.x <- local_object$sigma.x
  dgamlss_object$nu.x <- local_object$nu.x
  dgamlss_object$tau.x <- local_object$tau.x

  dgamlss_prediction <- list(y = local_object$y,
                             mu.fv = gamlss_family$mu.linkinv(dgamlss_object$mu.x %*% dgamlss_object$mu.coefficients),
                             sigma.fv = gamlss_family$sigma.linkinv(dgamlss_object$sigma.x %*% dgamlss_object$sigma.coefficients),
                             nu.fv = gamlss_family$nu.linkinv(dgamlss_object$nu.x %*% dgamlss_object$nu.coefficients),
                             tau.fv = gamlss_family$tau.linkinv(dgamlss_object$tau.x %*% dgamlss_object$tau.coefficients))


  dFamily <- match.fun(paste0("d", gamlss_family$family[1]))
  pFamily <- match.fun(paste0("p", gamlss_family$family[1]))
  qFamily <- match.fun(paste0("q", gamlss_family$family[1]))
  dgamlss_prediction$median <- qFamily(0.5, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                                       dgamlss_prediction$nu.fv, dgamlss_prediction$tau.fv)
  dgamlss_prediction$quantile <- pFamily(dgamlss_prediction$y, dgamlss_prediction$mu.fv, dgamlss_prediction$sigma.fv,
                                       dgamlss_prediction$nu.fv, dgamlss_prediction$tau.fv)

  return(dgamlss_prediction)
}
