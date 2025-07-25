% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgamlss_get_inits.R
\name{dgamlss_get_inits}
\alias{dgamlss_get_inits}
\title{Get initial values based on running gamlss on the coordinating site data.}
\usage{
dgamlss_get_inits(local_gamlss)
}
\arguments{
\item{local_gamlss}{gamlss output from local site}
}
\value{
List of initial coefficients pulled from the local_gamlss object. Since the local_gamlss design matrix may be singular even if the pooled design matrix would not be, NAs are substituted by 0s.
}
\description{
Get initial values based on running gamlss on the coordinating site data.
}
\examples{
\dontrun{
data(abdom)
site1_data <- abdom[1:110, ]
site2_data <- abdom[111:610, ]
site_data <- list(site1_data, site2_data)
site1_gamlss <- gamlss(mu.formula = y ~ x1 + x2, sigma.formula = ~ x1,
nu.formula = ~ 1, tau.formula = ~ 1, family = BCT, data = site1_data)
initial_values <- dgamlss_get_inits(site1_gamlss)
}
}
