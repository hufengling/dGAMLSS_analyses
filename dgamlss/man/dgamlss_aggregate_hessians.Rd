% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgamlss_aggregate_hessian.R
\name{dgamlss_aggregate_hessians}
\alias{dgamlss_aggregate_hessians}
\title{Aggregate local Hessians into pooled Hessian}
\usage{
dgamlss_aggregate_hessians(hessian_list)
}
\arguments{
\item{hessian_list}{List of local Hessians produced by dgamlss_get_hessian()}
}
\value{
Pooled hessian
}
\description{
Aggregate local Hessians into pooled Hessian
}
\examples{

\dontrun{
# Assuming site1_gamlss, site2_gamlss, and pooled_coefs are already defined

site1_hess <- dgamlss_get_hessian(site1_gamlss, pooled_coefs)
site2_hess <- dgamlss_get_hessian(site2_gamlss, pooled_coefs)
pooled_hess <- dgamlss_aggregate_hessian(list(site1_hess, site2_hess))
}
}
