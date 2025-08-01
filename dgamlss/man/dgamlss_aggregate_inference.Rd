% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dgamlss_aggregate_inference.R
\name{dgamlss_aggregate_inference}
\alias{dgamlss_aggregate_inference}
\title{Aggregate inference across sites}
\usage{
dgamlss_aggregate_inference(
  site_list,
  pooled_coefs,
  lambda_list = NULL,
  penalty_matrix_list = NULL
)
}
\arguments{
\item{site_list}{List of joined outputs from dgamlss_aggregate_inference}

\item{pooled_coefs}{Named list of final pooled coefficients. Names must be "mu", "sigma", "nu", "tau".}

\item{lambda_list}{Named list of fixed lambdas, if fixed penalty smooth terms are desired. Names must be "mu", "sigma", "nu", and "tau".}

\item{penalty_matrix_list}{Named list of penalty matrices, output from dgamlss_bs(). Names must be "mu", "sigma", "nu", and "tau".}
}
\value{
List containing pooled Hessian, effective degrees of freedom, and overall number of samples.
}
\description{
Aggregate inference across sites
}
\examples{
\dontrun{
dgamlss_aggregate_inference(site_list, pooled_coef, lambda_list, penalty_matrix_list)
}
}
