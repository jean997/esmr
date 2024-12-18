% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grad_log_py.R
\name{optimize_lpy2}
\alias{optimize_lpy2}
\title{One step update}
\usage{
optimize_lpy2(
  fit,
  max_steps = 1,
  tol = 1e-05,
  calc_hess = FALSE,
  max_components = Inf,
  max_prob = 1,
  sub_size = fit$n
)
}
\arguments{
\item{fit}{ESMR fit object}

\item{max_steps}{Number of steps to take (default is 1)}

\item{tol}{Tolerance for convergence if max_steps is > 1}

\item{calc_hess}{Calculate the hesssian? (this can be slow)}

\item{max_components}{Maximum number of likelihood components to use in approximation}

\item{max_prob}{Maximum integral of approximate likelihood}

\item{sub_size}{Size of subset to use to calculate the gradient of the likelihood.}
}
\description{
One step update
}
\details{
This function computes a one step update from the variational solution provided by ESMR. By default
the gradient of the likelihood is computed exactly. This should work well if the number of traits is < 15.
For larger numbers of traits, you may want to use an approximation to the likelihood controlled by parameters max_components,
max_prob and sub_size.
}
