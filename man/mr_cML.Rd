% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_cML.R
\name{mr_cML}
\alias{mr_cML}
\title{MRcML method}
\usage{
mr_cML(
  b_exp,
  b_out,
  se_exp,
  se_out,
  K_vec = 0:(length(b_exp) - 2),
  random_start = 0,
  maxit = 100,
  random_seed = 0,
  n
)
}
\arguments{
\item{b_exp}{Vector of estimated effects for exposure.}

\item{b_out}{Vector or estimated effects for outcome.}

\item{se_exp}{Vector of standard errors for exposure.}

\item{se_out}{Vector of standard errors for outcome.}

\item{K_vec}{Sets of candidate K's, the constraint parameter representing number of invalid IVs.}

\item{random_start}{Number of random start points for cML, default is 0.}

\item{maxit}{Maximum number of iterations for each optimization.}

\item{random_seed}{Random seed, an integer. Default is
0, which does not set random seed; user could specify a positive integer
as random seed to get replicable results.}

\item{n}{Sample size.}
}
\value{
A list contains full results of cML methods.
MA_BIC_theta, MA_BIC_se, MA_BIC_p:
Estimate of theta,
its standard error and p-value from cML-MA-BIC.
Similarly for BIC_theta, BIC_se, BIC_p from cML-BIC;
for MA_AIC_theta, MA_AIC_se, MA_AIC_p from cML-MA-AIC;
for AIC_DP_theta, AIC_DP_se, AIC_DP_p from cML-AIC.
BIC_invalid is the set of invalid IVs selected by cML-BIC,
AIC_invalid is the set of invalid IVs selected by cML-AIC.
BIC_vec is the BIC vector.
}
\description{
This is the main function of MRcML method, without data perturbation.
}
