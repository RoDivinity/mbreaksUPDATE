% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation.R
\name{nldat}
\alias{nldat}
\title{Computation of global minimizer for partial structural change model}
\usage{
nldat(
  y,
  z,
  x,
  h,
  m,
  p,
  q,
  bigT,
  fixb = 0,
  eps = 1e-05,
  maxi = 10,
  betaini = NULL,
  printd = 0
)
}
\arguments{
\item{y}{dependent variable in matrix form}

\item{z}{matrix of regressors which coefficients are allowed to change across regimes}

\item{x}{matrix of regressors which coefficients are constant across regime}

\item{h}{minimum segment length}

\item{m}{number of breaks}

\item{p}{number of `z` regressors}

\item{q}{number of `x` regressors}

\item{bigT}{the sample size T}

\item{fixb}{option to use initial \eqn{\beta} If \code{1}, procedure requires \code{betaini}.
If \code{0}, procedure will not use initial beta values}

\item{eps}{Convergence criterion (For partial change model ONLY)}

\item{maxi}{Maximum number of iterations (For partial change model ONLY)}

\item{betaini}{initial beta values. Required when use with option \code{fixb}}

\item{printd}{option to print output from iterated estimations. If \code{1}, the results
for each iteration will be printed in console log. If \code{0}, no output will be printed}
}
\value{
A list containing the following components:
\item{glb}{minimum global SSR}
\item{datevec}{Vector of dates (optimal minimizers)}
\item{bigvec}{Associated SSRs}
}
\description{
`nldat()` computes the break dates of a partial structural change model
for a pre-specified number of breaks `m`. The procedure iterates between
estimating the invariant and changing coefficients of `x` and `z` regressors
until convergence, by noting that the residuals from linear regression model between
`y` and `x` regressors is a pure structural change model,
while the residuals from pure structural change model between `y` and `z` regressors
is a linear regression
}
\references{
Bai J, Perron P (1998). \emph{"Estimating and Testing Linear Models with Multiple Structural
Changes"} Econometrica, 66, 47-78.
Bai J, Perron P (2003). \emph{"Computation and Analysis of Multiple Structural Change Models"}
Journal of Applied Econometrics 18, 1-22
}
