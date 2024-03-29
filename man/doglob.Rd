% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Procedures.R
\name{doglob}
\alias{doglob}
\title{Global SSR minimizer for structural change model}
\usage{
doglob(y, z, x, m, eps, h, maxi, fixb, betaini, printd, eps1 = 0.15)
}
\arguments{
\item{y}{matrix of dependent variable}

\item{z}{matrix of independent variables with coefficients allowed to change across
regimes}

\item{x}{matrix of independent variables with coefficients constant across regimes}

\item{m}{number of breaks in the structural change model}

\item{eps}{convergence criterion for iterative recursive computation. (For partial change model ONLY)}

\item{h}{Minimum segment length of regime considered in estimation. If users want to specify a particular value, please set `eps1=0`}

\item{maxi}{maximum number of iterations. (For partial change model ONLY)}

\item{fixb}{option to use fixed initial input \eqn{\beta}. If \code{1},
the model will use values given in \code{betaini}. If \code{0}, betaini is skipped}

\item{betaini}{Initial \eqn{beta_0} to use in estimation (Must be a `p x 1` matrix, where p is number of x variables)}

\item{printd}{Print option for model estimation. \code{default} = 0, to
suppress intermediate outputs printing to console}

\item{eps1}{trimming level. The default value for `eps1` is set at 0.15.}
}
\value{
A list containing the following components:
\itemize{
\item{glb} {Minimum global SSR}
\item{datevec} {Vector of dates (optimal minimizers)}
\item{bigvec} {Associated SSRs with possible break dates combination}}
}
\description{
`doglob()` identify if the structural change model is i) pure or ii)
partial change model. The procedure then calls appropriate functions \code{\link[mbreaks]{dating}} to estimate
the pure change model and \code{\link[mbreaks]{nldat}} to estimate the partial change model.
}
