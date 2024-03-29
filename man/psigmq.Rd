% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{psigmq}
\alias{psigmq}
\title{Construct diagonal matrix of estimated variance}
\usage{
psigmq(res, b, q, m, nt)
}
\arguments{
\item{res}{big residual vector of the model}

\item{b}{Estimated date of changes}

\item{q}{Number of `z` regressors}

\item{m}{Number of breaks}

\item{nt}{The size of `z` regressors}
}
\value{
sigmat (`m`+1)x(`m`+1) diagonal matrix with i-th entry
equal to estimated variance of regime i
}
\description{
Function computes a diagonal matrix of dimension m+1 by m+1
with i-th entry is the estimated variance of residuals of segment i
}
