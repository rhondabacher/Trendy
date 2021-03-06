% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{breakpointFit}
\alias{breakpointFit}
\title{break point fits}
\usage{
breakpointFit(J, tVectIn, lmLinear, numTry)
}
\arguments{
\item{J}{number of breakpoints in the model}

\item{tVectIn}{a numerical vector indicating the time-points or the 
order of samples. If it is NULL (default), then the time/order 
will be assumed to be equaly spaced from 1:N (N is number of samples).}

\item{lmLinear}{the linear model fit; no breakpoints}

\item{numTry}{the number of different seeds to try. If all numTry runs 
fail, then the linear regression (no breakpoints, one segment) model 
will be returned.}
}
\description{
break point fits
}
