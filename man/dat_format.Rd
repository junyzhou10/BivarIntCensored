% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dat_format.R
\name{dat_format}
\alias{dat_format}
\title{Format the input data}
\usage{
dat_format(Ti, ni = NULL)
}
\arguments{
\item{Ti}{Data from one outcome}

\item{ni}{Naming the outcome number, i.e., return the column names as tu/tv_<ni>. If not provided, it takes the name from input data.}
}
\value{
A formatted data with censoring type indicators
}
\description{
Format the input data
}
\details{
See Detail section for \code{BiIntCensd}
}
