% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/A_bound.R
\name{A_bound}
\alias{A_bound}
\title{Efficiency bounds}
\usage{
A_bound(n, v, b)
}
\arguments{
\item{n}{the total number of plots in the design.}

\item{v}{the total number of treatments in the design.}

\item{b}{the total number of blocks in the design.}
}
\description{
Finds upper A-efficiency bounds for regular block designs.
}
\details{
Upper bounds for the A-efficiency factors of regular block designs with equi-replicate
treatments and equal block sizes (see Chapter 2.8 of John and Williams 1995). 
Non-trivial A-efficiency upper bounds are calculated for regular block designs with equal block sizes 
and equal replication only. All other designs return NA.
}
\examples{

# 50 plots, 10 treatments and 10 blocks for a design with 5 replicates and blocks of size 5 
A_bound(n=50,v=10,b=10)

}
\references{
John, J. A. and Williams, E. R. (1995). Cyclic and Computer Generated Designs. Chapman and Hall, London.
}
