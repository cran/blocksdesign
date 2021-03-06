% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fraction.R
\name{fraction}
\alias{fraction}
\title{Optimum treatment set from a candidate set of treatments}
\usage{
fraction(
  treatments,
  size,
  treatments_model = NULL,
  restriction_model = NULL,
  searches = 50
)
}
\arguments{
\item{treatments}{is a data frame or a list containing a candidate set of factorial treatments.}

\item{size}{is the required number of treatments in the fractional set of treatments.}

\item{treatments_model}{is a model formula for the required treatments design.}

\item{restriction_model}{is a model formula which is a subset of the \code{treatments_model}
formula and which fixes those treatment factors contained in the restriction model formula.}

\item{searches}{are the maximum number of searches for selecting the best optimization.}
}
\value{
A list containing:
 \itemize{
 \item{"TF"} {An optimized treatment fraction of the required \code{size}}
 \item{"fullTF"} {The full candidate set of treatments but with the first \code{size} 
 rows containing the optimized treatment fraction}
 \item{"Efficiency"} {The D-efficiency of the optimized treatment fraction relative to the full candidate set of treatments}
}
}
\description{
Finds an optimum set of treatments from a candidate set of treatments for any arbitrary 
treatments design formula.
}
\details{
The candidate set \code{treatments} will normally contain one or more complete sets of 
treatment replicates. The algorithm re-arranges the rows of the \code{treatments} set to ensure 
that the first \code{size} rows of the optimized \code{treatments} set contains the optimized treatment
fraction. The maximum replication of any treatment is the number of times that treatment occurs in the 
candidate treatment set and for a polynomial response surface design extra replication of the candidate
set may be necessary to allow for differential replication of the design points. The design is
optimized with respect to the \code{treatments_model} conditional on the treatment factors
in the \code{restriction_model} being held constant. The \code{restriction_model} must be a subset
of the full \code{treatments_model} otherwise the design will be fully fixed and no further optimization
will be possible. Fitting a non-null \code{restriction_model} allows sequential optimization
with each successively \code{treatments_model} optimized conditional on all previously optimized models.
The D-optimal efficiency of the design for the optimized treatment set is calculated relative to the 
D-optimal efficiency of the design for the candidate treatment set. 

The default \code{treatments_model} parameter is an additive model for all treatment factors.
}
\examples{
#' ## Plackett and Burman (P&B) type design for eleven 2-level factors in 12 runs 
## NB. The algorithmic method is unlikely to succeed for larger P&B type designs. 

GF = list(F1 = factor(1:2,labels=c("a","b")), F2 = factor(1:2,labels=c("a","b")), 
                 F3 = factor(1:2,labels=c("a","b")), F4 = factor(1:2,labels=c("a","b")),
                 F5 = factor(1:2,labels=c("a","b")), F6 = factor(1:2,labels=c("a","b")),
                 F7 = factor(1:2,labels=c("a","b")), F8 = factor(1:2,labels=c("a","b")), 
                 F9 = factor(1:2,labels=c("a","b")), F10= factor(1:2,labels=c("a","b")), 
                 F11= factor(1:2,labels=c("a","b")) )
model = ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9 + F10 + F11
Z=fraction(GF,size=12,treatments_model=model,searches=100)
print(Z$TF)
print(Z$Efficiency)
round(crossprod(scale(data.matrix(Z$TF))),6)

## Factorial treatment designs defined by sequentially fitted factorial treatment models
## 4 varieties by 3 levels of N by 3 levels of K assuming degree-2 treatment model in 24 plots.
## The single stage model gives an unequal split for the replication of the four varieties
## whereas the two stage model forces an equal split of 6 plots per variety.
## The single stage model is slightly more efficient overall (about 1.052045 versus 1.043662)
## but unequal variety replication is undesirable if all varieties are equally important.

## model terms
treatments = list(Variety = factor(1:4), N = 1:3, K = 1:3)
variety_model = ~ Variety
full_model = ~ (Variety + N + K)^2  + I(N^2) + I(K^2)

## single stage model
opt_full_treatments = fraction(treatments,24,full_model,searches=10)
opt_full_treatments$Efficiency
table(opt_full_treatments$TF[,1]) # variety replication

## two stage model
opt_var_treatments  = fraction(treatments,24,variety_model,searches=10)
opt_full_treatments = fraction(opt_var_treatments$fullTF,24,full_model,variety_model,searches=10)
opt_full_treatments$Efficiency
table(opt_full_treatments$TF[,1]) # variety replication

}
\seealso{
\code{\link{design}}
}
