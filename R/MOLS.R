#' @title Mutually orthogonal latin square arrays of zize v
#'
#' @description
#' Constructs a list of orthogonal Latin square arrays of size v. 
#'
#' @details
#'  Returns at least 3 orthogonal latin square arrays of dimension v. The first and second arrays are the
#'  rows and columns of a v x v square and the third is a Latin square of size v x v. 
#'  If v is prime or prime power in the set  4, 8, 16, 32, 64, 128, 9, 27, 81, 25, 49 there
#'  are v-1 MOLS and \code{orthogLS(v)} will return a total of v+1 arrays. If v = 10 there are two MOLS
#'  and \code{orthogLS(10)} will give a total of 4 arrays ofsize 10 x 10.
#'  
#'  NB If S1 and T1 are mutually orthogonal Latin squares of order n1 and
#'  S2 and T2 are mutually orthogonal Latin squares of order n2 then the product squares S1xS2
#'  and T1xT2 are orthogonal to each other and have order n1n2 (not yet implemented)
#' 
#' @param v  the dimension of the required MOLS
#'
#' @return
#' \item{Treatments}{A table showing the replication number of each treatment in the design.}
#'
#' @references
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. 
#' https://CRAN.R-project.org/package=crossdes
#'
#' @examples
#' ## Set of 5 x 5 MOLS
#' orthogLS(5)
#' 
#' @export
#' @importFrom crossdes MOLS
  orthogLS=function(v) {
    if ((v*v)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(v*v))
      mols=MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])
      mols=lapply(seq(dim(mols)[3]), function(x){ mols[ , , x]-1})
      mols[[length(mols)+1]]=do.call(rbind, lapply(0:(v-1), function(j){ rep(0:(v-1))} ))
      mols[[length(mols)+1]]=t(mols[[length(mols)]])
      mols=mols[c(length(mols),    length(mols)-1 ,  1:(length(mols)-2))]
    } else if (v==10) {
        square1=matrix(c(1,    8,    9,    2,    0,    5,    3,    6,    4,     7,
                         9,    1,    0,    4,    2,    7,    8,    5,    3,     6,
                         0,    3,    1,    6,    8,    9,    7,    4,    2,     5,
                         3,    7,    4,    1,    5,    2,    6,    9,    8,     0,
                         8,    9,    5,    0,    1,    6,    4,    2,    7,     3,
                         2,    6,    3,    7,    4,    1,    5,    8,    0,     9,
                         5,    2,    6,    3,    7,    4,    1,    0,    9,     8,
                         4,    0,    7,    5,    3,    8,    9,    1,    6,     2,
                         7,    5,    8,    9,    6,    0,    2,    3,    1,     4,
                         6,    4,    2,    8,    9,    3,    0,    7,    5,     1),
                     nrow=10,ncol=10)
      square2 = matrix(c(1,    3,    5,    4,    2,    6,    7,    8,    9,     0,
                         3,    4,    8,    6,    7,    9,    1,    5,    0,     2,
                         5,    6,    7,    0,    9,    1,    4,    3,    2,     8,
                         9,    8,    1,    2,    3,    0,    5,    6,    7,     4,
                         2,    0,    4,    1,    6,    7,    8,    9,    5,     3,
                         8,    1,    2,    3,    0,    5,    9,    4,    6,     7,
                         0,    5,    9,    8,    1,    2,    3,    7,    4,     6,
                         4,    9,    6,    7,    5,    8,    2,    0,    3,     1,
                         7,    2,    0,    9,    4,    3,    6,    1,    8,     5,
                         6,    7,    3,    5,    8,    4,    0,    2,    1,     9),
                     nrow=10,ncol=10)
      mols=list(square1,square2)
      mols[[length(mols)+1]]=do.call(rbind, lapply(0:(v-1), function(j){ rep(0:(v-1))} ))
      mols[[length(mols)+1]]=t(mols[[length(mols)]])
      mols=mols[c(length(mols), length(mols)-1, 1:(length(mols)-2))]
    } else if (isPrime(v)) {
      mols=lapply(0:(v-1),function(z){do.call(rbind, lapply(0:(v-1), function(j){ (rep(0:(v-1))*z +j)%%v} ))})
      mols[[v+1]]=t(mols[[1]])
      mols=mols[c(1, length(mols) ,  2:(length(mols)-1))]
    } else {
      mols1= sapply(0:(v-1), function(j){ rep(0:(v-1))})
      mols2=t(mols1)
      mols3=  sapply(0:(v-1), function(j){ (rep(0:(v-1))+j)%%v} )
      mols=list(mols1,mols2,mols3)
    }
    return(mols)
  }