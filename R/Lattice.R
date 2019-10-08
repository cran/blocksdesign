#' @title Mutually orthogonal latin square arrays of zize v
#'
#' @description
#' Constructs a lattice treatment design for internal use by the \code{blocksdesign} algorithm
#' Use \code{blocksdesign::MOLS} for stand-alone construction of MOLS.
#'
#' @details
#'  Returns a data frame for V*v treatments allocating treatments to blocks assuming a lattice design with r
#'  complete replicate blocks and blocks of size v. Returns a simple lattice for r=2 or a triple lattice for 
#'  r=3 for any size of v. Returns a lattice for any r < v + 2 if v is prime or any prime power where 
#'  primes=c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97) 
#'  powers=c(12,7,5,4,3,3,3,3,rep(2,17))
#'  Returns a triple lattice for v=10 and r=4 
#'  
#' @param v the dimension of the required MOLS
#' 
#' @param u equals r-2 where r is the replication number
#'
#' @return
#' \item{Treatments}{A table showing the replication number of each treatment in the design.}
#'
#' @export
  lattices=function(v,u) {
    PP=isPrimePower(v)
    p=PP$base
    q=PP$power
    primes=c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97) 
    powers=c(12,7,5,4,3,3,3,3,rep(2,17))
    ppLat = ( p%in%primes)
    if (ppLat) qmax=powers[which(primes == p)] else qmax=0
    
    if (ppLat & q>1 & q<=qmax & u>1)   {
      mols=MOLS(p,q,u)    #  lattice designs for prime powers where q > 1 and r > 3  
    } else if ( (isPrime(v) & u<v) | u < 2) {
      # cyclic designs for any prime v with r < v+2 or any v with r<4 
      mols=sapply(0:u,function(z){ sapply(0:(v-1), function(j){ (rep(0:(v-1))*z +j)%%v}) })
      mols=data.frame(row=rep(1:v,v),mols+1) 
    } else if (v==10 & u<3 & u>0) {
      #  lattice design for v = 10 and r =3 or r= 4 
      square1=       c(1,    8,    9,    2,    0,    5,    3,    6,    4,     7,
                       9,    1,    0,    4,    2,    7,    8,    5,    3,     6,
                       0,    3,    1,    6,    8,    9,    7,    4,    2,     5,
                       3,    7,    4,    1,    5,    2,    6,    9,    8,     0,
                       8,    9,    5,    0,    1,    6,    4,    2,    7,     3,
                       2,    6,    3,    7,    4,    1,    5,    8,    0,     9,
                       5,    2,    6,    3,    7,    4,    1,    0,    9,     8,
                       4,    0,    7,    5,    3,    8,    9,    1,    6,     2,
                       7,    5,    8,    9,    6,    0,    2,    3,    1,     4,
                       6,    4,    2,    8,    9,    3,    0,    7,    5,     1)
      square2 =      c(1,    3,    5,    4,    2,    6,    7,    8,    9,     0,
                       3,    4,    8,    6,    7,    9,    1,    5,    0,     2,
                       5,    6,    7,    0,    9,    1,    4,    3,    2,     8,
                       9,    8,    1,    2,    3,    0,    5,    6,    7,     4,
                       2,    0,    4,    1,    6,    7,    8,    9,    5,     3,
                       8,    1,    2,    3,    0,    5,    9,    4,    6,     7,
                       0,    5,    9,    8,    1,    2,    3,    7,    4,     6,
                       4,    9,    6,    7,    5,    8,    2,    0,    3,     1,
                       7,    2,    0,    9,    4,    3,    6,    1,    8,     5,
                       6,    7,    3,    5,    8,    4,    0,    2,    1,     9)
      
      mols=data.frame(row=rep(1:v,v),col=rep(1:v,each=v),square1+1)
      if (u==2) mols=cbind(mols,square2+1)
    } else stop("Lattice design unavailable for the current choice of model parameters")
    TF=factor(sapply(1:(u+2),function(i){seq_len(v*v)[order(as.numeric(mols[[i]]))]}))
    df=data.frame(Reps=rep(sample(1:(u+2)),each=(v*v)),Blocks=rep(sample(1:(v*(u+2))),each=v),Plots=sample(1:(v*v*(u+2))),TF)
    TF=df[ do.call(order, df), ][,4]
    return(TF)
  }