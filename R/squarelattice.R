#' @title Square lattice designs
#' @description
#' Internal function constructs square lattice designs for resolvable nested block designs with 
#' v*v treatments, r complete replicate blocks and nested blocks of size v. 
#' 
#' Returns a simple lattice with r = 2 or a triple lattice with
#' r = 3 for any size of v. 
#' 
#' Returns a quadruple lattice with r = 4 for any v <= 30.
#' 
#' Returns a lattice for any r < v + 2 if v is a prime or prime power with p^q less than or equal to:
#' 
#'\itemize{
#' \item{2**12} 
#' \item{3**7} 
#' \item{5**5} 
#' \item{7**4} 
#' \item{(11,13,17,19)**3}  
#' \item{(23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97)**2}
#' }
#' 
#' @param v is the required block size and is the integer square root of the required number of treatments
#' @param r is is the required number of replicates
#' @keywords internal
#'
squarelattice=function(v,r) {
  mols=NULL
  Z=isPrimePower(v)
  if (!isFALSE(Z) & r>3 ) {
    mols=MOLS(Z$base,Z$power,r-2)
  } else if (r==4) { # 4 replicate Lattices 
    mols=GraecoLatin(v)
  } else if (r<4) { # simple Lattices (Latin squares for non-prime v)
    z=0:(v-1)
    s1=sapply(1:v,function(i) {(z+i-1)%%v + 1})
    mols=data.frame(Row=rep(1:v,each=v),Col=rep(1:v,v),T1=as.numeric(s1))
  }
  if (!is.data.frame(mols)) return(NULL) 
  TF=factor(sapply(1:r,function(i){seq_len(v*v)[order(as.numeric(mols[[i]]))]}))
  TF=data.frame(Reps=rep(sample(1:r),each=(v*v)),Blocks=rep(sample(1:(v*r)),each=v),plots=sample(1:(v*v*r)),TF)
  TF=TF[ do.call(order, TF), ]
  names(TF)[4] = "treatments"
  rownames(TF)=NULL
  return(TF[,4])
}