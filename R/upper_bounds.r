#' @title Efficiency bounds 
#' 
#' @description
#' Finds upper A-efficiency bounds for regular block designs.
#' 
#' @details
#' Upper bounds for the A-efficiency of regular nested block designs 
#' (see Chapter 2.8 of John and Williams 1995). Non-trivial bounds
#' are calculated for regular block designs with equal block sizes 
#' and equal replication. All other designs return NA.    
#' 
#' @param nplots the total number of plots in the design.
#' 
#' @param ntrts the total number of treatments in the design.
#' 
#' @param nblocks the total number of blocks in the design.
#' 
#' @references
#' 
#' John, J. A. and Williams, E. R. (1995). Cyclic and Computer Generated Designs. Chapman and Hall, London.
#' 
#' @examples 
#' 
#' # 50 plots, 10 treatments and 10 blocks for a design with 5 replicates and blocks of size 5 
#' upper_bounds(nplots=50,ntrts=10,nblocks=10)
#'
#' @export
#'  
upper_bounds=function(nplots,ntrts,nblocks){
  if (nplots%%ntrts != 0 | nplots%%nblocks != 0 | ntrts == 1 | nblocks == 1 | (ntrts+nblocks-1)>nplots ) return(NA) 
  if (nplots%%(nblocks*ntrts) == 0 ) return(1)  
  nreps = nplots/ntrts #replication
  bsize = nplots/nblocks #block size	
  # this bound is for non-binary designs where bsize>ntrts and can be improved - see John and Williams page 44
  if (bsize > ntrts) return(round( 1 - (bsize%%ntrts)*(ntrts - bsize%%ntrts)/(bsize*bsize*(ntrts - 1)) , 5))			
  # binary designs with bsize<=ntrts
  dual=ntrts>nblocks
  if (dual)
  {
    temp = nblocks
    nblocks = ntrts
    ntrts = temp
    nreps = nplots/ntrts
    bsize = nplots/nblocks
  }	
  ebar =  ntrts*(bsize - 1)/(bsize*(ntrts - 1))
  lambda = nreps*(bsize - 1)/(ntrts - 1)
  if (isTRUE(all.equal(lambda,floor(lambda)))) 
    bound=ebar
  else
  {
    alpha = lambda - floor(lambda) # fractional part of lambda
    s2=ntrts*(ntrts-1)*alpha*(1-alpha)/((nreps*bsize)**2)
    s=sqrt(s2/((ntrts-1)*(ntrts-2)))
    if ( alpha < ntrts/(2*(ntrts - 1)) )
      z = alpha*((ntrts + 1)*alpha - 3)
    else
      z = (1 - alpha)*(ntrts - (ntrts + 1)*alpha)
    
    s31 = alpha*ntrts*(ntrts - 1)*z/((nreps*bsize)**3)
    if (floor(lambda)==0)  # integer part of lambda
      s32 = alpha*ntrts*(ntrts - 1)*(  (ntrts + 1)*alpha*alpha - 3*alpha - bsize + 2)/((nreps*bsize)**3)
    else
      s32=s31
    
    U1= ebar - (ntrts - 2)*s*s/(ebar + (ntrts - 3)*s)
    U2= ebar - (1 - ebar)*s2/((1 - ebar)*(ntrts - 1) - s2)
    U3= ebar - s2*s2/((ntrts - 1)*(s31+ ebar*s2))	
    U4= ebar - s2*s2/((ntrts - 1)*(s32+ ebar*s2)) 	
    bound=min(U1,U2,U3,U4,na.rm = TRUE)	
  }	
  if (dual) {
    temp = nblocks
    nblocks = ntrts
    ntrts = temp
    bound = (ntrts - 1)/((ntrts - nblocks) + (nblocks - 1)/bound)
  }			
  round(bound,5)
}	
