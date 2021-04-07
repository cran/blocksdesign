#' @title Block designs for unstructured treatment sets
#'
#' @description
#'
#' Constructs randomized multi-level nested block designs for unstructured treatment sets.
#'
#' @details
#'
#' Constructs randomized multi-level nested block designs for any arbitrary
#' number of unstructured treatments and any arbitrary feasible depth of nesting.
#' 
#' \code{treatments} is a partition of the number of treatments into equi-replicate treatment sets.
#' 
#' \code{replicates} is a set of replication numbers for the equi-replicate treatment sets.
#' 
#' \code{blocks} are the nested blocks levels in decreasing order of block size where
#' each level defines the number of blocks nested within the blocks of the preceding level.
#' The top-level block is assumed to be single super-block containing a full set of plots. 
#' The algorithm finds block sizes automatically for each level of nesting and the block sizes within
#' each level of nesting will never differ by more than a single plot.
#'
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number 
#' of unreplicated treatments using the \code{treatments} and \code{replicates} formula. 
#' However, it may be preferable to find an efficient blocked design
#' for the replicated treatment sets and then to add the unreplicated treatments heuristically. 
#'
#' The \code{blocks} function constructs all block designs algorithmically except for certain special block designs
#' for r replicates of \code{v x v} treatments or r-1 replicates of \code{v x (v-1)} treatments in blocks of size v.
#' Provided that a set of r-1 mutually orthogonal Latin squares of size \code{v x v} exists, these designs are
#' constructed algebraically and are guaranteed to achieve optimality. See \code{\link[blocksdesign]{squarelattice}} 
#' and \code{\link[blocksdesign]{rectlattice}} for information about which design sizes are constructed algebraically.
#'
#' 
#' @param treatments  the total required number of treatments partitioned into equally replicated treatment sets.
#' @param replicates  the replication numbers of the equally replicated treatment sets.
#' @param blocks the number of nested blocks in each level of nesting from the top level down.
#' @param seed an integer initializing the random number generator.
#' @param searches the maximum number of local optima searched for a design optimization. 
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima.
#' @return
#' \item{Replication}{A table showing the replication number of each treatment in the design.}
#' \item{Design}{Data frame giving the optimized block and treatment design in plot order.}
#' \item{Plan}{Data frame showing a plan view of the treatment design in the bottom level of the design.}
#' \item{Blocks_model}{The D-efficiencies and the A-efficiencies of the blocks in each nested level of the 
#'  design together with A-efficiency upper-bounds, where available.}
#' \item{seed}{Numerical seed used for random number generator.}
#' \item{searches}{Maximum number of searches used for each level.}
#' \item{jumps}{Number of random treatment swaps used to escape a local maxima.}
#' @references
#' Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
#' @examples
#' 
#' ## The number of searches in the following examples have been limited for fast execution.  
#' ## In practice, the number of searches may need to be increased for optimum results.
#' ## Designs should be rebuilt several times to check that a near-optimum design has been found.  
#' 
#' # Completely randomized design for 6 treatments with 2 replicates and 1 control with 4 replicates 
#' blocks(treatments=list(6,1),replicates=list(2,4))
#' 
#' # 12 treatments x 4 replicates in 4 complete blocks with 4 sub-blocks of size 3
#' # rectangular lattice see Plan 10.10 Cochran and Cox 1957.
#' \donttest{blocks(treatments=12,replicates=4,blocks=list(4,4))}
#'
#' # 3 treatments x 2 replicates + 2 treatments x 4 replicates in two complete randomized blocks
#' blocks(treatments=list(3,2),replicates=list(2,4),blocks=2)
#'
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block
#' blocks(treatments=50,replicates=4,blocks=list(4,5))
#'
#' # as above but with 20 additional single replicate treatments, one single treatment per sub-block
#' \donttest{blocks(treatments=list(50,20),replicates=list(4,1),blocks=list(4,5))}
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,blocks=4)
#'
#' # 128 treatments x 2 replicates with two main blocks and 3 levels of nesting
#'  \donttest{blocks(128,2,list(2,2,2,2))}
#' 
#' # 64 treatments x 4 replicates with 4 main blocks, 8 nested sub-blocks of size 8
#' # (lattice), 16 nested sub-sub blocks of size 4 and 32 nested sub-sub-sub blocks of size 2
#'  \donttest{blocks(64,4,list(4,8,2,2))}
#' 
#' # 100 treatments x 4 replicates with 4 main blocks nested blocks of size 10 (lattice square)
#' blocks(100,4,list(4,10)) 
#' 
#' @export
#' @importFrom stats coef anova lm model.matrix as.formula setNames 
#' @importFrom plyr count
#'
 blocks = function(treatments,replicates,blocks=NULL,searches=NULL,seed=NULL,jumps=1) {
   options(contrasts=c('contr.treatment','contr.poly'))
   options(warn=0)
   tol = .Machine$double.eps ^ 0.5
   
   # ***********************************************************************************************
   # finds n sub-blocks nested within one main block where all sub-block sizes are as equal as possible
   # *********************************************************************************************** 
     subB=function(mainSize,n) {
       subBlocks=rep(mainSize%/%n,n) 
       if (mainSize%%n>0) subBlocks[1:(mainSize%%n)]=subBlocks[1:(mainSize%%n)]+1
       return(subBlocks)
     }
     
  # ***********************************************************************************************
  #  finds a set of equal or near-equal sub-blocks nested within a set of main blocks
  # each sub-set must contain the same number b of sub-blocks 
  # sub-blocks will never differ by more than a single plot in size
  # *********************************************************************************************** 
     nestedSizes=function(mainBlocks,b) {
       subBlocks=unlist(lapply(1:length(mainBlocks), function(j) {subB(mainBlocks[j],b)}))
       return(subBlocks)
     }
  # ***********************************************************************************************  
    
  if (is.list(treatments)) treatments=unlist(treatments)
  if (is.list(blocks)) blocks=unlist(blocks)
  if (is.list(replicates)) replicates=unlist(replicates)
  if (missing(treatments)|is.null(treatments)) stop(" Treatments missing or not defined ")
  if (missing(replicates)|is.null(replicates)) stop(" Replicates missing or not defined ")
  if (!is.null(seed)) set.seed(seed)
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>25) 
    stop(" maximum number of jumps is 25 ")
  if (any(replicates%%1!=0)|any(replicates<1)) stop(" replication numbers must be integers")
  if (anyNA(treatments)|any(is.nan(treatments))|any(!is.finite(treatments))|any(treatments<1)) 
    stop(" treatments parameter invalid")
  if (length(replicates)!=length(treatments)) 
    stop("the treatments parameter and the replicates parameter must be of equal length so
    that each treatment set has a matching replication number")
  if (is.null(blocks)) blocks=1
  if (anyNA(blocks)|any(is.nan(blocks))|any(!is.finite(blocks))|any(blocks%%1!=0)|any(blocks<1)) 
    stop(" blocks invalid")
  blocks=blocks[blocks!=1]
  if (length(blocks)==0) blocks=1
  nplots=sum(treatments*replicates)
  if (is.null(searches)) searches=1+5000%/%nplots
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  if (is.null(names(blocks))) 
    names(blocks)=unlist(lapply(1:length(blocks), function(j) {paste0("Level_",j)}))
  if (prod(blocks)*2 > nplots ) stop("Too many blocks for the available plots  - blocks must contain at least two plots")
  if ((sum(treatments)+prod(blocks)-1) > nplots) stop("Too many parameters for the available plots")
  # Finds nested block sizes where all block sizes are as equal as possible within each level of nesting
  blocksizes=nplots
  for (i in 1:length(blocks))
    blocksizes=nestedSizes(blocksizes,as.numeric(blocks[i]))
  blocksGrid=function(blocks){expand.grid(lapply(length(blocks):1,function(i) {
    factor(seq(blocks[i]),labels=lapply(1:blocks[i], function(j){paste0("B",j)}))}))[length(blocks):1]}
  blkDesign=blocksGrid(blocks)[rep(1:length(blocksizes),blocksizes),,drop=FALSE]

  blkDesign=data.frame(lapply(1:ncol(blkDesign),function(i) { droplevels(interaction(blkDesign[,1:i], lex.order = TRUE))}))
  colnames(blkDesign)=labels(blocks)
 
  hcf=HCF(replicates) 
  TF=data.frame(Treatments=factor(unlist(lapply(1:hcf,function(i){sample(rep(1:sum(treatments),rep(replicates/hcf,treatments)))}))))
  Z=nestedBlocks(TF[,1,drop=FALSE],blkDesign,searches,seed,jumps) # nestedBlocks function optimizes the blocks design
  list(Replication=count(TF),Blocks_model=Z$Effic,Design=Z$Design,Plan=Z$Plan,seed=seed,searches=searches,jumps=jumps)
 }
 