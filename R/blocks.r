#' @title Block designs for unstructured treatment sets
#'
#' @description
#'
#' Constructs randomized multiply nested block designs for unstructured treatment sets.
#'
#' @details
#'
#' Constructs randomized nested block designs for any arbitrary
#' number of unstructured treatments and any arbitrary feasible depth of nesting.
#' 
#' \code{treatments} is a set of numbers that partitions the total number of treatments into
#' equi-replicate treatment sets.
#' 
#' \code{replicates} is a set of treatment replication numbers that defines the replication number for each equi-replicate
#' treatment set.
#' 
#' \code{blocks} defines the levels of a sequence of nested blocks in decreasing order of block size.
#' Each level defines the number of blocks nested within the blocks of the preceding level
#' where the top-level block is assumed to be single super-block containing the full set of plots. 
#' The algorithm finds block sizes automatically for each level of nesting and the block sizes within
#' any one level of nesting will never differ in size by more than a single plot.
#'
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number 
#' of unreplicated treatments using the \code{treatments} and \code{replicates} formula. 
#' Usually, however, it will be preferable to find an efficient blocked design
#' for the replicated treatment sets and then to add unreplicated treatments individually, by hand. 
#'
#' Incomplete block designs are constructed algorithmically except for certain designs for r replicates
#' of \code{v x v} treatments or r-1 replicates of \code{v x (v-1)} treatments in blocks of size v.
#' Provided that a set of r-1 mutually orthogonal Latin squares of size \code{v x v} exists, these designs
#' are constructed algebraically and are guaranteed to achieve optimality. See 
#' \code{\link[blocksdesign]{squarelattice}} and \code{\link[blocksdesign]{rectlattice}} for information about
#' which design sizes are constructed algebraically.
#'
#' Outputs:
#'
#' \itemize{
#' \item  A data frame allocating treatments to blocks with successive nested strata in standard block order.\cr
#' \item  A table showing the replication number of each treatment in the design. \cr
#' \item  A table showing block levels and the achieved D-efficiency and A-efficiency factor for each nested level
#'    together with A-efficiency upper bounds, where available. \cr
#' \item  A plan showing the allocation of treatments to blocks in the bottom level of the design.\cr
#' }
#' @param treatments  the required number of treatments partitioned into equally replicated treatment sets.
#' @param replicates  the replication number for each partitioned treatment set.
#' @param blocks the number of nested blocks in each level of nesting from the top level down.
#' @param seed an integer initializing the random number generator.
#' @param searches the maximum number of local optima searched for a design optimization. 
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima.
#' @return
#' \item{Treatments}{A table showing the replication number of each treatment in the design.}
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
#' \donttest{blocks(64,4,list(4,8,2,2))}
#' 
#' # 100 treatments x 4 replicates with 4 main blocks nested blocks of size 10 (lattice square)
#' blocks(100,4,list(4,10)) 
#' 
#' @export
#' @importFrom stats coef anova lm model.matrix as.formula setNames 
#'
 blocks = function(treatments,replicates,blocks=NULL,searches=NULL,seed=NULL,jumps=1) {
   options(contrasts=c('contr.treatment','contr.poly'))
   options(warn=0)
   tol = .Machine$double.eps ^ 0.5
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
    stop("treatments and replicates parameters must both be the same length")
  hcf=HCF(replicates) 
  if (is.null(blocks)) blocks=1
  if (anyNA(blocks)|any(is.nan(blocks))|any(!is.finite(blocks))|any(blocks%%1!=0)|any(blocks<1)) 
    stop(" blocks invalid")
  blocks=blocks[blocks!=1]
  if(length(blocks)==0) blocks=1
  ntrts=sum(treatments)
  TF=unlist(lapply(1:hcf,function(i){sample(rep(factor(1:ntrts), rep(replicates/hcf,treatments)))}))
  if (is.null(searches)) searches=1+5000%/%length(TF)
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  if (is.null(names(blocks))) 
    names(blocks)=unlist(lapply(1:length(blocks), function(j) {paste0("Level_",j-1)}))
  if (prod(blocks)*2>length(TF)) stop("Too many blocks for the available plots  - blocks must contain at least two plots")
  if (  (sum(treatments)+prod(blocks))>length(TF)) stop("Too many parameters for the available plots")
  # Finds nested block sizes where all block sizes are as equal as possible within each level of nesting
  blocksizes=length(TF)
  for (i in 1:length(blocks))
    blocksizes=unlist(lapply(1:length(blocksizes), function(j) {
      nestblocksizes=rep(blocksizes[j]%/%blocks[i],blocks[i])
      blocksizes=nestblocksizes+rep(c(1,0),c(blocksizes[j]-sum(nestblocksizes), blocks[i]-blocksizes[j]+sum(nestblocksizes)))
    }))
  blocksGrid=function(blocks){expand.grid(lapply(length(blocks):1,function(i) {
    factor(seq(blocks[i]),labels=lapply(1:blocks[i], function(j){paste0("Blocks_",j)}))}))[length(blocks):1]}
  blkdesign=blocksGrid(blocks)
  # the Null column is needed for 'restriction' of the first (main) set of blocks 
  blkdesign=cbind("Null"=factor(rep("Blocks_1",nrow(blkdesign))),blkdesign)
  blkDesign=blkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  blkDesign=data.frame(lapply(1:ncol(blkDesign),function(i) { droplevels(interaction(blkDesign[,1:i], lex.order = TRUE))}))
  colnames(blkDesign)=c("Null",labels(blocks))
  Z=buildblocks(TF,blkDesign,searches,seed,jumps)
  
  list(Treatments=Z$Treatments,Blocks_model=Z$Blocks_model,Design=Z$Design,Plan=Z$Plan,seed=seed,searches=searches,jumps=jumps)
 }
 