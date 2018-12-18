#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction
#'  of nested or crossed block designs for general linear model treatment designs.
#'  
#' @details
#' 
#' Block designs group experimental units into homogeneous blocks to provide maximum 
#' precision of estimation of treatment effects within blocks. 
#' The most basic type of block design is a complete randomized blocks design where each 
#' block contains one or more complete replicate
#' sets of treatments. Complete randomized block designs estimate all treatment effects 
#' fully within individual blocks and are usually the 
#' best choice for small experiments. However, for large experiments, the 
#' variability within complete blocks can be large and then it may be beneficial to sub-divide 
#' each complete block into smaller more homogeneous incomplete blocks. 
#' 
#' Block designs with a single level of nesting are widely used in practical 
#' research but sometimes for very large experiments a single set
#' of nested blocks may still be too large to give good control of intra-block variability. 
#' In this situation, a second set of incomplete blocks can be
#' nested within the first set to reduce the intra-block variability still further. This process of 
#' recursive nesting can be repeated as often as required until the bottom set of blocks is sufficiently
#' small to give adequate control of intra-block variability.
#'  
#' Sometimes it can be advantageous to use a double blocking system in which one set of blocks, 
#' usually called row blocks, is crossed with a second set of blocks, 
#' usually called column blocks. Double blocking systems can be valuable for controlling 
#' block effects in two dimensions simultaneously.
#' 
#' The \code{blocksdesign} package provides functionality for the construction of general 
#' multi-level block designs with nested or crossed blocks
#' for any feasible depth of nesting. The design algorithm proceeds recursively with each nested set of 
#' blocks optimized conditionally within the levels of each preceding
#' set of blocks. The analysis of incomplete block designs is complex but the availability 
#' of modern computers and modern software, for example the R mixed model
#' software package \code{lme4} (Bates et. al. 2014), makes the analysis of any feasible 
#' nested block designs with any depth of nesting practicable. 
#' 
#' The \code{blocksdesign} package has two design functions:
#' 
#' i) \code{\link[blocksdesign]{blocks}}: This is a simple recursive function for nested block designs for unstructured treatment sets.
#'  The function generates designs for treatments with arbitrary levels of replication and with arbitrary depth of nesting 
#'  where each successive set of blocks is optimized within the levels of each preceding set of blocks using conditional D-optimality. 
#'  Special block designs such as lattice designs or latin or Trojan square designs are constructed algebraically. 
#'  The outputs from the \code{blocks} function include a data frame showing the allocation of treatments to blocks for each plot of the design and a 
#'  table showing the achieved D- and A-efficiency factors for each set of nested blocks together with A-efficiency upper bounds, where available. 
#'  A plan showing the allocation of treatments to blocks in the bottom level of the design is also included in the output.
#' 
#' i) \code{\link[blocksdesign]{design}}: This is a general purpose function  for unstructured or general qualitative or quantitative
#'  factorial treatment sets. The function first finds a D-optimal or near D-optimal treatment design of the required size, possibly a simple 
#'  unstructured treatment set. The function then finds a D-optimal or near D-optimal block design for that treatment design based on a set of 
#'  defined block factors, if present. The blocks \code{design} algorithm builds the blocks design by sequentially adding
#'  \code{blocks} factors where each block factor is optimized conditional on all previous block factors. Sequential optimization allows the
#'   blocking factors to be fitted in order of importance with the largest and most important blocks fitted first and the smaller and less important
#'   blocks fitted subsequently. If there are no defined bock factors, the algorithm assumes a completely randomised treatment design.
#'  The outputs include a data frame of the block and treatment factors for each plot and a table showing the achieved D-efficiency 
#'  factors for each set of nested or crossed blocks. Fractional factorial efficiency factors based on
#'  the generalized variance of the complete factorial design are also shown (see the \code{design} documentation for more details) 
#'
#' 
#' @references
#' 
#' BATES D., MAECHLER M., BOLKER B., WALKER S. (2015). 
#' Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 
#' 67(1), 1-48. doi:10.18637/jss.v067.i01.
#' 
NULL


