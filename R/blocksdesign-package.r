#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction
#'  of block and treatment designs for general linear models.
#'  
#' @details
#' 
#' Randomized complete block designs are usually the designs of choice for small experiments but for
#' large experiments, sub-division of complete blocks into smaller incomplete blocks 
#' can provide improved precision on within-block treatment comparisons. Incomplete block designs with a single
#' level of nesting are widely used but for large experiments with many plots, a single level of nesting
#' may not be adequate for the control of variability over a range of scales of measurement. Multi-level 
#' block designs with a hierarchy of nested block sizes can be valuable for accommodating block variability 
#' over the full range of scales of measurement and the \code{blocksdesign} algorithm provides for 
#' the recursive nesting of block factors down to any feasible depth of nesting. 
#' 
#' Sometimes a factorial blocking system
#' with two or more sets of crossed factorial block factors is needed for the control of
#' variability in two or more dimensions simultaneously and the \code{blocksdesign} algorithm can also
#' be used to build crossed block designs by sequential addition of block factors in a defined order of fitting. 
#' 
#' Functionality
#' 
#' \code{blocksdesign} has two main functions:
#' 
#' i) \code{\link[blocksdesign]{blocks}}: This is a simple recursive function for nested block designs for unstructured treatment sets.
#'  The function generates designs for treatments with arbitrary levels of replication and with arbitrary depth of nesting 
#'  where each successive set of blocks is optimized within the levels of each preceding set of blocks using conditional D-optimality. 
#'  The input requires the number of blocks for each level of nesting and the algorithm automatically finds block sizes that are 
#'  as equal as possible for each level of nesting.
#'  Special block designs including square and rectangular lattice designs (see Cochran and Cox 1957) are constructed algebraically from MOLS
#'  or Graeco-Latin squares. The outputs from the \code{blocks} function include a data frame showing the allocation of treatments
#'  to blocks for each plot of the design and a table showing the achieved D- and A-efficiency factors for each set of nested blocks
#'  together with A-efficiency upper bounds, where available. A plan showing the allocation of treatments to blocks in the bottom
#'   level of the design is also included in the output.
#' 
#' i) \code{\link[blocksdesign]{design}}: This is a general purpose function for linear models for qualitative
#'  and quantitative level treatment factors and for qualitative level block factors. 
#'  The function first finds a D-optimal or near D-optimal treatment design 
#'  and then finds a D-optimal or near D-optimal block design conditional on that choice of treatment design. 
#'  The blocks \code{design} algorithm builds the blocks design by sequentially adding
#'  \code{blocks} factors where each block factor is optimized conditional on all previously added block factors remaining constant.
#'  The outputs include a data frame of the block and treatment factors for each plot and a table showing the achieved D-efficiency 
#'  factors for each set of nested or crossed blocks. Fractional factorial efficiency factors based on
#'  the generalized variance of the complete factorial design are also shown.
#'  
#'  For further explanation see Edmondson (2020) and \code{vignette(package = "blocksdesign")}.
#' 
#' @references
#' 
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' Edmondson, R.N. Multi-level Block Designs for Comparative Experiments. JABES (2020). 
#' https://link.springer.com/article/10.1007/s13253-020-00416-0 
#'
NULL


