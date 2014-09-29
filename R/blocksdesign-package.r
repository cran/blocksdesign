#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction of block designs for unstructured 
#' treatment sets with arbitrary levels of replication and arbitrary depth of nesting.
#'  
#' @details
#' 
#' Block designs group experimental units into homogeneous blocks to provide maximum precision for treatment comparisons. 
#' The most basic type of block design are complete randomised blocks where each block contains one or more complete sets of treatments. 
#' Complete randomized blocks are excellent for small designs but for larger designs, the variability within blocks may become too large 
#' for reliable treatment comparison and then it may be desirable to sub-divide the complete blocks into smaller incomplete blocks.
#'
#' Traditionally, nested block designs in large experiments have used a single stratum of nested blocks contained within a set of complete main blocks.
#' The complexity of design and analysis of designs with more than a single stratum of nesting has made multi-stratum nesting infeasible for 
#' practical experiments. However, modern software such as the \code{lme4} 
#' mixed model  package (Bates et al 2014) and the availability of modern design algorithms have largely eliminated these restrictions
#' and multi-stratum nesting for large block designs is now entirely feasible. 
#' 
#' The advantage of multi-stratum nesting is that random variability can be captured across
#' a range of block sizes and this allows for more realistic modelling of block effects compared with single stratum nesting. 
#' The \code{blocksdesign} package is a general purpose tool that provides for the construction of general block designs where
#' treatments can have any number of levels of replication and blocks can be nested to any feasible depth of nesting. 
#' Where designs have one or more levels of nesting, blocks are optimized hierarchically with each successive set of nested blocks
#' optimized within the blocks of the preceding set. 
#'
#' The main function is \code{\link[blocksdesign]{blocks}} which is used to generate the actual required design. The output from 
#' \code{blocks} includes a data frame of the block and treatment factors for each plot, 
#' a data frame of the allocation of treatments to plots for each block in the design,
#' blocks-by-treatments incidence matrices for each stratum in the design and an A-efficiency factor for each stratum in the design,
#'  together with an efficiency upper bound, where available.
#'  
#' The secondary function \code{\link[blocksdesign]{efficiencies}} takes the design output from the \code{blocks} 
#' function and uses it to construct tables of efficiency factors for each pairwise treatment difference in each stratum, as required.
#' 
#' The subsidiary function \code{\link[blocksdesign]{upper_bounds}} estimates
#' A-efficiency upper bounds for regular block designs with equally replicated treatments and equal block sizes. 
#'
#' Further discussion of multi-stratum nesting can be found in the package vignette at: vignette("blocksdesign")
#' 
#' @references
#' 
#' Bates, D., Maechler, M., Bolker, B. and Walker, S. (2014). lme4: Linear mixed-effects models using Eigen and S4. 
#' R package version 1.1-6. http://CRAN.R-project.org/package=lme4
#' 
NULL


