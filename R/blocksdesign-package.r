#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction of nested or crossed block designs for factorial or 
#' unstructured treatment sets with arbitrary levels of replication and arbitrary depth of nesting.
#'  
#' @details
#' 
#' Block designs aim to group experimental units into homogeneous blocks to provide maximum precision of estimation of treatment effects within blocks. 
#' The most basic type of block design is a complete randomized blocks design where every block contains one or more complete replicate
#' sets of treatments. Complete randomized block designs estimate all treatment effects fully within individual blocks and are usually the 
#' best choice for small experiments. However, for large experiments, the average variability within complete replicate 
#' blocks can be large and then it may be beneficial to sub-divide each complete replicate block into smaller incomplete blocks which give 
#' improved precision of comparison on inter-block treatment effects.
#' 
#' Block designs with a single level of nesting are widely used in practical research but sometimes the blocks of large designs with a single set
#' of nested blocks may still be too large to give good control of intra-block variability. In this situation, a second set of incomplete blocks can be
#' nested within the first set to reduce the intra-block variability still further. This process of recursive nesting of blocks can be repeated indefinitely 
#' as often as required until the bottom set of blocks are sufficiently small to give good control of intra-block variability.
#'  
#' Sometimes it can also abe advantageous to use a double blocking system in which one set of blocks, usually called row blocks, is crossed with a second set of blocks, usually called column blocks. 
#' Double blocking systems can be valuable for controlling block effects in two dimensions simultaneously.
#' 
#' The \code{blocksdesign} package provides functionality for the construction of general block designs with simple or crossed blocks that can be nested repeatedly 
#' to any feasible depth of nesting. The design algorithm proceeds recursively with each
#' nested set of blocks optimized conditionally within each preceding set of blocks. Block sizes within 
#' any nested level are as equal as possible and never differ by more than a single plot. The analysis of incomplete block designs is complex but 
#' the availability of modern computers and modern software, for example the R mixed model software package \code{lme4} (Bates et al 2014), 
#' makes the analysis of any feasible nested block designs with any depth of nesting practicable. 
#' 
#'Currently, the \code{blocksdesign} package has two main block design functions:
#' 
#' i) \code{\link[blocksdesign]{blocks}} : The \code{blocks} function is used to generate block designs for any arbitrary number of unstructured treatments where each treatment can have
#'  any arbitrary number of replicates. This function generates arbitrary nested block designs with arbitrary depth of nesting 
#'  where each succesive set of blocks is optimized within the levels of each preceding set of blocks using a conditional D-optimality design criterion. 
#'  Special block designs such as lattice designs or latin or Trojan square designs are constructed algebraically. 
#'  The outputs from the \code{blocks} function includes a data frame showing the allocation of treatments to blocks for each plot of the design and a table showing
#'   the achieved D- and A-efficiency factors for each set of nested blocks together with A-efficiency upper bounds, where available. 
#'   A plan showing the allocation of treatments to blocks in the bottom level of the design is also included in the output.
#' 
#' ii) \code{\link[blocksdesign]{factblocks}} : The \code{factblocks} function is used to generate designs for factorial treatment sets.
#'   The \code{factblocks} function finds a D-optimal or near D-optimal treatment design
#'   of the required size for any required factorial model and then finds a D-optimal or near D-optimal block design
#'   for the fitted treatment design using the same algorithm as in the \code{blocks} function.
#'   The output from \code{factblocks} includes a data frame of the block and treatment factors for 
#'  each plot and a table showing the achieved D-efficiency factors for each set of nested blocks. Fractional factorial efficiency factors based on
#'  the generalized variance of the complete factorial design are also shown (see the \code{factblocks} documentation for details) 
#'
#' Further discussion of designs with repeatedly nested strata can be found in the package vignette at: vignette("blocksdesign")
#' 
#' 
#' @references
#' 
#' Bates, D., Maechler, M., Bolker, B. and Walker, S. (2014). lme4: Linear mixed-effects models using 'Eigen' and S4. 
#' R package version 1.1-12. https://CRAN.R-project.org/package=lme4
#' 
NULL


