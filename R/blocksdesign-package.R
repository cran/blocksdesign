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
#' Randomized complete blocks are the designs of choice for small experiments with few treatments.
#'  For large experiments with many treatments, however, a single set of complete blocks may not be adequate  
#' and then sub-division into smaller nested blocks may be required. Block designs with a single level of nesting are
#' widely used but a single level of nesting may be inadequate for very large experiments with many treatments.
#' \code{blocksdesign} provides for the construction of designs with multiple levels of nesting down to any
#' feasible depth of nesting. 
#' 
#' Sometimes block designs for the control of variability in two or more dimensions are required and
#' \code{blocksdesign} can also build crossed block designs allowing for both additive and interactive crossed
#' block effects simultaneously.
#' 
#' The \code{blocksdesign} package has two main functions:
#' 
#' i) \code{\link[blocksdesign]{blocks}}: This is a simple recursive function for nested blocks for
#'  unstructured treatments. The function generates designs for treatments with arbitrary levels of replication
#'  and with arbitrary depth of nesting where blocks sizes are assumed to be as equal as possible for each level of nesting.
#'  Special square and rectangular lattice designs (see Cochran and Cox 1957) are constructed
#'  algebraically from mutually orthogonal Latin squares (MOLS). The outputs from the \code{blocks} function include a data
#'  frame showing the allocation of treatments to blocks and a table showing the achieved D- and A-efficiency factors for each
#'  set of nested blocks together with A-efficiency upper bounds, where available. A plan showing the allocation of treatments
#'  to blocks for the bottom level of the design is also included in the output.
#' 
#' ii) \code{\link[blocksdesign]{design}}: This is a general purpose function for linear models with qualitative
#'  or quantitative level treatment factors and qualitative level block factors. The function finds a D-optimal
#'  or near D-optimal design for a specified treatment model and then finds a conditional D-optimal or 
#'  near D-optimal block design for that choice of treatment design. The \code{design} algorithm builds the blocks design
#'  by sequentially adding \code{blocks} factors where each blocks factor is optimized conditional on all previously 
#'  added \code{blocks} factors. The outputs include a data frame of the block and treatment factors for each plot and a table
#'  showing the achieved D-efficiency factors for each set of nested or crossed blocks. 
#'  Fractional factorial efficiency factors based on the generalized variance of the complete factorial design are also shown.
#'  
#'  Other available functions are \code{\link[blocksdesign]{A_bound}}, which finds upper A-efficiency bounds for regular
#'  block designs, \code{\link[blocksdesign]{MOLS}}, which constructs sets of mutually orthogonal prime-power 
#'  Latin squares (MOLS), \code{\link[blocksdesign]{GraecoLatin}}, which constructs mutually orthogonal Graeco-Latin 
#'  squares not necessarily prime-power, \code{\link[blocksdesign]{isPrime}}, which tests an integer for primality, 
#'  \code{\link[blocksdesign]{isPrimePower}}, which factorizes prime powers and \code{\link[blocksdesign]{HCF}},
#'  which finds the highest common factor (hcf) for a set of positive integer numbers.  
#'  
#'  For further explanation see Edmondson (2020) and \code{vignette(package = "blocksdesign")}.
#' 
#' @references
#' 
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' Edmondson, R.N. Multi-level Block Designs for Comparative Experiments. JABES 25, 500–522 (2020).
#' 
NULL


