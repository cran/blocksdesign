#' @name blocksdesign-package
#' @title Blocks design package
#' @aliases blocksdesign
#' @docType package
#' 
#' @description The \code{blocksdesign} package provides functionality for the construction
#'  of block designs for general linear model treatment designs.
#'  
#' @details
#' 
#' Block designs group experimental units into homogeneous blocks to provide maximum 
#' precision of estimation of treatment effects within blocks. 
#' The most basic type of block design is a complete randomized blocks design where each 
#' block contains one or more complete replicate
#' sets of treatments. Complete randomized blocks designs estimate all treatment effects 
#' fully within individual blocks and are usually the 
#' best choice for small experiments. However, for large experiments, the 
#' variability within complete blocks can be large and then it may be beneficial to sub-divide 
#' each complete block into smaller, more homogeneous, incomplete blocks. 
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
#' usually called column blocks. Double blocking systems are commonly used for the control of 
#' block effects in two dimensions simultaneously.
#' 
#' Recursive nesting builds stucture into the block design at the design stage and an appropriate set of 
#' block effects for the control of positional effects can be selected at the analusis stage by the use
#' of selection criteria such as the AIC statistic. See Burnham and Anderson (2002).  
#' 
#' Funtionality
#' 
#' \code{blocksdesign} has two main functions:
#' 
#' i) \code{\link[blocksdesign]{blocks}}: This is a simple recursive function for nested block designs for unstructured treatment sets.
#'  The function generates designs for treatments with arbitrary levels of replication and with arbitrary depth of nesting 
#'  where each successive set of blocks is optimized within the levels of each preceding set of blocks using conditional D-optimality. 
#'  The input requires the number of blocks for each level of nesting and the algorithm automatically finds block sizes that are 
#'  as equal as possible for each level of nesting.
#'  Special block designs including square and rectangular lattice designs (see Cochran and Cox 1957) are constructed algebraically. 
#'  The outputs from the \code{blocks} function include a data frame showing the allocation of treatments to blocks for each plot of the design and a 
#'  table showing the achieved D- and A-efficiency factors for each set of nested blocks together with A-efficiency upper bounds, where available. 
#'  A plan showing the allocation of treatments to blocks in the bottom level of the design is also included in the output.
#' 
#' i) \code{\link[blocksdesign]{design}}: This is a general purpose function for arbitrary linear treatment designs and arbitrary linear block designs with 
#'  qualitative factor levels. The function first finds a D-optimal or near D-optimal treatment design 
#'  and then finds a D-optimal or near D-optimal block design for that treatment design. 
#'  The blocks \code{design} algorithm builds the blocks design by sequentially adding
#'  \code{blocks} factors where each block factor is optimized conditional on all previously added block factors.
#'  The outputs include a data frame of the block and treatment factors for each plot and a table showing the achieved D-efficiency 
#'  factors for each set of nested or crossed blocks. Fractional factorial efficiency factors based on
#'  the generalized variance of the complete factorial design are also shown.
#'  
#'  For more details see the 'blocksdesign" vignette:
#'  \code{vignette(package = "blocksdesign")}  
#' 
#' @references
#' 
#' Bates D., Maechler M., Bolker B., Walker S. (2015). 
#' Fitting Linear Mixed-Effects Models Using lme4. Journal of Statistical Software, 
#' 67(1), 1-48. doi:10.18637/jss.v067.i01.
#' 
#' Burnham, K. P., Anderson, D. R. (2002). Model Selection and Multimodel Inference:  
#' A Practical Information-Theoretic Approach (2nd ed.), Springer-Verlag.
#' 
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' 
NULL


