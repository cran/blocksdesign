#' @title Rectangular lattice designs
#' @description
#' Internal function constructs rectangular lattice designs for resolvable nested block designs with 
#' v*(v-1) treatments, r complete replicate blocks and nested blocks of size v. 
#' 
#' Returns a rectangular lattice for v*(v-1) treatments and r replicates whenever a square lattice with v*v  
#' treatments and r+1 replicates exists; see \code{\link[blocksdesign]{squarelattice}} for square lattice designs.
#' 
#' See Cochran and Cox, Experimental Designs, 2nd Edition, Page 417 (Shrikhande method).
#'
#'  
#' @param v is the required block size and must be the integer square root of the required number of treatments
#' @param r is the required number of replicates
#' @references Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
#' @keywords internal
#'
# **************************************************************************************************
# Tests for and constructs rectangular lattice designs in top 2 levels of a rectangular lattice.
# Further nested levels for levels 3... etc can be nested within the levels of the rectangular lattice blocks
# ***************************************************************************************************
rectlattice=function(v,r) {
  LT=squarelattice(v,r+1)
  if (is.null(LT)) return(NULL)
  LT=split(LT, factor(rep(1:((r+1)*v), each=v))) 
  drop=factor((v*(v-1)+1) : (v*v))
  dropblock=which(sapply (1:length(LT), function(i) all(drop%in%LT[[i]]))   )
  droprep=(dropblock-1)%/%v + 1
  omitblocks=((droprep-1)*v + 1):(droprep*v)
  LT[omitblocks]=NULL
  LT=unlist(LT)
  TF=data.frame(droplevels(LT[!LT%in%drop]))
  names(TF)[1] = "treatments"
  rownames(TF)=NULL
  return(TF[,1])
}
