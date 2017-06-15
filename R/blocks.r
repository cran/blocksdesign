#' @title Block designs
#'
#' @description
#'
#' Constructs randomized nested block designs for factorial or fractional factorial or unstructured treatment designs with any
#' feasible depth of nesting and up to two crossed block structures in each level of nesting.
#'
#' @details
#'
#' Constructs randomized nested block designs with arbitrary depth of nesting for factorial or fractional
#' factorial  or unstructured treatment designs. The treatment model can be any arbitrary combination of quantitative
#' or qualitative factorial model terms or can be a single set of unstructured treatments. 
#'  
#' The \code{treatments} parameter defines the treatment factors of the design and can be either a data frame with
#' a column for each factor and a row for each factor combination (see examples) or a set of P cardinal numbers for an unstructured treatment set
#' where each cardinal represents a set of equally replicated treatments and the sum of the cardinals is the total number of treatments
#' (see examples). 
#' 
#' If the \code{treatments} parameter is a data frame, the treatment factors can be any mixture of qualitative or quantitative level factors and the 
#' treatment model can be any feasible model defined by the \code{models} formula of the \code{\link[stats]{model.matrix}} package (see examples). 
#' 
#' Quantitative factors can be modelled either by raw or by orthogonal polynomials. Orthogonal polynomials are numerically more stable
#' than raw polynomials and are usually the best choice at the design stage. Polynomial models can be fitted at the analysis stage either by raw or
#' by orthogonal polynomials regardless of the type of polynomial fitted at the design stage.
#'   
#' The \code{replicates} parameter defines the required replication for the treatments design.  
#' If the \code{treatments} parameter is a data frame, the replication number must be a single number, not necessarily integral, representing any
#' required multiple or any required fraction of the \code{treatments} data frame. If the \code{treatments} parameter is a set of cardinal numbers, the 
#' \code{replicates} parameter must be a matching set of replication numbers, one for each equally replicated treatments sets (see examples).
#' 
#' If the \code{treatments} parameter is a data frame and if the replication number is non-integral, the algorithm will find a
#' D-optimal or near D-optimal fraction of the required size for the fractional part of replication number, assuming the required design is non-singular.
#' 
#' The \code{rows} parameter, if any, defines the nested row blocks for each level of nesting taken in order from the highest to the lowest. The
#' first number, if any, is the number of nested row blocks in the first-level of nesting, the second number, if any, is the number of nested row blocks in
#' the second-level of nesting and so on for the required feasible depth of nesting.
#' 
#' The \code{columns} parameter, if any, defines the nested column blocks for each level of nesting taken in order from the highest to the lowest.
#' The first number, if any, is the number of nested column blocks in the first-level of nesting, the second, if any, is the number of nested column blocks in
#' the second-level of nesting and so on for the required feasible depth of nesting.
#' 
#' The \code{rows} and \code{columns} parameters, if defined, must be of equal length and if a simple nested blocks design is required for
#' any particular level of nesting, the number of columns for that level should be set to unity. Any required combination of simple or
#' crossed blocks can be obtained by appropriate choice of the levels of the \code{rows} and \code{columns} parameters.
#'  
#' If the \code{columns} parameter is undefined, a single crossed block is assumed for each level of nesting. 
#'
#' If both the \code{rows} parameter and the \code{columns} parameter are null, the default block design will be a set of orthogonal
#' main blocks equal in number to the highest common factor of the replication numbers. If the \code{rows} parameter is defined but the \code{columns}
#' parameter is null, the design will have simple nested blocks for each level of nesting, as defined by the \code{rows} parameter.
#'
#' Block sizes are always as nearly equal as possible and will never differ by more than a single plot for any particular block classification. 
#' Row blocks and column blocks must always contain at least two plots per block and this restriction will constrain the permitted numbers of 
#' rows and columns in the various nested levels of a block design.
#'
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number of single unreplicated treatments 
#' to give augmented blocks that never differ in size by more than a single plot. General crossed block designs are more complex and currently 
#' the algorithm will only accommodate single unreplicated treatments in a crossed block design if the block sizes of the replicated part of 
#' the design are all equal in each nested level of the design.
#'
#' For any particular level of nesting, the algorithm first optimizes the row blocks conditional on any higher-level blocks
#' and then optimizes the columns blocks, if any, conditional on the rows blocks.
#' 
#' The efficiency factor for a fractional factorial design is the generalized variance of the complete factorial design divided by the generalized variance of
#' the fractional factorial design where the generalized variance of a design is the (1/p)th power of the determinant of the crossed-product of the p-dimensional
#' model matrix divided by the number of observations in the design. 
#'
#' Special designs:
#'
#' Trojan designs are row-and-column designs for p replicates of v*p treatments arranged in p-rows and p-columns where v < p and 
#' where every row x column intersection contains v plots. Trojan designs have orthogonal rows and columns and optimal rows x columns
#' blocks and exist whenever p is prime or prime-power. The \code{blocksdesign} constructs these designs algebraically from mutually
#' orthogonal Latin squares (MOLS).
#'
#' Square lattice designs are efficient resolvable incomplete block designs for r replicates of p*p treatments arranged in blocks of size p where
#' r < p+2 for prime or prime power p or r < 4 for general p. \code{blocksdesign} constructs lattice designs algebraically from Latin squares or MOLS.
#'
#' Lattice designs and Trojan designs based on prime-power MOLS require the \code{\link[crossdes]{MOLS}} package.
#'
#' All other designs are constructed algorithmically.
#'
#' Comment:
#'
#' Row-and-column designs may contain useful treatment information in the individual row-by-column intersection blocks but \code{blocksdesign} does not currently
#' optimize the efficiency of these blocks except for the special case of Trojan designs.
#'
#' Row-and-column design with 2 complete treatment replicates, 2 complete rows and 2 complete columns will always confound one treatment contrast in the
#' rows-by-columns interaction. For these designs, it is impossible to nest a non-singular block design in the rows-by-columns intersections and instead
#' we suggest a randomized nested blocks design with four incomplete main blocks.
#'
#' Outputs:
#'
#' The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
#'  \item  A table showing the replication number of each treatment in the design. \cr
#'  \item  An efficiency factor for fractional factorial treatment designs. \cr
#'  \item  A table showing the block levels and the achieved D-efficiency and A-efficiency (unstructured treatment designs only) factors for each stratum together
#'   with A-efficiency upper bounds, where available. \cr
#'  \item  A plan showing the allocation of treatments to blocks or to rows and to columns in the bottom stratum of the design (unstructured treatment
#'  designs only).\cr
#' }
#'
#' @param treatments  either a data frame with columns for individual treatment factors or a partition
#' of the total required number of treatments into sets of equally replicated treatments.
#'
#' @param replicates  either a single replication number, not necessarily integral, if the \code{treatments} parameter is a data frame, or a set of 
#' replication numbers, one per replication set, if the \code{treatments} parameter is a partition.
#'
#' @param rows the numbers of rows nested in each higher-level block for each level of nesting from the top block downwards. The top-level block is a
#' single super-block which does not need to be defined unless a completely randomized design is required. The default number of blocks is the hcf
#' of the replication numbers, which gives a maximal set of orthogonal row blocks.
#' 
#' @param columns the numbers of columns nested in each higher-level block for each level of nesting from the top block downwards. The \code{rows} and 
#' \code{columns} parameters must be of equal length unless the \code{columns} parameter is null, in which case the columns block design degenerates to a
#' single column block for each level of nesting, which gives a simple nested row blocks design.
#'
#' @param model a model equation for the treatment factors in the design where the equation is defined using the model.matrix notation
#' in the {\link[stats]{model.matrix}} package. If undefined, the model is a full factorial treatment design.
#'
#' @param seed  an integer initializing the random number generator. The default is a random seed.
#'
#' @param searches  the maximum number of local optima searched for a design optimization. The default is 1 plus the floor of 10000 divided by the number of plots.
#'
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#'
#' @return
#' \item{Treatments}{The treatment factors defined by the \code{treatments} inputs in standard factorial order.}
#' \item{model.matrix}{The model.matrix used to define the \code{treatments} design.}
#' \item{Design}{Data frame giving the optimized block and treatment factors in plot order.}
#' \item{Plan}{Data frame for single factor designs showing a plan view of the treatment design in the bottom stratum of the design. A NULL plan is returned for multi-factor designs.}
#' \item{BlocksEfficiency}{The D-efficiencies and the A-efficiencies (unstructured designs) of the blocks in each stratum of the design together with A-efficiency upper-bounds, where available.}
#' \item{DesignEfficiency}{The generalized variance of the complete factorial design divided by the generalized variance of the fractional factorial design.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#'
#'
#' @references
#'
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. https://CRAN.R-project.org/package=crossdes
#'
#' Edmondson R. N. (1998). Trojan square and incomplete Trojan square designs for crop research. Journal of Agricultural Science, Cambridge, 131, pp.135-142
#'
#' Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
#'
#'
#' @examples
#' 
#' ## The number of searches in the following examples have been limited for fast execution.  
#' ## In practice, the number of searches may need to be increased for optimum results.
#' ## Designs should be rebuilt several times to check that a near-optimum design has been found.  
#' 
#' 
#' ## Factorial designs defined by a treatments data frame and a factorial model equation.
#' 
#' # Main effects of five 2-level factors in a half-fraction of a 4 x 4 row-and column design.
#' GF = expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' blocks(treatments=GF,model="~ F1+F2+F3+F4+F5",replicates=.5,rows=4,columns=4,searches=20)
#' 
#' # Quadratic regression for one 6-level numeric factor in 2 randomized blocks
#' blocks(treatments=expand.grid(X=1:6),model=" ~ poly(X,2)",rows=2,searches=5) 
#' 
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' blocks(treatments=GF,model=" ~ (F1+F2+F3+F4+F5)*(F1+F2+F3+F4+F5)",rows=4,searches=5)
#' 
#' # First-order model for 1/3rd fraction of four qualitative 3-level factors in 3  blocks
#' GF=expand.grid(F1=factor(1:3),F2=factor(1:3),F3=factor(1:3),F4=factor(1:3))
#' blocks(treatments=GF,model=" ~ (F1+F2+F3+F4)",replicates=(1/3),rows=3,searches=5)
#' 
#' # Second-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' GF=expand.grid( F1=factor(1:3), F2=factor(1:3), F3=factor(1:3), F4=factor(1:3), F5=factor(1:3) )
#' modelform=" ~ (F1+F2+F3+F4+F5)*(F1+F2+F3+F4+F5)"
#' blocks(treatments=GF,model=modelform,rows=3,replicates=(1/3),searches=1)
#' 
#' # Second-order model for two qualitative and two quantitative level factors in 4 randomized blocks
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:3),V1=1:3,V2=1:4)
#' modelform=" ~ F1 + F2 + poly(V1,2) +  poly(V2,2) + (poly(V1,1)+F1+F2):(poly(V2,1)+F1+F2) "
#' blocks(treatments=GF,model=modelform,rows=4,searches=5)
#' 
#' # Plackett and Burman design for eleven 2-level factors in 12 runs (needs large number of searches)
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2),
#' F6=factor(1:2),F7=factor(1:2),F8=factor(1:2),F9=factor(1:2),F10=factor(1:2),F11=factor(1:2))
#' \dontrun{blocks(treatments=GF,model="~ F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11",replicates=(12/2048))}
#' 
#'
#' ## Unstructured treatments partitioned into equally replicated treatment sets
#'
#' # 3 treatments x 2 replicates + 2 treatments x 4 replicates 
#' blocks(treatments=c(3,2),replicates=c(2,4),searches=10)
#'
#' # 4 treatments x 4 replicates with 2 main rows each containing two complete replicates
#' blocks(treatments=4,replicates=4,rows=2)
#'
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block
#' blocks(treatments=50,replicates=4,rows=c(4,5))
#'
#' # as above but with 20 single replicate treatments giving one extra treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),rows=c(4,5))
#'
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,rows=4)
#'
#' # 4 replicates of 13 treatments arranged in a 13 x 4 Youden rectangle
#' blocks(treatments=13,replicates=4,rows=13,columns=4)
#'
#' # 64 treatments x 2 replicates with nested 8 x 8 row-and-column designs in two main blocks
#' blocks(treatments=64,replicates=2,rows=c(2,8),columns=c(1,8),searches=10)
#'
#' # 64 treatments x 2 replicates with two main blocks and a 4 x 4 row-and-column in each main block
#' blocks(treatments=64,replicates=2,rows=c(2,4),columns=c(1,4),searches=10)
#' 
#' # 128 treatments x 2 replicates with two main blocks and 3 levels of nesting
#' \dontrun{blocks(128,2,c(2,2,2,2))}
#'
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames
#'
blocks = function(treatments,replicates=1,rows=NULL,columns=NULL,model=NULL,searches=NULL,seed=sample(10000,1),jumps=1) {

  # ********************************************************************************************************************************************************
  # Finds the highest common factor (hcf) of a set of numbers omitting any zero values (Euclidean algorithm)
  # ********************************************************************************************************************************************************
  HCF=function(replicates)  {
    if (any(ceiling(replicates)!=floor(replicates))) return(1)
    replicates=sort(replicates[replicates>0])
    for (i in  seq_len(length(replicates)))
      while (!isTRUE(all.equal(replicates[i]%%replicates[1],0))) replicates[c(1,i)] = c(replicates[i]%%replicates[1], replicates[1])
    return(replicates[1])
  }
  # ********************************************************************************************************************************************************
  # Tests a given number for primality and returns TRUE or FALSE
  # ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v < 4) return(TRUE)
    if ( isTRUE(all.equal(v %% 2,0)) |  isTRUE(all.equal(v %% 3,0)) ) return(FALSE)
    if (v<25) return(TRUE)
    for(i in  6*seq_len(length(floor((sqrt(v)+1)/6)))        )
      if ( isTRUE(all.equal(v %% (i-1) , 0)) |   isTRUE(all.equal(v %% (i+1) , 0)) ) return(FALSE)
    return(TRUE)
  }
 
  # ********************************************************************************************************************************************************
  # Finds row and column sizes in each stratum of a design
  # ********************************************************************************************************************************************************
  Sizes=function(blocksizes,stratum) {
    nblocks=length(blocksizes)
    newblocksizes=NULL
    for (j in 1:nblocks) {
      rowsizes=rep(blocksizes[j]%/%rows[stratum],rows[stratum])
      resid=blocksizes[j]-sum(rowsizes)
      if (resid>0)
        rowsizes[1:resid]=rowsizes[1:resid]+1
      rowcolsizes=vector(mode = "list", length =rows[stratum])
      for ( z in 1:rows[stratum])
        rowcolsizes[[z]]=rep(rowsizes[z]%/%columns[stratum] , columns[stratum])
      shift=0
      for (z in seq_len(rows[stratum])) {
        resid=rowsizes[z]-sum(rowcolsizes[[z]])
        if (resid>0) {
          rowcolsizes[[z]][(shift:(shift+resid-1))%%columns[stratum]+1]=rowcolsizes[[z]][(shift:(shift+resid-1))%%columns[stratum]+1]+1
          shift=shift+resid
        }
      }
      newblocksizes=c(newblocksizes,unlist(rowcolsizes))
    }
    return(newblocksizes)
  }
  # ********************************************************************************************************************************************************
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  FactEstEffics=function(TF,MF,BF) {
    if (nlevels(MF)==nlevels(BF)) return(1)
    dd=data.frame(TF,BF)
    TM=model.matrix(as.formula(model),dd)[,-1,drop=FALSE] # drops mean contrast
    BM=model.matrix(~BF,dd)[,-1,drop=FALSE] # drops mean contrast
    TM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    BM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(BM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    BM=BM[, -seq( nlevels(BF)/nlevels(MF), nlevels(BF) , by=nlevels(BF)/nlevels(MF) ) ,drop=FALSE]
    TB=crossprod(TM,BM)
    RI=backsolve(  chol(crossprod(TM)) ,diag(ncol(TM)))
    QI=backsolve(chol(crossprod(BM)),diag(ncol(BM)))
    U=crossprod(t(crossprod(RI,TB)),QI)
    return(round(det( diag(ncol(TM))-tcrossprod(U))**(1/ncol(TM)),6))
  }
  # ********************************************************************************************************************************************************
  # Finds efficiency factors for row-and-column designs
  # ********************************************************************************************************************************************************
  FactRowColEffics=function(Design) {
    TF=Design[, c( (ncol(Design)-ncol(treatments)+1):ncol(Design)),drop=FALSE]
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata))
      for (j in 1:3)
        effics=c(effics,FactEstEffics(TF,Design[,3*(i-1)+1],Design[,3*(i-1)+1+j]))
    names =unlist(lapply(1:strata, function(j) {c(paste("Rows",j),paste("Columns",j),paste("Rows x Columns",j))}))
    blocklevs=unlist(lapply(1:strata, function(j) {c( nlevels(Design[,3*(j-1)+2]),nlevels(Design[,3*(j-1)+3]),nlevels(Design[,3*(j-1)+4]))}))
    efficiencies=data.frame(cbind(names,blocklevs,effics))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies")
    return(efficiencies)
  }
  # ********************************************************************************************************************************************************
  # Calculates D and A-efficiency factors for treatment factors TF assuming block factor BF
  # ********************************************************************************************************************************************************
  FactBlocksEffics=function(Design) {
    TF=Design[, c( (ncol(Design)-ncol(treatments)+1):ncol(Design)),drop=FALSE]
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata))
      effics=c(effics,FactEstEffics(TF,Design[,i],Design[,i+1]))
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocklevs=unlist(lapply(2:(strata+1), function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocklevs,effics))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies")
    return(efficiencies)
  }
  # ********************************************************************************************************************************************************
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  EstEffics=function(TF,BF) {
    k=nlevels(BF)
    if (k==1) return(c(1,1))
    if (nlevels(TF)<=k)
      e=eigen( (diag(nlevels(TF))-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(nlevels(TF)-1)] else
        e=c(rep(1,(nlevels(TF)-k)),
            eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])
      return(round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6))
  }
  # ********************************************************************************************************************************************************
  # Finds efficiency factors for block designs
  # ********************************************************************************************************************************************************
  BlockEfficiencies=function(Design) {
    effics=matrix(NA,nrow=strata,ncol=2)
    for (i in seq_len(strata))
      effics[i,]=EstEffics(Design[,ncol(Design)],Design[,i])
    bounds=rep(NA,strata)
    if (regReps)
      for (i in seq_len(strata))
        if (nunits%%nlevels(Design[,i])==0 )
          bounds[i]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,i]) )
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocklevs=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocklevs,effics,bounds))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    return(efficiencies)
  }
  # ********************************************************************************************************************************************************
  # displays efficiency factors for Latin squares
  # ********************************************************************************************************************************************************
  LatinEfficiencies=function(Design) {
    names = c(paste("Rows",1),paste("Columns",1))
    blocklevs=c( nlevels(Design[,1]),nlevels(Design[,2]))
    effics=c(1,1)
    bounds=c(1,1)
    efficiencies=data.frame(cbind(names,blocklevs,effics,effics,bounds))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    return(efficiencies)
  }
  # ********************************************************************************************************************************************************
  # Finds efficiency factors for row-and-column designs
  # ********************************************************************************************************************************************************
  RowColEfficiencies=function(Design) {
    Design=Design[,-(ncol(Design)-1)]
    effics=matrix(NA,nrow=(3*strata),ncol=2)
    for (i in 1:strata) {
      effics[3*(i-1)+1,]=EstEffics(Design[,ncol(Design)],Design[,3*(i-1)+1])
      effics[3*(i-1)+2,]=EstEffics(Design[,ncol(Design)],Design[,3*(i-1)+2])
      effics[3*(i-1)+3,]=EstEffics(Design[,ncol(Design)],Design[,3*(i-1)+3])
    }
    bounds=rep(NA,(3*strata))
    if (max(replicates)==min(replicates)) {
      for (i in seq_len(strata))  {
        if (nunits%%nlevels(Design[,3*(i-1)+1])==0)
          bounds[3*(i-1)+1]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+1]))
        if (nunits%%nlevels(Design[,3*(i-1)+2])==0)
          bounds[3*(i-1)+2]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+2]))
        if (nunits%%nlevels(Design[,3*(i-1)+3])==0)
          bounds[3*(i-1)+3]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+3]))
      }
    }
    names =unlist(lapply(1:strata, function(j) {c(paste("Rows",j),paste("Columns",j),paste("Rows x Columns",j))}))
    blocklevs=unlist(lapply(1:strata, function(j) {c( nlevels(Design[,3*(j-1)+1]),nlevels(Design[,3*(j-1)+2]),nlevels(Design[,3*(j-1)+3]))}))
    efficiencies=data.frame(cbind(names,blocklevs,effics,bounds))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    return(efficiencies)
  }
  # *******************************************************************************************************************************************************
  # Returns v-1 cyclically generated v x v Latin squares plus the rows array (first) and the columns array (last). If v is prime, the squares are MOLS
  # ********************************************************************************************************************************************************
  cMOLS=function(v) {
    mols=lapply(0:(v-1),function(z){do.call(rbind, lapply(0:(v-1), function(j){ (rep(0:(v-1))*z +j)%%v} ))})
    mols[[v+1]]=t(mols[[1]])
    mols=mols[c(1,length(mols),2:(length(mols)-1))]
    mols=lapply(1:length(mols),function(z){ (z-1)*v+mols[[z]] })
    return(mols)
  }
  # *******************************************************************************************************************************************************
  # Tests for and constructs  balanced lattice designs
  # ********************************************************************************************************************************************************
  lattice=function(v,r) {
    TF=vector(length=r*v*v)
    if ( r<4 | (isPrime(v) & r<(v+2)) ) {
      TF=rep(seq_len(v*v),r)[order(unlist(cMOLS(v))[1:(r*v*v)])]
    } else if (r<(v+2)  & (v*v)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      TF[1:(2*v*v)]=c(seq_len(v*v),seq_len(v*v)[order(rep(0:(v-1),v))]   )
      if (r>2) {
        index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(v*v))
        mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])
        for (i in 1:(r-2)  )
          TF[ c( (v*v*(i-1)+1) : (v*v*i)) + 2*v*v  ] = seq_len(v*v)[order(as.numeric(mols[,,i]))]
      }
    } else if (v==10  & r==4) {
      square1=c(1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7,
                5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8)
      square2=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1,
                6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6)
      TF=c(seq_len(100),seq_len(100)[order(rep(0:9,10))],seq_len(100)[order(square1)],seq_len(100)[order(square2)])
    } else TF=NULL
    if (!is.null(TF)) TF=data.frame(factor(TF))
    return(TF)
  }
  # *******************************************************************************************************************************************************
  # Tests for balanced trojan designs and constructs available designs
  # ********************************************************************************************************************************************************
  trojan=function(r,k) {
    TF=vector(length=r*r*k)
    if (isPrime(r)) {
      for (z in 1:k)
        for (y in 0:(r-1))
          for (x in 0:(r-1))
            TF[(x + y*r)*k + z]=(y+x*z)%%r + (z-1)*r +1
    } else if ((r*r)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(r*r))
      mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])
      for (i in 1:r)
        for (j in 1:r)
          for (x in 1:k)
            TF[x+(j-1)*k+(i-1)*k*r]=mols[i,j,x]+(x-1)*r
    } else TF=NULL
    if (!is.null(TF)) TF=data.frame(factor(TF))
    return(TF)
  }
     # ******************************************************************************************************************************************************** 
     # Maximises the blocks design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
     # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
     # ********************************************************************************************************************************************************
     gDMax=function (VTT,VBB,VTB,TF,MF,TM,BM,restrict)  {
       locrelD=1
       mainSizes=tabulate(restrict)
       nSamp=pmin( rep(8,nlevels(restrict)), mainSizes)
       repeat {
         kmax=1
         for (k in 1: nlevels(restrict)) {
           s=sort(sample( seq_len(nrow(TF)) [restrict==levels(restrict)[k]], nSamp[k])) 
           TMB=crossprod(t(crossprod(t(TM[s,]),VTB)),t(BM[s,]))
           TMT=crossprod(t(crossprod(t(TM[s,]),VTT)),t(TM[s,]))
           BMB=crossprod(t(crossprod(t(BM[s,]),VBB)),t(BM[s,]))
           TMB=sweep(TMB,1,diag(TMB))
           TMT=sweep(TMT,1,diag(TMT))
           BMB=sweep(BMB,1,diag(BMB))
           dMat=(1+TMB+t(TMB))**2 - (TMT + t(TMT))*(BMB + t(BMB))
           sampn=which.max(dMat)
           i=1+(sampn-1)%%nrow(dMat)
           j=1+(sampn-1)%/%nrow(dMat)
           if (dMat[i,j]>kmax) {kmax=dMat[i,j]; pi=s[i]; pj=s[j]} 
         }
         if (kmax>(1+tol)) {
           locrelD=locrelD*kmax
           t=TM[pi,]-TM[pj,]
           b=BM[pj,]-BM[pi,]
           up=gUpDate(VTT,VBB,VTB,t,b)
           VTT=up$VTT
           VBB=up$VBB
           VTB=up$VTB
           TF[c(pi,pj),]=TF[c(pj,pi),]
           TM[c(pi,pj),]=TM[c(pj,pi),]
         } else if (sum(nSamp) == nrow(TF)) break
         else nSamp=pmin(mainSizes,2*nSamp)
       }
       list(VTT=VTT,VBB=VBB,VTB=VTB,TF=TF,locrelD=locrelD,TM=TM)
     }
     # ******************************************************************************************************************************************************** 
     # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
     # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
     # ********************************************************************************************************************************************************
     DMax=function(VTT,VBB,VTB,TF,fBF,TM,restrict) {  
       locrelD=1
       mainSizes=tabulate(restrict)
       nSamp=pmin(rep(8,nlevels(restrict)), mainSizes)
       repeat {
         kmax=1
         for (k in 1: nlevels(restrict)) {
           s=sort(sample(seq_len(length(TF))[restrict==levels(restrict)[k]], nSamp[k])) 
           TB=VTB[TF[s],fBF[s],drop=FALSE]
           TT=VTT[TF[s],TF[s],drop=FALSE]
           BB=VBB[fBF[s],fBF[s],drop=FALSE]
           TB=sweep(TB,1,diag(TB))
           TT=sweep(TT,1,diag(TT))
           BB=sweep(BB,1,diag(BB))
           dMat=(TB+t(TB)+1)**2-(TT+t(TT))*(BB+t(BB))
           sampn=which.max(dMat)
           i=1+(sampn-1)%%nrow(dMat)
           j=1+(sampn-1)%/%nrow(dMat)
           if (dMat[i,j]>kmax) {kmax=dMat[i,j]; pi=s[i]; pj=s[j]} 
         }
         if (kmax>(1+tol)) {
           locrelD=locrelD*kmax
           up=UpDate(VTT,VBB,VTB,TF[pi],TF[pj],fBF[pi],fBF[pj])
           VTT=up$VTT
           VBB=up$VBB
           VTB=up$VTB
           TF[c(pi,pj)]=TF[c(pj,pi)]
           TM[c(pi,pj),]=TM[c(pj,pi),]
         } else if (sum(nSamp) == length(TF)) break
         else nSamp=pmin(mainSizes,2*nSamp)
       }
       TF=as.data.frame(TF)
       colnames(TF)="Treatments"
       list(VTT=VTT,VBB=VBB,VTB=VTB,TF=TF,locrelD=locrelD,TM=TM)
     } 
   # ********************************************************************************************************************************************************
   # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
   # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
   # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
   # ********************************************************************************************************************************************************
   gUpDate=function(VTT,VBB,VTB,t,b) {
     VTTt=crossprod(VTT,t)
     VBBb=crossprod(VBB,b)
     VTBt=crossprod(VTB,t)
     VTBb=crossprod(t(VTB),b)
     tMt=as.numeric(crossprod(t,VTTt))
     bMb=as.numeric(crossprod(b,VBBb))
     tMb=as.numeric(crossprod(b,VTBt))
     f1=(VTTt+VTBb)/sqrt(2)
     f2=(VBBb+VTBt)/sqrt(2)
     g1=(VTBb-VTTt)/sqrt(2)
     g2=(VBBb-VTBt)/sqrt(2)
     a=(tMt+bMb+2*tMb)/2
     b=(tMt+bMb-2*tMb)/2
     c=(bMb-tMt)/2
     d=(1+a)*(1-b)+c*c
     VTT=VTT- (tcrossprod(f1)*(1-b) - tcrossprod(g1)*(1+a) + (tcrossprod(g1,f1)+tcrossprod(f1,g1))*c)/d
     VBB=VBB- (tcrossprod(f2)*(1-b) - tcrossprod(g2)*(1+a) + (tcrossprod(g2,f2)+tcrossprod(f2,g2))*c)/d
     VTB=VTB- (tcrossprod(f1,f2)*(1-b) - tcrossprod(g1,g2)*(1+a) + (tcrossprod(g1,f2)+tcrossprod(f1,g2))*c)/d
     list(VTT=VTT,VBB=VBB,VTB=VTB)
   }
  # ******************************************************************************************************************************************************** 
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb   
  # ********************************************************************************************************************************************************
  UpDate=function(VTT,VBB,VTB,ti,tj,bi,bj) { 
    mtt=VTT[ti,ti]+VTT[tj,tj]-2*VTT[tj,ti]
    mbb=VBB[bi,bi]+VBB[bj,bj]-2*VBB[bi,bj]
    mtb=1-VTB[ti,bi]+VTB[tj,bi]+VTB[ti,bj]-VTB[tj,bj]  
    TBbij=VTB[,bi]-VTB[,bj]
    TBtij=VTB[ti,]-VTB[tj,]
    TTtij=VTT[,ti]-VTT[,tj]
    BBbij=VBB[bi,]-VBB[bj,]
    Z1 = (TBbij-TTtij)/sqrt(2*mtb+mtt+mbb)   
    Z2 = (BBbij-TBtij)/sqrt(2*mtb+mtt+mbb)
    W1 = ( sqrt(2*mtb+mtt+mbb)*(TTtij+TBbij) - (mbb-mtt)*Z1) /(2*sqrt(mtb**2-mtt*mbb))
    W2 = ( sqrt(2*mtb+mtt+mbb)*(TBtij+BBbij) - (mbb-mtt)*Z2) /(2*sqrt(mtb**2-mtt*mbb))
    VTT = VTT - tcrossprod(Z1) + tcrossprod(W1)
    VBB = VBB - tcrossprod(Z2) + tcrossprod(W2)
    VTB = VTB - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(VTT=VTT,VBB=VBB,VTB=VTB)
  }  
  # ********************************************************************************************************************************************************
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,MF,fBF,restrict,VTT,VBB,VTB,TM,BM) {
    globrelD=0
    relD=1
    globTF=TF
    blocksizes=tabulate(fBF)
    if (regReps & max(blocksizes)==min(blocksizes) & ncol(treatments)==1 & is.factor(treatments[,1]))
      rowsEffBound=upper_bounds(nrow(TF),nlevels(TF[,1]),nlevels(fBF)) 
    else if ( ncol(treatments)==1 & nlevels(fBF)>1 & is.factor(treatments[,1])) 
      rowsEffBound=1
    else rowsEffBound=NA
    for (r in 1:searches) {
      if (unstructured) dmax=DMax(VTT,VBB,VTB,TF[,1],fBF,TM,restrict)
      else dmax=gDMax(VTT,VBB,VTB,TF,MF,TM,BM,restrict) 
      if (dmax$locrelD>(1+tol)) {
        relD=relD*dmax$locrelD
        TF=dmax$TF
        TM=dmax$TM
        VTT=dmax$VTT
        VBB=dmax$VBB
        VTB=dmax$VTB
        if (relD>globrelD) {
          globTF=TF 
          globrelD=relD
          if (!is.na(rowsEffBound)) {
            reff=EstEffics(globTF[,1],fBF)[2]
            if (isTRUE(all.equal(rowsEffBound,reff))) return(globTF)
          }  
        }
      }  
      if (r<searches) { 
        for (iswap in 1:jumps) {
          counter=0
          repeat {
            counter=counter+1
            s1=sample(seq_len(nunits),1)
            available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
            z= seq_len(nunits)[MF==MF[s1] & restrict==restrict[s1] & fBF!=fBF[s1] & available]   
            if (length(z)==0) next
            if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z)
            if (unstructured) {
              v11=VTT[c(TF[s[1],],TF[s[2],]),c(TF[s[1],],TF[s[2],])]
              v22=VBB[c(fBF[s[1]],fBF[s[2]]),c(fBF[s[1]],fBF[s[2]])]
              v12=VTB[c(TF[s[1],],TF[s[2],]),c(fBF[s[1]],fBF[s[2]])]
              Dswap=(1-2*sum(diag(v12))+sum(v12))**2-(2*sum(diag(v22))-sum(v22))*(2*sum(diag(v11))-sum(v11))
               } else {
              TMT=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],VTT)),TM[s[2],]-TM[s[1],])
              BMB=crossprod(t(crossprod(BM[s[1],]-BM[s[2],],VBB)),BM[s[2],]-BM[s[1],])
              TMB=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],VTB)),BM[s[2],]-BM[s[1],] )
              Dswap=(1+TMB)**2-TMT*BMB
            } 
            if (Dswap>.1| counter>1000) break
          }
          if (counter>1000) return(globTF) # finish with no non-singular swaps
          relD=relD*Dswap 
          if (unstructured) 
            up=UpDate(VTT,VBB,VTB,TF[s[1],],TF[s[2],], fBF[s[1]], fBF[s[2]])
          else
            up=gUpDate(VTT,VBB,VTB, TM[s[1],]-TM[s[2],], BM[s[2],]-BM[s[1],] )
          VTT=up$VTT
          VBB=up$VBB
          VTB=up$VTB
          TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]  
          TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
        } #end of jumps
      }
    }
    return(globTF)
  }
  # ******************************************************************************************************************************************************** 
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,restrict,pivot,rank,nunits) {
    candidates=NULL
    while (is.null(candidates)) {
      if (rank<(nunits-1)) s1=sample(pivot[(1+rank):nunits],1) else s1=pivot[nunits]
      available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
      candidates = (1:nunits)[ MF==MF[s1] & restrict==restrict[s1] & BF!=BF[s1] & available==TRUE]
    }
    if ( length(candidates)>1 ) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  }
  # ********************************************************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************
  NonSingular=function(TF,MF,BF,restrict,TM,BM) {
    fullrank=ncol(cbind(TM,BM))
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:1000) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      s=Swaps(TF,MF,BF,restrict,pivot,rank,nunits)
      tindex=1:nrow(TF)
      tindex[c(s[1],s[2])]=tindex[c(s[2],s[1])]
      Q=qr(t(cbind(BM,TM[tindex,])))
      if (Q$rank>rank) {
        TF=TF[tindex,,drop=FALSE]
        TM=TM[tindex,,drop=FALSE]
        rank=Q$rank
        pivot=Q$pivot
      }
    }
    stop(" Unable to find an initial non-singular choice of treatment design for this choice of block design")
  }
  # *******************************************************************************************************************************************************
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************
  blocksOpt=function(TF,MF,fBF,BF,restrict) {
      TM=model.matrix(as.formula(model),TF)[,-1,drop=FALSE] # drops mean contrast
      TM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(TM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)})) #centres treatments within main blocks
      if (nlevels(MF)==1) BM=model.matrix(as.formula(~BF))[,-1,drop=FALSE] else BM=model.matrix(as.formula(~MF+MF:BF))[,-c(1:nlevels(MF)),drop=FALSE]
      BM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(BM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)})) #centres sub-blocks within main blocks
      BM=BM[ ,as.numeric(matrix(1:(nlevels(BF)*nlevels(MF)),nrow=nlevels(BF),ncol=nlevels(MF),byrow=TRUE)[-nlevels(BF),]) ,drop=FALSE] # reorder within main blocks
      if ( (ncol(cbind(TM,BM))+1) > nrow(TM)  ) stop( paste("Too many parameters: plots df = ",nrow(TM)-1," model:df = ",ncol(cbind(TM,BM)) )) 
      nonsing=NonSingular(TF,MF,fBF,restrict,TM,BM)
      TF=nonsing$TF
      TM=nonsing$TM
      V=chol2inv(chol(crossprod(cbind(TM,BM))))
      VTT = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
      VBB = V[ (ncol(TM)+1):ncol(V) , (ncol(TM)+1):ncol(V),drop=FALSE]
      VTB = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
      if (unstructured) {
        VTT=rbind ( cbind( VTT, rep(0,ncol(TM) )) , rep(0,(1+ncol(TM)) ))
        VBB=rbind(  cbind( VBB, matrix(0, nrow=ncol(BM),ncol=nlevels(MF)) ),  matrix(0, nrow=nlevels(MF),ncol=(ncol(BM)+nlevels(MF))))
        VTB=rbind(  cbind( VTB, matrix(0, nrow=ncol(TM),ncol=nlevels(MF)) ),  rep(0,(ncol(BM)+nlevels(MF)  )))
        reorder=c(rbind( matrix(seq_len(ncol(BM)), nrow=(nlevels(BF)-1), ncol=nlevels(MF)), seq(ncol(BM)+1, nlevels(fBF) )))
        VBB=VBB[reorder,reorder] 
        VTB=VTB[,reorder] 
      }
      TF=Optimise(TF,MF,fBF,restrict,VTT,VBB,VTB,TM,BM)
    return(TF)
  }
  # *******************************************************************************************************************************************************
  # Model names for complete factorial models
  # ********************************************************************************************************************************************************
  modelnames=function (treatments) {
    unlist(lapply(1:ncol(treatments), function(i) {
      if (!is.factor(treatments[,i])) paste0("poly(",colnames(treatments)[i],",",length(unique(treatments[,i]))-1,")") else colnames(treatments)[i]
    })) 
  }
  # ********************************************************************************************************************************************************
  # Updates variance matrix for swapped rows where mi is swapped out and mj is swapped in
  # ********************************************************************************************************************************************************
  fractUpDate=function(D,mi,mj) {
    f=crossprod(D,mj)
    g=crossprod(D,mi)
    a=as.numeric(crossprod(mj,f))
    b=as.numeric(crossprod(mi,g))
    c=as.numeric(crossprod(mi,f))
    W=g*sqrt(1+a)-f*c/sqrt(1+a)
    d=(1+a)*(1-b)+c*c
    V=f*sqrt(d/(1+a))
    D=D-(tcrossprod(V)-tcrossprod(W))/d
    return(D)
  }
  # ********************************************************************************************************************************************************
  # Fractional factorials
  # ********************************************************************************************************************************************************
  factorial=function(TF,model,replicates) {
    allfactors=all(unlist(lapply(TF, class))=="factor")
    munits=nrow(TF)*floor(replicates)
    nunits=munits+floor(nrow(TF)*replicates%%1)
    fTF=TF[rep(seq_len(nrow(TF)),ceiling(replicates)), ,drop=FALSE] # expand by indexing
    fTM=model.matrix(as.formula(model),fTF)
    if ( ncol(fTM)>nunits) stop("Fractional factorial design too small to estimate all the required model parameters ")
    maxDeff=(det(crossprod(fTM))**(1/ncol(fTM)))/nrow(fTM) 
    gfracDeff=0
    geff=0
    gTF=NULL
    if (is.null(searches)) searches=1+10000%/%nunits
    for (i in 1:searches) {
      counter=0
      repeat {
        fTF=fTF[unlist(lapply(1:ceiling(replicates),function(j){(j-1)*nrow(TF)+sample(1:nrow(TF))})),,drop=FALSE]
        fTM=model.matrix(as.formula(model),fTF)
        if (counter>99 |  qr(fTM[1:nunits,,drop=FALSE])$rank==ncol(fTM)) break else counter=counter+1
      }
      if (counter>99)  stop("Unable to find non-singular fractional design ")
      D=chol2inv(chol(crossprod(fTM[1:nunits,,drop=FALSE])))
      repeat {
        M1DM1=1-diag(tcrossprod(tcrossprod(fTM[(munits+1):nunits,,drop=FALSE] ,D),fTM[(munits+1):nunits,,drop=FALSE] ))
        M2DM2=1+diag(tcrossprod(tcrossprod(fTM[(nunits+1):nrow(fTM),,drop=FALSE],D),fTM[(nunits+1):nrow(fTM),,drop=FALSE]))
        Z=tcrossprod(M1DM1,M2DM2)+(tcrossprod(  tcrossprod(fTM[(munits+1):nunits,,drop=FALSE] ,D),  fTM[(nunits+1):nrow(fTM),,drop=FALSE]))**2
        maxindex=which.max(Z)
        i=1+(maxindex-1)%%nrow(Z)
        j=1+(maxindex-1)%/%nrow(Z)
        if(Z[i,j]<(1+tol)) break
        tswapout=fTF[munits+i,,drop=FALSE]
        tswapin =fTF[nunits+j,,drop=FALSE]
        fTF[munits+i,]=tswapin
        fTF[nunits+j,]=tswapout
        mswapout=fTM[munits+i,,drop=FALSE]
        mswapin =fTM[nunits+j,,drop=FALSE]
        fTM[munits+i,]=mswapin
        fTM[nunits+j,]=mswapout
        D=fractUpDate(D,as.numeric(mswapout),as.numeric(mswapin))
      }
      fracDeff=det(crossprod(fTM[(1:nunits),,drop=FALSE]))**(1/ncol(fTM))/nunits
      if (fracDeff>gfracDeff) {
        gfracDeff=fracDeff
        geff=gfracDeff/maxDeff
        gTF=fTF[(1:nunits),,drop=FALSE]
        if (geff>(1-tol) & allfactors) break
      }
    }
    return(list(TF=gTF,eff=geff,fraction=nunits/nrow(TF)))
  }
  # ******************************************************************************************************************************************************** 
  # Design data frames for rows columns and row.column blocks
  # ********************************************************************************************************************************************************     
  dataframes=function(rows,columns) {
    if (max(rows)==1 & max(columns)==1) return(list(nblkdesign=NULL,nrowdesign=NULL,ncoldesign=NULL,fblkdesign=NULL,frowdesign=NULL,fcoldesign=NULL))
    fblkdesign=as.data.frame(lapply(1:(strata+1),function(i) {gl(cumblocks[i],cumblocks[strata+1]/cumblocks[i],
                                                             labels=unlist(lapply(1:cumblocks[i], function(j) {paste0("Blocks_",j)})))}))
    nblkdesign=as.data.frame(lapply(1:(strata+1),function(i) {gl(nestblocks[i],cumblocks[strata+1]/cumblocks[i],
                                                             labels=unlist(lapply(1:cumblocks[i], function(j) {paste0("Blocks_",j)})))})) 
    nrowdesign=data.frame(lapply(1:strata,function(i) {gl(rows[i],cumblocks[strata+1]/cumblocks[i]/rows[i],cumblocks[strata+1],
                                                          labels=unlist(lapply(1:rows[i], function(j) {paste0("Rows_",j)}))) }))
    ncoldesign=data.frame(lapply(1:strata,function(i) {gl(columns[i],cumblocks[strata+1]/cumblocks[i+1],cumblocks[strata+1],
                                                          labels=unlist(lapply(1:columns[i], function(j) {paste0("Cols_",j)}))) }))
    frowdesign=data.frame(lapply(1:ncol(nrowdesign), function(i){ interaction(fblkdesign[,i], nrowdesign[,i], sep = ":", lex.order = TRUE) }))
    fcoldesign=data.frame(lapply(1:ncol(ncoldesign), function(i){ interaction(fblkdesign[,i], ncoldesign[,i], sep = ":", lex.order = TRUE) }))
    colnames(fblkdesign)=unlist(lapply(1:ncol(fblkdesign), function(j) {paste0("Blocks_",j-1)}))
    colnames(frowdesign)=unlist(lapply(1:ncol(frowdesign), function(j) {paste0("Rows_",j)}))
    colnames(fcoldesign)=unlist(lapply(1:ncol(fcoldesign), function(j) {paste0("Cols_",j)}))
    colnames(nblkdesign)=unlist(lapply(1:ncol(fblkdesign), function(j) {paste0("Blocks_",j-1)}))
    colnames(nrowdesign)=unlist(lapply(1:ncol(frowdesign), function(j) {paste0("Rows_",j)}))
    colnames(ncoldesign)=unlist(lapply(1:ncol(fcoldesign), function(j) {paste0("Cols_",j)}))
    list(nblkdesign=nblkdesign,nrowdesign=nrowdesign,ncoldesign=ncoldesign,fblkdesign=fblkdesign,frowdesign=frowdesign,fcoldesign=fcoldesign)
  }
  # ********************************************************************************************************************************************************
  # Main design program which tests input variables, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************
  options(contrasts=c('contr.SAS','contr.poly'))
  tol=.Machine$double.eps^0.5
  if (missing(treatments)|is.null(treatments)) stop(" Treatments missing or not defined ")
  if (is.null(replicates)|anyNA(replicates)|any(is.nan(replicates))|any(!is.finite(replicates))) stop(" replicates invalid")
  if (is.na(seed) | !is.finite(seed) | is.nan(seed) | seed%%1!=0 | seed<0 ) stop(" seed parameter invalid  ")
    set.seed(seed)
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) stop(" number of jumps parameter is invalid (max is 10) ")
  fractionalEff=data.frame("Full",1)
  hcf=HCF(replicates)
  # constructs treatment  data frame possibly a fractional factorial design 
  if (!is.data.frame(treatments)) {
    if (any(replicates%%1!=0)|any(replicates<1)) stop(" If the treatment numbers are defined by cardinals then the replication numbers must be defined by integers")
    if (anyNA(treatments)|any(is.nan(treatments))|any(!is.finite(treatments))|any(treatments%%1!=0)|any(treatments<1)) stop(" treatments parameter invalid")
    if (length(replicates)!=length(treatments)) stop("treatments and replicates parameters must both be the same length")
    fulltrts=treatments
    redreps=replicates[replicates>1]
    redtrts=fulltrts[replicates>1]
    treatments=data.frame(Treatments=factor(unlist(lapply(1:hcf,function(i){sample(rep(1:sum(treatments), rep(replicates/hcf,treatments)))}))))
  } else if (ceiling(replicates)==floor(replicates)) {
    ntrts=nrow(treatments)
    treatments=treatments[rep(seq_len(ntrts),replicates), ,drop=FALSE] # expand by indexing
    treatments=treatments[unlist(lapply(1:replicates,function(j){(j-1)*ntrts+sample(1:ntrts)})),,drop=FALSE] #randomize
    fractionalEff=data.frame(replicates,1)
  } else {
    Z=factorial(treatments,model,replicates)
    treatments=Z$TF
    fractionalEff=data.frame(Z$fraction,Z$eff)
  }
  nunits=nrow(treatments)
  
  colnames(fractionalEff)=c("Fraction","D-Efficiency")
  if (is.null(model)) model=paste0("~ ",paste0( unlist(lapply( 1:ncol(treatments), 
       function(i) {if (!is.factor(treatments[,i])) paste0("poly(",colnames(treatments)[i],",",length(unique(treatments[,i]))-1,")") else colnames(treatments)[i]})), collapse="*"))
  if (is.null(rows)) rows=hcf
  if (anyNA(rows)|any(is.nan(rows))|any(!is.finite(rows))|any(rows%%1!=0)|any(rows<1)|is.null(rows)) stop(" rows invalid")
  if (anyNA(columns)|any(is.nan(columns))|any(!is.finite(columns))|any(columns%%1!=0)|any(columns<1)) stop(" columns parameter invalid")
  if (is.null(columns)) columns=rep(1,length(rows))
  if (length(columns)!=length(rows)) stop("rows and columns vectors must be the same length ")
  if (max(rows*columns)==1) { rows=1; columns=1} else {index=rows*columns>1; rows=rows[index]; columns=columns[index]}
  for (i in 1:ncol(treatments))
    if (isTRUE(all.equal(treatments[,i], rep(treatments[1,i], length(treatments[,i]))))) stop("One or more treatment factors is a constant which is not valid")
  fnames=colnames(treatments)
  strata=length(rows)
  
  unstructured=(ncol(treatments)==1 & is.factor(treatments[,1]))  
 
  if (max(rows)>1 | max(columns)>1) {    
  # omit any single replicate treatments for unstructured factorial designs and find hcf for factor replicates
  if ( unstructured & min(replicates)==1 & max(replicates)>1) {
    replications=rep(replicates,fulltrts)
    newlevs=(1:sum(fulltrts))[replications>1]
    hcf=HCF(redreps)
    treatments=data.frame(Treatments=factor(unlist(lapply(1:hcf,function(i){sample(rep(1:sum(redtrts), rep(redreps/hcf,redtrts)))}))))
    levels(treatments[,1]) = newlevs
    nunits=nrow(treatments)
  }

  blocksizes=nunits
  for (i in 1:strata)
    blocksizes=Sizes(blocksizes,i)
  regBlocks=isTRUE(all.equal(max(blocksizes), min(blocksizes)))
  if ( ncol(treatments)==1 & is.factor(treatments[,1]) & max(columns)>1 & min(replicates)==1 & max(replicates)>1 & regBlocks==FALSE )
    stop("The algorithm does not deal with irregular row-and-column designs containing single replicate treatments ")
  
  if (is.null(searches)) 
    if (nunits<1000) searches=10000%/%nunits else if (nunits<5000) searches=5000%/%nunits else searches=1
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  regReps=isTRUE(all.equal(max(replicates), min(replicates)))
  # tests for viable design sizes
  
  nestblocks=c(1,rows*columns)
  cumrows=cumprod(rows)
  cumcols=cumprod(columns)
  cumblocks=c(1,cumprod(rows*columns))
  
  if (cumrows[strata]*2>nunits) stop("Too many row blocks for the available plots  - every row block must contain at least two plots")
  if (cumcols[strata]*2>nunits) stop("Too many column blocks for the available plots  - every column block must contain at least two plots")
  if (cumblocks[strata+1]>nunits & cumrows[strata]>1 & cumcols[strata]>1) stop("Too many blocks - every row-by-column intersection must contain at least one plot")
  Plots=factor(1:nunits)
 # nested factor level data frames

  df1=dataframes(rows,columns)
  nblkdesign=df1$nblkdesign
  nrowdesign=df1$nrowdesign
  ncoldesign=df1$ncoldesign
  nblkDesign=nblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  nrowDesign=nrowdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  ncolDesign=ncoldesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE] 
  fblkdesign=df1$fblkdesign
  frowdesign=df1$frowdesign
  fcoldesign=df1$fcoldesign
  fblkDesign=fblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  frowDesign=frowdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  fcolDesign=fcoldesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE] 

  # for s orthogonal Latin squares of dimension r x r there are r x kr Trojan designs for r replicates of kr treatments in blocks of size k where k<=s
  v=sqrt(nlevels(treatments[,1]))  # dimension of a lattice square
  k=nunits/cumblocks[strata+1]  # average block size
  orthoMain=(regReps & identical(replicates[1],rows[1]))
  Lattice=(regReps & regBlocks & orthoMain & max(columns)==1 & identical(v,floor(v)) & identical(k,v) & length(rows)==2)
  Trojan=(regReps & regBlocks & orthoMain & max(columns)>1 & identical(columns[1],replicates[1]) & length(rows)==1 & length(columns)==1 & k<replicates[1] & k>1) 
  Latin=(regReps & regBlocks & orthoMain & max(columns)>1 & identical(columns[1],replicates[1]) & length(rows)==1 & length(columns)==1 & k==1)

  TF=NULL
  r1=replicates[1]
  if (ncol(treatments)==1 & is.factor(treatments[,1])) {
    if (Lattice)
      TF=lattice(v,r1)
    else if (Trojan)
      TF=trojan(r1,k)
    else if (Latin) 
    TF=data.frame(factor(unlist(lapply(1:r1,function(i){(i-2+1:r1)%%r1+1}))))
    
    if (!is.null(TF)) colnames(TF)=fnames
  }
  attempts=0
  CRB= (max(columns)==1 & ncol(treatments)==1 & is.factor(treatments[,1]) & length(rows)==1 & hcf%%rows[1]==0)
  while (is.null(TF) & attempts<10) {
    attempts=attempts+1
    TF=treatments
    colnames(TF)=fnames
    if (!CRB) {
      for ( i in 1:strata) {
        if (hcf%%cumrows[i]!=0) TF=blocksOpt(TF,fblkDesign[,i],frowDesign[,i],nrowDesign[,i],fblkDesign[,i])
        if (columns[i]>1)       TF=blocksOpt(TF,fblkDesign[,i],fcolDesign[,i],ncolDesign[,i],frowDesign[,i])
      }
    }
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  rownames(TF)=NULL 

  # add back single rep treatments for nested stratum blocks only
  if (unstructured & min(replicates)==1 & max(replicates)>1) {
    nunits=sum(replicates*fulltrts)
    fullblocksizes=nunits
    for (i in 1:strata)
      fullblocksizes=Sizes(fullblocksizes,i)
    replications=rep(replicates,fulltrts)
    TrtsInBlocks= split( levels(TF[,1]) [TF[,1]] , rep(1:length(blocksizes),blocksizes))
    singleTF=split((1:sum(fulltrts))[replications==1], rep(1:length(fullblocksizes),(fullblocksizes-blocksizes)))

    for (i in names(singleTF))
      TrtsInBlocks[[i]]=sample(c(TrtsInBlocks[[i]] ,singleTF[[i]]))
    TF=data.frame(unlist(TrtsInBlocks))
    blocksizes=fullblocksizes
    Plots=factor(1:nunits)
    frowDesign=frowdesign[ rep( 1:length(blocksizes),  blocksizes ),,drop=FALSE]
  }

  if (max(columns)==1) {
    frowDesign=cbind(frowDesign,Plots)
    frowDesign=data.frame(lapply(1:ncol(frowDesign), function(r){sample(nlevels(frowDesign[,r]))[frowDesign[,r]]})) # Randomize labels - NB gives numeric columns
    frowDesign=cbind(frowDesign,TF)
    frowDesign=frowDesign[do.call(order, frowDesign), ] # re-order
    blocksizes=table(frowDesign[,ncol(frowDesign)-ncol(TF)-1])[unique(frowDesign[,ncol(frowDesign)-ncol(TF)-1])]
    TF=frowDesign[,c((ncol(frowDesign)-ncol(TF)+1):ncol(frowDesign)),drop=FALSE]
    for (i in 1 : ncol(treatments))
      if (is.factor(treatments[,i])) TF[,i]=as.factor(TF[,i])
    Design  = data.frame(frowdesign[rep(1:length(blocksizes),blocksizes ),],Plots, TF)
    rownames(Design)=NULL
    colnames(Design)=c(colnames(frowdesign),"Plots",fnames)
    V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))],Design[,(ncol(Design)-ncol(TF)-1)])
    V=lapply(V, function(x){ length(x) =max(blocksizes); x })
    Blocks.Plots=rep("",length(V))
    if (ncol(treatments)==1) {
      Plan=data.frame(frowdesign,Blocks.Plots,matrix( unlist(V), nrow=length(V),byrow=TRUE))
      colnames(Plan)= c(colnames(Plan[,1:(strata+1)]), c(1:(ncol(Plan)-strata-1)))
      rownames(Plan)=NULL
    }
    else Plan=NULL
  }

  if (max(columns)>1) {
    rdf = data.frame(frowdesign,fcoldesign,fblkdesign[,-1,drop=FALSE])[c(t(matrix(1:(3*strata),nrow=strata)))] # reorder columns
    rcDesign=data.frame(rdf[rep(1:length(blocksizes), blocksizes),],Plots)
    rcDesign=data.frame(lapply(1:ncol(rcDesign), function(r){ sample(nlevels(rcDesign[,r]))[rcDesign[,r]]})  ) # Randomize
    rcDesign=data.frame(rcDesign, TF)
    rcDesign=rcDesign[do.call(order,rcDesign), ] # re-order
    colnames(rcDesign)=NULL
    TF=rcDesign[,c((ncol(rcDesign)-ncol(TF)+1):ncol(rcDesign)),drop=FALSE]
    blocksizes=table(rcDesign[,ncol(rcDesign)-ncol(TF)-1])[unique(rcDesign[,ncol(rcDesign)-ncol(TF)-1])]
    Design  = data.frame( rdf[rep(1:length(blocksizes), blocksizes),],Plots,TF)  # rebuild factor levels
    colnames(Design)=c(colnames(rdf),"Plots",fnames)
    rownames(Design)=NULL
    if (max(blocksizes)>1 & ncol(treatments)==1)   {
      V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))],Design[,(ncol(Design)-ncol(TF)-1)]) # split on blocks
      V=lapply(V, function(x){ length(x) =max(blocksizes); x })
      Blocks.Plots=rep("",length(V))
      Plan=data.frame(rdf,Blocks.Plots,matrix( unlist(V), nrow=length(V),byrow=TRUE))
      colnames(Plan)=c(colnames(rdf),"Plots_In_Blocks",1:max(blocksizes))
      Plan=Plan[,-(ncol(Plan)-max(blocksizes)-1)]
      sets=nlevels(Plan[,ncol(rdf)-3])
      sizes=nrow(Plan)/sets
    } else if (max(blocksizes)==1 & ncol(treatments)==1) {
      V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))], Design[,(ncol(Design)-ncol(TF)-3)] ) # split on rows
      plan = rdf[seq(1,nrow(rdf)-columns[strata]+1,columns[strata]),-c(ncol(rdf)-1,ncol(rdf)) ,drop=FALSE]
      Columns=rep("",nrow(plan))
      Plan=data.frame(plan,Columns,matrix( unlist(V), nrow=length(V),byrow=TRUE))
      colnames(Plan)=c(colnames(plan),paste0("Stratum_",strata,".Cols"),1:columns[strata])
      rownames(Plan)=NULL
      sets=nlevels(Plan[,ncol(plan)-1])
      sizes=nrow(Plan)/sets
    } else Plan=NULL
    if (!is.null(Plan) & (strata>1)) {
      reorder=order(c(1:nrow(Plan)+rep(0:(sets-1),each=sizes),(1:sets)*(sizes+1)))
      Plan=rbind(Plan, setNames(data.frame( matrix(" ", nrow=sets,ncol=ncol(Plan)) ), names(Plan)))[reorder,]
    }
  }
  # efficiencies
  if (Latin) 
    Efficiencies=LatinEfficiencies(Design)
  else if (ncol(treatments)==1 & is.factor(treatments[,1]) & max(columns)>1) 
    Efficiencies=RowColEfficiencies(Design) 
  else if (ncol(treatments)==1 & is.factor(treatments[,1]) & max(replicates)==1)
    Efficiencies=data.frame(Strata = "Stratum_1", Blocks = 1, D_Efficiencies = 1, A_Efficiencies = 1, A_Bounds = 1)
  else if (ncol(treatments)==1 & is.factor(treatments[,1]))  
    Efficiencies=BlockEfficiencies(Design) 
  else if ( max(columns)>1)
    Efficiencies=FactRowColEffics(Design) 
  else if ( max(columns)==1)
    Efficiencies=FactBlocksEffics(Design)
  row.names(Efficiencies)=NULL

  # treatment replications
  if(ncol(treatments)==1 & is.factor(treatments[,1])) {
  TreatmentsTable=data.frame(table(Design[,ncol(Design)]))
  TreatmentsTable=TreatmentsTable[order( as.numeric(levels(TreatmentsTable[,1])) [TreatmentsTable[,1]]),]
  TreatmentsTable[]=lapply(TreatmentsTable, as.factor)
  colnames(TreatmentsTable)=c("Treatments","Replicates")
  rownames(TreatmentsTable)=1:nrow(TreatmentsTable)
  } else if (replicates%%1==0) {
    TreatmentsTable=data.frame(treatments,rep(replicates,nrow(treatments)))
    colnames(TreatmentsTable)=c(colnames(treatments),"Replicates")
    TreatmentsTable=TreatmentsTable[ do.call(order, TreatmentsTable), ]
  } else {
    TreatmentsTable=data.frame(treatments)
    colnames(TreatmentsTable)=colnames(treatments)
    TreatmentsTable=TreatmentsTable[ do.call(order, TreatmentsTable), ]
  }
  rownames(TreatmentsTable)=NULL
  } else {
      TF=data.frame(treatments[sample(1:nrow(treatments)), ,drop=FALSE])
      colnames(TF)=fnames
      rownames(TF)=NULL
      if(ncol(treatments)==1 & is.factor(treatments[,1]) ) {
        TreatmentsTable=data.frame(table(TF))
        TreatmentsTable[]=lapply(TreatmentsTable, as.factor)
        colnames(TreatmentsTable)=c("Treatments","Replicates")
      } else if (all(replicates%%1==0)) {
        TreatmentsTable=data.frame(treatments,rep(replicates,nrow(treatments)))
        colnames(TreatmentsTable)=c(colnames(treatments),"Replicates")
        TreatmentsTable=TreatmentsTable[ do.call(order, TreatmentsTable), ]
      } else {
        TreatmentsTable=data.frame(treatments)
        colnames(TreatmentsTable)=colnames(treatments)
        TreatmentsTable=TreatmentsTable[ do.call(order, TreatmentsTable), ]
      }
      rownames(TreatmentsTable)=1:nrow(TreatmentsTable)
      Efficiencies=data.frame(Strata = "Stratum_1", Blocks = 1, D_Efficiencies = 1, A_Efficiencies = 1, A_Bounds = 1)
      Plan=NULL
      Design=TF
    }
  list(Treatments=TreatmentsTable,model=model,DesignEfficiency=fractionalEff,BlocksEfficiency=Efficiencies,Plan=Plan,Design=Design,seed=seed,searches=searches,jumps=jumps)
}
