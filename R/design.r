#' @title General block and treatment designs.
#'
#' @description
#' Constructs block and treatment designs for any feasible combination of nested or crossed block 
#' factors and any feasible combination of  qualitative or quantitative level treatment factors.   
#' 
#' @param treatments a single treatment factor or a data frame containing one or more qualitative or 
#' quantitative level treatment factors.  
#' 
#' @param blocks a single blocks factor or a data frame containing one or more qualitative level block
#' factors in the required order of fitting.
#' 
#' @param treatments_model a model formula for the treatments design. 
#' 
#' @param weighting a weighting factor between 0 and 1 for the relative importance of the
#' block factor interaction effects versus the block factor main effects for a crossed blocks factor design. 
#'
#' @param seed an integer initializing the random number generator. 
#'
#' @param searches the maximum number of local optima searched at each stage of a
#'  design optimization
#'
#' @param jumps the number of pairwise random treatment swaps used to escape a local maxima. 
#'
#' @details
#' 
#' \code{treatments} is a factor or a data frame for one or more qualitative or quantitative
#' level treatment factors. The \code{treatments} object is the candidate set of treatments or treatment
#' factor combinations from which the treatment design is selected and should include
#' all the possible feasible treatment combinations that could or should be included in the final design. 
#' 
#' \code{blocks} is a factor or data frame for one or more qualitative block factors where the length or number of rows
#' of \code{blocks} is the number of plots in the design. If the design is a completely randomized design, 
#' the \code{blocks} object should be a blocks factor of the required length but with only a single factor level.
#' 
#' The ratio of the length of the \code{blocks} object to the length of the \code{treatments} object is the replication 
#' number and he integer part of the replication number, if any, defines the number of complete replications of the
#' \code{treatments} object while the fractional part, if any, defines a sample fraction of that size drawn from the
#' candidate set of \code{treatments} combinations.
#' 
#' The blocks design can comprise any ordered set of crossed or nested block factors provided only that the
#' nested blocks are arranged in a hierarchical order of nesting from largest to smallest. 
#' The block factors are added in sequence and each successively added blocks factor is optimised by maximising 
#' the current blocks model but with all the previously added blocks factors held constant.
#' 
#' N.B. crossed block designs can fail due to inherent aliasing. For example, \code{blocksdesign}
#' will always try to fit orthogonal row and column block effects which means that a row-and-column
#' design with two rows, two columns and two treatment replicates will automatically confound one treatment contrast
#' with the rows-by-columns interaction effect, which means that a full second-order block design will always be singular. 
#' 
#' \code{treatments_model} is a design formula for the \code{treatments} model based on the
#' \code{models} formula of the \code{\link[stats]{model.matrix}} package. The model fits factorial treatment contrasts
#' for qualitative level factors and polynomial treatment contrastst for quantitative level factors. The  default treatment
#' model assumes fully crossed factor effects for qualitative level factors and first-order effects for for quantitative 
#' level factors.
#' 
#' \code{weighting} is a weighting parameter between 0 and 1 which differentially weights the higher-order effects of
#' a crossed blocks design relative to the the first-order effects.
#' Designs with crossed blocks are usually assumed to fit an additive main block effects model but this is a very
#' strong assumption and wherever possible, crossed blocks designs should allow for higher-order interaction
#' effects. Setting the weighting parameter between 0 and 1 gives increasing importance to the higher-order 
#' effects in a blocks model relative to the first-order effects as the weighting increases from 0 to 1. 
#' The default weighting is 0.5 (only applies when a design has crossed blocks).
#' 
#' \code{seed} is an integer which can be used to set the initial seed for the \code{blocksdesign} random number generator. 
#' Normally it is best to use the \code{NULL} seed setting which allows \code{blocksdesign} to find its own random seed. 
#' 
#' The \code{blocks_model} output shows the overall achieved D and A-efficiency factors for each sequentially fitted blocks 
#' model. Efficiency factors are shown for a first-order and a second-order model for each sequentially added
#' block factor. For a fully nested blocks design, the two models are eqivalent and the two sets of efficiency factors will be
#' equal but for a general crossed blocks model the two sets of efficiency factors will be different and 
#' may provide some guidance on the best choice of weighting parameter for an efficient design. 
#' 
#' The definition of efficiency used by the \code{design} algorithm is the ratio of the generalized variance of the 
#' full treatment design relative to the generalized variance of the optimized block and treatment design. 
#' Using this definition, it is possible for quantitative level treatment designs to have efficiency factors greater
#' than one. Therefore the main use of efficiency factors is to compare the relative efficiencies of different
#'  optimizations of the same design.  
#' 
#' Outputs:
#'
#' The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing a randomized allocation of treatments to blocks. \cr
#'  \item  A table showing the fractional size of the treatment design and the 
#'  D-efficiency factors of that fraction. \cr
#'  \item  A table showing the blocks sub-model design and the D- and A-efficiency factor 
#'  of each successively fitted blocks sub-model. \cr
#' }
#' 
#' @return
#' \item{design}{The design layout showing the allocation of treatment and block design 
#' factors to individual plots.}
#' \item{treatments_model}{The fractional size of the treatment design together with the 
#' D-efficiency of that fraction.}
#' \item{blocks_model}{The blocks sub-model design together with the D and A-efficiency factor 
#' of each successively fitted blocks sub-model.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#' 
#' @references
#'
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' DURBAN, M., HACKETT, C., MCNICOL, J., NEWTON, A., THOMAS, W., & CURRIE, I. (2003). The practical use of semi-parametric models
#'  in field trials, Journal of Agric. Biological and Envir. Stats., 8, 48-66.
#'  
#' @examples
#' 
#' ## For optimum results, the number of searches may need to be increased in practice.
#' ## Designs can be rebuilt repeatedly to check that a near-optimum design has been found.  
#' 
#' ## 4 replicates of 12 treatments with 16 nested blocks of size 3
#' treatments = factor(1:12)
#' Blocks = factor(rep(1:4,each=12))
#' subBlocks = factor(rep(1:16,each=3))
#' blocks = data.frame(Blocks,subBlocks)
#' design(treatments,blocks)$blocks_model
#' 
#' ## 4 x 12 design for 4 replicates of 12 treatments with 16 nested blocks of size 3
#' ## only the default weighting (0.5) will ensure an optimal Trojan design 
#' treatments = factor(1:12)
#' MainCols = factor(rep(rep(1:4,each=3),4))
#' MainRows = factor(rep(1:4,each=12))
#' Cols = factor(rep(1:12,4))
#' blocks = data.frame(MainRows,MainCols,Cols)
#' \donttest{design(treatments,blocks,searches=200,weighting=0)$blocks_model
#' design(treatments,blocks,searches=200)$blocks_model
#' design(treatments,blocks,searches=200,weighting=1)$blocks_model}
#' 
#' ## 4 x 13 Row-and-column design for 4 replicates of 13 treatments 
#' ## Youden design Plan 13.5 Cochran and Cox (1957).
#' treatments=factor(1:13)
#' Rows =factor(rep(1:4,each=13))
#' Cols =factor(rep(1:13,4))
#' blocks =data.frame(Rows,Cols)
#' \donttest{design(treatments,blocks,searches=700)}
#' 
#' ## Durban - 272 treatments in a 16 x 34 design with nested rows-and-columns
#' data(durban) 
#' durban=durban[c(3,1,2,4,5)]
#' durban=durban[ do.call(order, durban), ]
#' treatments=data.frame(gen=durban$gen)
#' Reps = factor(rep(1:2,each=272))
#' Rows = factor(rep(1:16,each=34))
#' Col1 = factor(rep(rep(1:4,c(9,8,8,9)),16))
#' Col2 = factor(rep(rep(1:8,c(5,4,4,4,4,4,4,5)),16))
#' Col3 = factor(rep(1:34,16))
#' blocks = data.frame(Reps,Rows,Col1,Col2,Col3)
#' \donttest{design(treatments,blocks,searches=10)$blocks_model
#' ## Compare with efficiency factors of original design; Durban et al (2003)
#' blockEfficiencies(treatments,blocks)
#' } 
#' 
#' ## differential replication including single replicate treatments (13 to 24)
#' treatments=factor(c(rep(1:12,2), rep(13:24,1)))
#' Main=factor(rep(1:2,each=18))
#' Sub =factor(rep(1:6,each=6))
#' blocks =data.frame(Main,Sub)
#' design(treatments,blocks,searches=5)
#' 
#' ## 48 treatments in 2 replicate blocks of size 4 x 12 with 2 main rows and 3 main columns
#' treatments=factor(1:48)
#' replicates=factor(rep(1:2,each=48))
#' rows=factor(rep(rep(1:2,each=24),2))
#' cols=factor(rep(rep(1:3,each=8),4))
#' blocks=data.frame(replicates,cols,rows)
#' design(treatments,blocks,searches=5)
#' 
#' ## Factorial treatment designs defined by a treatments data frame and a factorial model equation.
#' ## For some examples, a repeat loop with a break based on an efficiency factor criterion is used
#' ## to search for a global optimal design from amongst a large number of local optimal designs 
#' 
#' ## Main effects of five 2-level factors in a half-fraction of a 2/2/2 nested blocks design
#' treatments = expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' blocks=data.frame(b1=factor(rep(1:2,each=8)),b2=factor(rep(1:4,each=4)),b3=factor(rep(1:8,each=2)))
#' treatments_model="F1 + F2 + F3 + F4 + F5"
#' \donttest{repeat {
#'  z=design(treatments,blocks,treatments_model,searches=5)
#' if ( z$blocks_model[3,3]==1  ) break }
#' z}
#' 
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' treatments=expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' blocks=factor(rep(1:4,each=8))
#' treatments_model="(F1+F2+F3+F4+F5)^2"
#' design(treatments,blocks,treatments_model,searches=5)
#' 
#' # Main effects of five 2-level factors in a half-fraction of a 4 x 4 row-and column design.
#' treatments = expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),
#' F5=factor(1:2))
#' blocks=data.frame( rows=factor(rep(1:4,each=4)), cols=factor(rep(1:4,4)))
#' treatments_model="~ F1+F2+F3+F4+F5"
#' design(treatments,blocks,treatments_model,searches=20)
#' 
#' # Quadratic regression for one 6-level numeric factor in 2 randomized
#' #  blocks assuming 10/6 fraction
#' treatments=expand.grid(X=1:6)
#' blocks=factor(rep(1:2,each=5))
#' treatments_model=" ~ poly(X,2)"
#' design(treatments,blocks,treatments_model,searches=5) 
#' 
#' # First-order model for 1/3rd fraction of four qualitative 3-level factors in 3  blocks
#' treatments=expand.grid(F1=factor(1:3),F2=factor(1:3),F3=factor(1:3),F4=factor(1:3))
#' blocks=factor(rep(1:3,each=9))
#' treatments_model=" ~ F1+F2+F3+F4"
#' \donttest{repeat {
#' z=design(treatments,blocks,treatments_model,searches=10)
#' if ( z$blocks_model[1,3]==1  ) break }
#' z}
#' 
#' # Second-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' treatments=expand.grid( F1=factor(1:3), F2=factor(1:3), F3=factor(1:3), F4=factor(1:3), 
#' F5=factor(1:3))
#' blocks=factor(rep(1:3,each=27))
#' treatments_model=" ~ (F1+F2+F3+F4+F5)^2"
#' \donttest{repeat {
#' z=design(treatments,blocks,treatments_model,searches=25)
#' if ( z$blocks_model[1,3]==1  ) break }
#' z}
#' 
#' # Second-order model for two qualitative and two quantitative level factors in 4 randomized blocks
#' treatments=expand.grid(F1=factor(1:2),F2=factor(1:3),V1=1:3,V2=1:4)
#' blocks=factor(rep(1:4,each=18))
#' treatments_model = " ~ F1 + F2 + poly(V1,2) +  poly(V2,2) + (poly(V1,1)+F1+F2):(poly(V2,1)+F1+F2) "
#'  \donttest{design(treatments,blocks,treatments_model,searches=5)}
#'  
#' # Plackett and Burman design for eleven 2-level factors in 12 runs (needs large number of searches)
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2),
#' F6=factor(1:2),F7=factor(1:2),F8=factor(1:2),F9=factor(1:2),F10=factor(1:2),F11=factor(1:2))
#' blocks=factor(rep(1,12))
#' model=model="~ F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11"
#' \donttest{design(GF,blocks,model,searches=25)}
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames 
#' @importFrom lme4 lmer
#' 
  design = function(treatments,blocks,treatments_model=NULL,weighting=0.5,searches=NULL,seed=NULL,jumps=1) {
  options(contrasts=c('contr.SAS','contr.poly'))
  tol=1e-6
  TF=treatments
  if (is.factor(TF)) TF=data.frame(TF)
  if (!is.data.frame(TF)) stop("treatments must be a data frame or a factor")
  BF=blocks
  if (is.factor(BF)) BF=data.frame(BF)
  if (!is.data.frame(BF)) stop("blocks must be a data frame or a factor")
  if (!is.numeric(weighting))  stop("weighting must be a number between 0 and 1")
  if (weighting<0 | weighting>1)  stop("weighting must be a number between 0 and 1")
  if (is.null(treatments_model)) treatments_model=paste("~",paste(colnames(TF),collapse="*"))
  if (!substring(trimws(treatments_model, "l"), 1,1)=="~") treatments_model=paste("~",treatments_model)
  if (!is.null(seed)) set.seed(seed)
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) 
    stop(" number of jumps parameter is invalid (max is 10) ")
  nunits=nrow(BF)
  if (is.null(searches)) 
    if (nunits<100) 
      searches=50 else if (nunits<1000) 
        searches=5000%/%nunits else if (nunits<5000) 
          searches=5000%/%nunits else 
            searches=1
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
 
   # *******************************************************************************************************************
   # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
   # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
   # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
   # *******************************************************************************************************************
   UpDate=function(V,t,b) {
     VTTt=crossprod(V[1:length(t),1:length(t),drop=FALSE],t)
     VBBb=crossprod(V[(1+length(t)):ncol(V),(1+length(t)):ncol(V),drop=FALSE]  ,b)
     VTBt=crossprod(V[1:length(t),(1+length(t)):ncol(V),drop=FALSE],t)
     VTBb=crossprod(V[(1+length(t)):ncol(V),1:length(t),drop=FALSE],b)
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
     V[1:length(t),1:length(t)] = V[1:length(t),1:length(t), drop=FALSE]-
       (tcrossprod(f1)*(1-b) - tcrossprod(g1)*(1+a) + (tcrossprod(g1,f1)+tcrossprod(f1,g1))*c)/((1+a)*(1-b)+c*c)
     V[(1+length(t)):ncol(V),(1+length(t)):ncol(V)] = V[(1+length(t)):ncol(V),(1+length(t)):ncol(V),drop=FALSE]- 
       (tcrossprod(f2)*(1-b) - tcrossprod(g2)*(1+a) + (tcrossprod(g2,f2)+tcrossprod(f2,g2))*c)/((1+a)*(1-b)+c*c)
     V[1:length(t),(1+length(t)):ncol(V)] = V[1:length(t),(1+length(t)):ncol(V),drop=FALSE]-
       (tcrossprod(f1,f2)*(1-b) - tcrossprod(g1,g2)*(1+a) + (tcrossprod(g1,f2)+tcrossprod(f1,g2))*c)/((1+a)*(1-b)+c*c)
     V[(1+length(t)):ncol(V), 1:length(t)]= t(V[1:length(t),(1+length(t)):ncol(V),drop=FALSE])
     return(V)
   }
  # ********************************************************************************************************************
  # Random swaps
  # ********************************************************************************************************************    
  Swaps=function(TF,pivot,rank,mainblocks) {
    candidates=NULL
    while (is.null(candidates)) {
        s1=sample(pivot[(1+rank):nunits],1)
        swaptrts=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
        candidates = (1:nunits)[swaptrts & mainblocks==mainblocks[s1]]
    }
    if ( length(candidates)>1 ) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  }
  # **********************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient, 
  # random swaps with positive selection are used to to increase design rank
  # **********************************************************************************************************************
  NonSingular=function(TF,TM,BM,mainblocks,searches) {
    if(searches==0)  return(list(TF=TF,TM=TM))
    if (ncol(TM)+ncol(BM)>(nunits-1)) stop(" Block design has too many parameters")
    fullrank=ncol(TM)+ncol(BM) # fullrank not more than (nunits-1)
    Q=qr(t(cbind(TM,BM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:500) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      s=Swaps(TF,pivot,rank,mainblocks)
      TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]
      Q=qr(t(cbind(TM,BM)))
      if (Q$rank>rank) {
        TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]
        rank=Q$rank
        pivot=Q$pivot
      } else TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]
    }
    stop(" Unable to find an initial non-singular design - try different initial seed or simplify the block design")
  }
  # *********************************************************************************************************************
  #  dMat is a matrix where the ith, jth matrix element is the change in the D-max criterion of the block design 
  # due to swapping ith and jth treatments of a sample s of the set of plots 
  # *********************************************************************************************************************
  dMat=function(TM,BM,V,s,dim0,dim1,dim2) {
    sTM=TM[s,,drop=FALSE]
    sBM=BM[s,,drop=FALSE]
    TMT=crossprod(t(crossprod(t(sTM),V[1:dim0,1:dim0,drop=FALSE])),t(sTM))
    TMB=crossprod(t(crossprod(t(sTM),V[1:dim0,(dim0+1):ncol(V), drop=FALSE])),t(sBM))
    BMB=crossprod(t(crossprod(t(sBM),V[(dim0+1):ncol(V),(dim0+1):ncol(V),drop=FALSE])),t(sBM))
    TMB=sweep(TMB,1,diag(TMB))
    TMT=sweep(TMT,1,diag(TMT))
    BMB=sweep(BMB,1,diag(BMB))
    dMat=(1+TMB+t(TMB))**2 - (TMT + t(TMT))*(BMB + t(BMB))
    return(dMat)
  }
  # ******************************************************************************************************************
  # Maximises the matrix dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially but later a full search is used to ensure steepest ascent optimization.
  # *******************************************************************************************************************
  DMax=function(TF,TM,BM,V,mainBlocks,dim0,dim1,dim2) {
    locrelD=1
    mainSets=tabulate(mainBlocks)
    nSamp=pmin(rep(8,nlevels(mainBlocks)), mainSets)
    repeat {
      kmax=1
      for (k in 1: nlevels(mainBlocks)) {
        s=sort(sample((1:nunits)[mainBlocks==levels(mainBlocks)[k]] , nSamp[k])) 
        dMat=dMat(TM,BM,V,s,dim0,dim1,dim2)
        z=which(dMat == max(dMat,na.rm=TRUE), arr.ind = TRUE)[1,]
        if (dMat[z[1],z[2]]>kmax) {
          kmax=dMat[z[1],z[2]]
          pi=s[z[1]]
          pj=s[z[2]]
        } 
      }
      if (kmax>(1+tol)) {
        locrelD=locrelD*kmax
        V=UpDate(V,(TM[pi,]-TM[pj,]),(BM[pj,]-BM[pi,])  )
        TM[c(pi,pj),]=TM[c(pj,pi),]
        TF[c(pi,pj),]=TF[c(pj,pi),]
      }  else if (sum(nSamp) == nunits) break
      else nSamp=pmin(mainSets,2*nSamp)
    }
    list(V=V,locrelD=locrelD,TF=TF,TM=TM)
  }
  # ************************************************************************************************************************
  # Optimize the nested blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ************************************************************************************************************************
   blocksOpt=function(TF,TM,BF,searches,nunits) {
     dim0=ncol(TM)
    mainBlocks=data.frame(factor(rep(1,nrow(BF))), matrix(nrow=nrow(BF),ncol=ncol(BF)))
    for (i in 1: ncol(BF)) mainBlocks[,i+1]= droplevels(interaction( mainBlocks[,i],BF[,i] ))
    colnames(mainBlocks)=c("mean", unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse=".")})) )
    for (u in 1:ncol(BF)) {
      BM1=scale(model.matrix(as.formula(paste("~", paste0(colnames(BF)[1:u],collapse="+"))),BF), center = TRUE, scale = FALSE)
      Q=qr(BM1)
      BM1 = BM1[,Q$pivot[1:Q$rank],drop=FALSE] # full rank
      BM1 = qr.Q(qr(BM1)) # orthogonal basis 
      dim1=ncol(BM1)
      if ((1+dim0+dim1)>nunits) stop("Additive block design has too many fixed effects for number of plots ")
      if(u>1 & weighting>0){
        BM2=scale(model.matrix(as.formula(paste("~",paste0("(",paste0(colnames(BF)[1:u],collapse="+"),")^2"))),BF), center = TRUE, scale = FALSE)
        Q=qr(BM2)
        BM2 = BM2[,Q$pivot[1:Q$rank],drop=FALSE] # full rank
        BM2 = qr.Q(qr(BM2)) # orthogonal basis 
        } else BM2=BM1
      dim2=ncol(BM2)
      deg2=(weighting>0 & dim2>dim1 & (dim0+dim2)<nunits)
      if (deg2) BM=BM2 else BM=BM1
      NS=NonSingular(TF,TM,BM,mainBlocks[,u],searches)
      TM=NS$TM
      TF=NS$TF  
      globrelD=0
      relD=1
      globTM=TM
      globTF=TF
      if (deg2) BM[,(dim1+1):dim2] = BM[,(dim1+1):dim2]*weighting
      Info=crossprod(cbind(TM,BM))
      if (deg2) 
        Info[(1+dim0+dim1):(dim0+dim2),(1+dim0+dim1):(dim0+dim2)] = Info[(1+dim0+dim1):(dim0+dim2),(1+dim0+dim1):(dim0+dim2)]/(weighting^2)
      V=chol2inv(chol(Info))
      for (r in 1:searches) {
        dmax=DMax(TF,TM,BM,V,mainBlocks[,u],dim0,dim1,dim2)
        if (dmax$locrelD>(1+tol)) {
          relD=relD*dmax$locrelD
          TM=dmax$TM
          TF=dmax$TF
          V=dmax$V
          if (relD>globrelD) {
            globTM=TM 
            globTF=TF 
            globrelD=relD
          }
        }
        if (r==searches) break
        for (iswap in 1:jumps) {
          counter=0
          repeat {
            counter=counter+1
            s1=sample(seq_len(nunits),1)
            available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
            z= seq_len(nunits)[mainBlocks[,u]==mainBlocks[s1,u] & mainBlocks[,u+1]!=mainBlocks[s1,u+1] & available]  
            if (length(z)==0 & counter<501) next 
            else if (length(z)==0 & counter>500) break
            if (length(z)>1) s=c(s1,sample(z,1)) else s=c(s1,z)
            testDswap=dMat(TM,BM,V,s,dim0,dim1,dim2)[2,1]
            if (testDswap<tol & counter<501) next
            else if (testDswap<tol & counter>500) break
            Dswap=testDswap
            break
          }
          if (counter>500) break
          relD=relD*Dswap 
          V=UpDate(V,(TM[s[1],]-TM[s[2],]),(BM[s[2],]-BM[s[1],]))
          TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
          TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]  
        } #jumps
      } #searches
      TM=globTM
      TF=globTF 
    } # list length
    Effics=blockEfficiencies(TF,BF,treatments_model)
    list(TF=TF,TM=TM,Effics=Effics)
  }
  # ******************************************************************************************************************
  # Updates variance matrix for swapped rows where mi is swapped out and mj is swapped in
  # ******************************************************************************************************************
  fractUpDate=function(V,mi,mj) {
    f=crossprod(V,mj)
    g=crossprod(V,mi)
    a=as.numeric(crossprod(mj,f))
    b=as.numeric(crossprod(mi,g))
    c=as.numeric(crossprod(mi,f))
    W=g*sqrt(1+a)-f*c/sqrt(1+a)
    U=f*sqrt(1-b+(c*c)/(1+a))
    V=V-(tcrossprod(U)-tcrossprod(W))/((1+a)*(1-b)+c*c)
    return(V)
  }
  

  # *********************************************************************************************************************
  # Fractional factorials
  # *********************************************************************************************************************
  factorial=function(TF,treatments_model,searches,nunits) {
    #TF is the supplied set of treatment factors
    allfact= (all(sapply(TF,is.factor))) 
    replicates=nunits/nrow(TF)
    # fully replicated base which remains unchanged
    nbase=nrow(TF)*floor(replicates) 
    TM=model.matrix(as.formula(treatments_model),TF)
    TM = qr.Q(qr(TM)) # orthogonal basis
    if (ncol(TM)>nunits) stop("too many treatment parameters to be estimated for the given number of experimental units")

    if (replicates==round(replicates)) {
      rerand=unlist(lapply(1:replicates,function(i){ sample(seq_len(nrow(TF)))}))
      return(list(TF=TF[rerand,,drop=FALSE],TM=TM[rerand,,drop=FALSE],eff=1,fraction=replicates))
    }
 
    rerand=vector("list", searches)
    for (t in 1:searches) {
      for (z in seq_len(1000)) {
        rerand[[t]]=unlist(lapply(1:ceiling(replicates),function(j){ sample(seq_len(nrow(TF)))}))
        if (qr(TM[ rerand[[t]][1:nunits],,drop=FALSE ])$rank==ncol(TM))break
      }
      if (z>999) stop("Unable to find a non-singular treatment start design of the required size by random search ")
    }
    
    gDfrac=0
    for (i in 1:searches) {
      fTM=TM[rerand[[i]],,drop=FALSE]
      fTF=TF[rerand[[i]],,drop=FALSE]
      V=chol2inv(chol(crossprod(fTM[1:nunits,,drop=FALSE]))) # V is the variance matrix of the the non-singular starting design
      repeat {
        # plot numbers of a random sample of the unused treatment combinations
        if ((nrow(fTM)-nunits)>1)
        s= sort(sample( (nunits+1):nrow(fTM),min( (nrow(fTM)-nunits),2000)))
        else s=nunits+1
        bfTM=fTM[(nbase+1):nunits,,drop=FALSE] # included plots that can be swapped out
        sfTM=fTM[s,,drop=FALSE] # excluded plots that can be swapped in
        M1VM2=       tcrossprod(tcrossprod( bfTM,V),sfTM)
        M1VM1=1-diag(tcrossprod(tcrossprod( bfTM,V),bfTM))
        M2VM2=1+diag(tcrossprod(tcrossprod( sfTM,V),sfTM))
        Z=M1VM2**2 + tcrossprod(M1VM1,M2VM2)
        z=which(Z == max(Z), arr.ind = TRUE)[1,]
        if (  Z[z[1],z[2]]<(1+tol)) break
        p1= nbase+z[1]
        p2= s[z[2]]
        V=fractUpDate(V, fTM[p1,] , fTM[p2,])   # parameters(V,row_swappedout,row_swappedin) 
        fTM[c(p1,p2),] = fTM[c(p2,p1),]
        fTF[c(p1,p2),] = fTF[c(p2,p1),]
      }
      Dfrac=exp(determinant(crossprod(fTM[1:nunits,,drop=FALSE]),logarithm = TRUE)$modulus[1]/ncol(fTM))/replicates

      if (Dfrac>gDfrac) {
        gTM=fTM[1:nunits,,drop=FALSE]
        gTF=fTF[1:nunits,,drop=FALSE]
        gDfrac=Dfrac
      }
      if (allfact & isTRUE(all.equal(gDfrac,1))) break
    }
    return(list(TF=gTF,TM=gTM,eff=gDfrac,fraction=replicates))
  }
  

  # *************************************************************************************************************************
  # Main design program tests, optimizes treatment design then optimizes block design sequentially for each added block factor
  # *************************************************************************************************************************
  Z=factorial(TF,treatments_model,searches,nunits)
  TM=scale(Z$TM, center = TRUE, scale = FALSE)[,-1,drop=FALSE]
  TF=Z$TF
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  treatsModel=data.frame("Treatment model"=treatments_model,"Treatment fraction"=Z$fraction,"D-Efficiency"=Z$eff)
  if (max(sapply(BF, nlevels))>1) {
  Opt=blocksOpt(TF,TM,BF,searches,nunits)
  TM=Opt$TM
  TF=Opt$TF
  blocksModel=Opt$Effics
  } else blocksModel=NULL
  Design=cbind(BF, TF)
  rownames(Design)=NULL
  list(design=Design,treatments_model=treatsModel,blocks_model=blocksModel,weighting=weighting,seed=seed,searches=searches,jumps=jumps)
  }
  