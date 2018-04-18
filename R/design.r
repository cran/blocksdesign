#' @title General block and treatment designs.
#'
#' @description
#' Constructs D-optimal block and treatment designs for feasible combinations of nested or crossed block 
#' factors and feasible linear treatment models.   
#' 
#' @param treatments a single treatment factor or a data frame containing multiple treatment factors. 
#' Treatment factors can be qualitative or quantitative (numeric) level factors.  
#' 
#' @param blocks a single blocks factor or a data frame containing multiple block factors in the
#'  required order of fitting. Block factors must be qualitative level factors.  
#' 
#' @param treatments_model a model formula for the required treatments design. The default 
#' is a fully crossed factorial model.
#' 
#' @param weighting a weighting factor between 0 and 1 for general crossed block designs. 
#' The default weighting is 0.5. 
#'
#' @param seed an integer initializing the random number generator. The default is a positive 
#' integer seed chosen randomly.
#'
#' @param searches the maximum number of local optima searched at each stage of a treatment 
#' and block design optimization. The default depends on the design size.
#'
#' @param jumps the number of pairwise random treatment swaps used to escape a local maxima. 
#' The default is a single swap.
#'
#' @details
#'
#' \code{treatments} is a treatments factor or a data frame for two or more columns of 
#' treatment factors where the factors can be qualitative or quantitative. 
#' 
#' \code{blocks} is a blocks factor or a data frame for two or more columns of qualitative block factors
#' where the default is a single level factor of length equal to the length of the \code{treatments} factors. 
#' 
#' \code{treatments_model} is a design formula for the \code{treatments} factors based on the
#' \code{models} formula of the \code{\link[stats]{model.matrix}} package. The default assumes
#'  a complete factorial design.
#' 
#' The total number of design plots is defined by the length of the \code{blocks} factors, if present, 
#' otherwise by the length of the \code{treatments} factors. 
#' 
#' The treatments are replicated in the ratio of the total number of design plots to the length
#' of the \code{treatments} factors. The integer part of the ratio, possibly zero, defines the number of
#' complete replications of the \code{treatments} factors while the fractional part of the ratio, if any, 
#' defines a sample fraction of that size drawn from the set of \code{treatments} factor combinations. 
#' Samples are drawn without replacement and are chosen to optimize the D-optimality of the \code{treatments_model}. 
#' 
#' The \code{design} algorithm builds the blocks design by sequentially adding \code{blocks} factors
#' in the column order of the \code{blocks} data frame. Each block factor is optimized  
#' conditional on all preceding block factors but ignoring all succeeding block factors.
#' This method of optimization allows the blocking factors to be fitted in order of importance with the 
#' largest and most important blocks fitted first and the smaller and less important blocks fitted subsequently. 
#' 
#' For fully nested block designs, the block effects are optimized recursively assuming a simple additive block effects model.
#' For general crossed block designs, the crossed blocks intersections may contain blocks of two or more plots and 
#' then both the main effects and the intersection blocks must be optimized simultaneously.
#' For this situation, the \code{design} algorithm 
#' provides an option for optimizing a weighted combination of the
#' additive crossed blocks information matrix and the multiplicative crossed blocks information matrix. If the \code{weighting} 
#' factor is zero, the design is a fully additive crossed blocks model, if the \code{weighting} factor is one the design is a 
#' fully multiplicative crossed blocks model while for any intermediate \code{weighting}, the design is a
#' compromise between a fully additive and a fully multiplicative crossed blocks model.
#'  
#' The \code{blocks_model} output shows the overall achieved D-efficiency for each sequentially fitted blocks models. 
#' For nested blocks, the \code{blocks_model} output shows the efficiency factors for each successively nested
#' blocks design whereas for crossed blocks the \code{blocks_model} output shows the efficiency 
#' factors for both the additive and, where available, for the multiplicative effects of each sequentially fitted crossed 
#' blocks design. Comparison of the efficiency factors of the weighted crossed block designs using different weightings 
#' will provide guidance on the best choice of weighting for an efficient design. 
#' 
#' The efficiency factor used here is the ratio of the generalized variance of the full treatment design
#' relative to the generalized variance of the optimized block and 
#' treatment design. Using this definition, it is possible for quantitative level treatment designs to have efficiency 
#' factors greater than one therefore the main use of efficiency factors is to compare different optimizations
#' of the same design.  
#' 
#'   
#' Outputs:
#'
#' The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing a randomized allocation of treatments to blocks. \cr
#'  \item  A table showing the fractional size of the treatment design and the 
#'  D-efficiency factor of that fraction. \cr
#'  \item  A table showing the blocks sub-model design and the D-efficiency factor 
#'  of each successively fitted blocks sub-model. \cr
#' }
#' 
#' @return
#' \item{design}{The design layout showing the allocation of treatment and block design 
#' factors to individual plots.}
#' \item{treatments_model}{The fractional size of the treatment design together with the 
#' D-efficiency of that fraction.}
#' \item{blocks_model}{The blocks sub-model design together with the D-efficiency factor 
#' of each successively fitted blocks sub-model.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#' 
#' @references
#'
#' Cochran W.G & Cox G.M (1957) Experimental Designs 2nd Edition John Wiley & sons
#'  
#' Edmondson R. N. (1998). Trojan square and incomplete Trojan square designs for crop research. 
#' Journal of Agricultural Science, Cambridge, 131, pp.135-142
#'
#' @examples
#' 
#' ## For optimum results, the number of searches may need to be increased in practice.
#' ## Designs should be rebuilt repeatedly to check that a near-optimum design has been found.  
#' 
#' ## 4 replicates of 12 treatments with 16 nested blocks of size 3
#' treatments = factor(1:12)
#' Blocks = factor(rep(1:4,each=12))
#' subBlocks = factor(rep(1:16,each=3))
#' blocks = data.frame(Blocks,subBlocks)
#' design(treatments,blocks)$blocks_model
#' 
#' ## 4 x 12 design for 4 replicates of 12 treatments with 16 nested blocks of size 3
#' ## only the intermediate weighting will give an optimal Trojan design 
#' treatments = factor(1:12)
#' MainCols = factor(rep(rep(1:4,each=3),4))
#' MainRows = factor(rep(1:4,each=12))
#' Columns = factor(rep(1:12,4))
#' blocks = data.frame(MainCols,MainRows,Columns)
#' \donttest{design(treatments,blocks,searches=100,weighting=0)$blocks_model
#' design(treatments,blocks,searches=100,weighting=0.5)$blocks_model
#' design(treatments,blocks,searches=100,weighting=1)$blocks_model}
#'  
#' ## 4 x 13 Row-and-column design for 4 replicates of 13 treatments 
#' ## Youden design Plan 13.5 Cochran and Cox (1957).
#' treatments=factor(1:13)
#' Rows =factor(rep(1:4,each=13))
#' Cols =factor(rep(1:13,4))
#' blocks =data.frame(Rows,Cols)
#' \donttest{design(treatments,blocks,searches=500)}
#' 
#' ## Two replicates of 272 treatments in a 16 x 34 design with nested row blocks
#' treatments=factor(1:272)
#' Reps = factor(rep(1:2,each=272))
#' Rows = factor(rep(1:16,each=34))
#' MainCols = factor(rep(rep(1:4,c(9,8,8,9)),16))
#' SubCols = factor(rep(1:34,16))
#' blocks = data.frame(Reps,Rows,MainCols,SubCols)
#' \donttest{design(treatments,blocks,searches=1)$blocks_model}
#' 
#' ## differential replication including single replicate treatments
#' treatments=factor(c(rep(1:12,2), rep(13:24,1)))
#' Main=factor(rep(1:2,each=18))
#' Sub =factor(rep(1:6,each=6))
#' blocks =data.frame(Main,Sub)
#' design(treatments,blocks,searches=5)
#' 
#' 
#' ## Factorial treatment designs defined by a treatments data frame and a factorial model equation.
#' 
#' ## Main effects of five 2-level factors in a half-fraction of a 2/2/2 nested blocks design
#' treatments = expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' blocks=data.frame(b1=factor(rep(1:2,each=8)),b2=factor(rep(1:4,each=4)),b3=factor(rep(1:8,each=2)))
#' treatments_model="F1 + F2 + F3 + F4 + F5"
#' design(treatments,blocks,treatments_model,searches=5)
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
#' design(treatments,blocks,treatments_model,searches=5)
#' 
#' # Second-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' treatments=expand.grid( F1=factor(1:3), F2=factor(1:3), F3=factor(1:3), F4=factor(1:3), 
#' F5=factor(1:3))
#' blocks=factor(rep(1:3,each=27))
#' treatments_model=" ~ (F1+F2+F3+F4+F5)^2"
#'  \donttest{design(treatments,blocks,treatments_model,searches=100)}
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
#' \donttest{design(GF,blocks,model)}
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames
#'
 design = function(treatments,blocks=NULL,treatments_model=NULL,weighting=0.5,searches=NULL,seed=NULL,jumps=1) {
  options(contrasts=c('contr.SAS','contr.poly'))
  tol=1e-6
  TF=treatments
  BF=blocks
  if (is.factor(TF)) TF=data.frame(TF)
  if (is.null(BF)) BF=factor(rep(1,nrow(TF))) 
  if (is.factor(BF)) BF=data.frame(BF)
  if (!is.data.frame(TF)) stop("treatments must be a treatments data frame or a treatments factor")
  if (!is.data.frame(BF)) stop("blocks must be a blocks data frame or a blocks factor")
  if (!all(sapply(BF,is.factor))) stop("all blocks must be factors")
  if (!is.numeric(weighting))  stop("weighting must be a number between 0 and 1")
  if (weighting<0 | weighting>1)  stop("weighting must be a number between 0 and 1")
  if (is.null(treatments_model)) treatments_model=paste("~",paste(colnames(TF),collapse="*"))
  if (!substring(trimws(treatments_model, "l"), 1,1)=="~") treatments_model=paste("~",treatments_model)
  if (is.null(seed)) seed=sample(10000,1)
  if (is.na(seed) | !is.finite(seed) | is.nan(seed) | seed%%1!=0 | seed<0 ) stop(" seed parameter invalid  ")
  set.seed(seed)
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
   UpDate=function(VTT,VBB,VTB,t,b) {
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
   
   # *********************************************************************************************************************
   #  dMat is a matrix where the ith, jth matrix element is the change in the D-max criterion of the block design 
   # due to swapping ith and jth treatments of a sample s of the set of plots 
   # *********************************************************************************************************************
   
   dMat=function(TM,BM,VTT,VBB,VTB,s) {
     TMB=crossprod(t(crossprod(t(TM[s,,drop=FALSE]),VTB)),t(BM[s,,drop=FALSE]))
     TMT=crossprod(t(crossprod(t(TM[s,,drop=FALSE]),VTT)),t(TM[s,,drop=FALSE]))
     BMB=crossprod(t(crossprod(t(BM[s,,drop=FALSE]),VBB)),t(BM[s,,drop=FALSE]))
     TMB=sweep(TMB,1,diag(TMB))
     TMT=sweep(TMT,1,diag(TMT))
     BMB=sweep(BMB,1,diag(BMB))
     dMat=(1+TMB+t(TMB))**2 - (TMT + t(TMT))*(BMB + t(BMB))
   }
   
   # ******************************************************************************************************************
   # Maximises the matrix dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
   # Sampling is used initially but later a full search is used to ensure steepest ascent optimization.
   # *******************************************************************************************************************
   DMax=function(TF,TM,aVTT,aVBB,aVTB,aBM,mVTT,mVBB,mVTB,mBM,main) {
    locrelD=1
    mainSets=tabulate(main)
    nSamp=pmin(rep(8,nlevels(main)), mainSets)
    repeat {
      kmax=1
      for (k in 1: nlevels(main)) {
        s=sort(sample((1:nunits)[main==levels(main)[k]] , nSamp[k])) 
        adMat=dMat(TM,aBM,aVTT,aVBB,aVTB,s)
         if (!is.null(mBM)) { 
          mdMat=dMat(TM,mBM,mVTT,mVBB,mVTB,s)
          adMat[adMat<1]=-99
          mdMat[mdMat<1]=-99
          dMat=(1-weighting)*adMat+weighting*mdMat
        } else dMat=adMat
        z=which(dMat == max(dMat), arr.ind = TRUE)[1,]
        if (dMat[z[1],z[2]]>kmax) {
          kmax=dMat[z[1],z[2]]
          pi=s[z[1]]
          pj=s[z[2]]
        } 
      }
      if (kmax>(1+tol)) {
        locrelD=locrelD*kmax
        t=TM[pi,]-TM[pj,]
        ab=aBM[pj,]-aBM[pi,]
        aup=UpDate(aVTT,aVBB,aVTB,t,ab)
        aVTT=aup$VTT
        aVBB=aup$VBB
        aVTB=aup$VTB
        if (!is.null(mBM)) { 
          mb=mBM[pj,]-mBM[pi,]
          mup=UpDate(mVTT,mVBB,mVTB,t,mb)
          mVTT=mup$VTT
          mVBB=mup$VBB
          mVTB=mup$VTB
        }
        TM[c(pi,pj),]=TM[c(pj,pi),]
        TF[c(pi,pj),]=TF[c(pj,pi),]
      }  else if (sum(nSamp) == nunits) break
      else nSamp=pmin(mainSets,2*nSamp)
    }
    list(aVTT=aVTT,aVBB=aVBB,aVTB=aVTB,mVTT=mVTT,mVBB=mVBB,mVTB=mVTB,locrelD=locrelD,TF=TF,TM=TM)
   }
 
  # ********************************************************************************************************************
  # Random swaps
  # ********************************************************************************************************************    
  Swaps=function(TF,pivot,rank,main) {
    candidates=NULL
    while (is.null(candidates)) {
        s1=sample(pivot[(1+rank):nunits],1)
        swaptrts=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
        candidates = (1:nunits)[swaptrts & main==main[s1]]
    }
    if ( length(candidates)>1 ) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  }
  # **********************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient, 
  # random swaps with positive selection are used to to increase design rank
  # **********************************************************************************************************************
  NonSingular=function(TF,TM,BM,main) {
    if (ncol(TM)+ncol(BM)>(nunits-1)) return(list(TF=NULL,TM=NULL))
    fullrank=ncol(TM)+ncol(BM) # fullrank not more than (nunits-1)
    Q=qr(t(cbind(TM,BM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:500) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      s=Swaps(TF,pivot,rank,main)
      TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]
      Q=qr(t(cbind(TM,BM)))
      if (Q$rank>rank) {
        TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]
        rank=Q$rank
        pivot=Q$pivot
      } else TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]
    }
    return(list(TF=NULL,TM=NULL))
  }
  # ************************************************************************************************************************
  # Optimize the nested blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ************************************************************************************************************************
  blocksOpt=function(TF,TM,BF) {
    Effics=matrix(ncol = 3, nrow = 0)
    Dbase=determinant(crossprod(TM),logarithm = TRUE)$modulus/ncol(TM)
    for (i in 1:ncol(BF)) {
      if (i>1) mainblocks=droplevels(interaction(BF[,1:(i-1)])) else  mainblocks=factor(rep(1,nunits))
      subblocks=nlevels(droplevels(interaction(BF[,1:i])))
      modBF=BF[,1:i,drop=FALSE]
      addform=paste("~", paste(colnames(modBF),collapse="+"))
      aBM=scale(model.matrix(as.formula(addform),modBF), center = TRUE, scale = FALSE)
      Q=qr(aBM)
      aBM = aBM[,Q$pivot[1:Q$rank],drop=FALSE]
      aBM = qr.Q(qr(aBM)) # orthogonal basis
      if ( (ncol(TM)+ncol(aBM)+1) > nunits) stop("Additive block design has too many fixed effects for assumed number of plots ")
      NS=NonSingular(TF,TM,aBM,mainblocks)
      if (is.null(NS$TM) | is.null(NS$TF)) stop(" Unable to find an initial non-singular design for this choice of block design")
      TM=NS$TM
      TF=NS$TF   
      mVTT=NULL
      mVBB=NULL
      mVTB=NULL
        if ( (ncol(TM) + subblocks ) < nunits & (subblocks-1)>ncol(aBM)) {
        multform=paste("~", paste(colnames(modBF),collapse="*"))
        mBM=scale(model.matrix(as.formula(multform),modBF), center = TRUE, scale = FALSE )
        Q=qr(mBM)
        mBM=mBM[,Q$pivot[1:Q$rank],drop=FALSE]
        mBM = qr.Q(qr(mBM)) # orthogonal basis
        NS=NonSingular(TF,TM,mBM,mainblocks)
          if (!is.null(NS$TM) & !is.null(NS$TF) ) { 
            TM=NS$TM
            TF=NS$TF
            V=chol2inv(chol(crossprod(cbind(TM,mBM))))
            mVTT = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
            mVBB = V[ (ncol(TM)+1):ncol(V) , (ncol(TM)+1):ncol(V),drop=FALSE]
            mVTB = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
          } else  mBM=NULL
        } else mBM=NULL
      V=chol2inv(chol(crossprod(cbind(TM,aBM))))
      aVTT = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
      aVBB = V[ (ncol(TM)+1):ncol(V) , (ncol(TM)+1):ncol(V),drop=FALSE]
      aVTB = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
      globrelD=0
      relD=1
      globTM=TM
      globTF=TF
      for (r in 1:searches) {
        dmax=DMax(TF,TM,aVTT,aVBB,aVTB,aBM,mVTT,mVBB,mVTB,mBM,mainblocks)
        if (dmax$locrelD>(1+tol)) {
          relD=relD*dmax$locrelD
          TM=dmax$TM
          TF=dmax$TF
          aVTT=dmax$aVTT
          aVBB=dmax$aVBB
          aVTB=dmax$aVTB
          mVTT=dmax$mVTT
          mVBB=dmax$mVBB
          mVTB=dmax$mVTB
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
              z= seq_len(nunits)[mainblocks==mainblocks[s1] & BF[,i]!=BF[s1,i] & available]  
              if (length(z)==0 & counter<501) next 
              else if (length(z)==0 & counter>500) break
              if (length(z)>1) s=c(s1,sample(z,1)) else s=c(s1,z)
              aDswap=dMat(TM,aBM,aVTT,aVBB,aVTB,s)[2,1]
              if (aDswap<tol & counter<501) next
              else if (aDswap<tol & counter>500) break
              if (!is.null(mBM))  {
                mDswap=dMat(TM,mBM,mVTT,mVBB,mVTB,s)[2,1]
                if (mDswap<tol & counter<501) next
                else if (mDswap<tol & counter>500) break
                Dswap=(1-weighting)*aDswap+weighting*mDswap
              } else Dswap=aDswap
              break
            }
            if (counter>500) break
            relD=relD*Dswap 
            t=TM[s[1],]-TM[s[2],] 
            ab= aBM[s[2],]-aBM[s[1],]
            up=UpDate(aVTT,aVBB,aVTB,t,ab)
            aVTT=up$VTT
            aVBB=up$VBB
            aVTB=up$VTB
            if (!is.null(mBM)) { 
              ab= mBM[s[2],]-mBM[s[1],]
              up=UpDate(mVTT,mVBB,mVTB,t,ab)
              mVTT=up$VTT
              mVBB=up$VBB
              mVTB=up$VTB
            }
          TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
          TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]  
        } #jumps
      } #searches
      TM=globTM
      TF=globTF 
      
      H=diag(rep(1,nrow(aBM)))-tcrossprod(tcrossprod(aBM,solve(crossprod(aBM))),aBM)
      Dmodel=determinant(crossprod(t(crossprod(TM,H)),TM),logarithm = TRUE)$modulus/ncol(TM)
      Effics=rbind(Effics,c(addform,ncol(aBM),round(exp(Dmodel-Dbase),6)))

      if(!is.null(mBM)) {
        H=diag(rep(1,nrow(mBM)))-tcrossprod(tcrossprod(mBM,solve(crossprod(mBM))),mBM)
        Dmodel=determinant(crossprod(t(crossprod(TM,H)),TM),logarithm = TRUE)$modulus/ncol(TM)
        Effics=rbind(Effics,c(multform,ncol(mBM),round(exp(Dmodel-Dbase),6)))
      }
    }
    Effics=setNames(data.frame(Effics), c("Blocks model", "Blocks df","D-efficiency"))
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
  factorial=function(TF,treatments_model ) {
    gDfrac=0
    nrowTF=nrow(TF)
    replicates=nunits/nrowTF
    baseunits=nrowTF*floor(replicates) 
    TM=model.matrix(as.formula(treatments_model),TF)
    TF=TF[rep(1:nrowTF,ceiling(replicates)),,drop=FALSE]
    TM=TM[rep(1:nrowTF,ceiling(replicates)),,drop=FALSE]
    if ( replicates!=round(replicates) ) 
      Dfull=exp( determinant(crossprod(TM),logarithm = TRUE)$modulus[1]/ncol(TM) )/nrow(TM)
    for (i in 1:searches) {
      for (z in seq_len(1000)) {
        rerand=unlist(lapply(1:ceiling(replicates),function(j){(j-1)*nrowTF + sample(nrowTF)}))
        if (qr(TM[rerand[1:nunits],,drop=FALSE])$rank==ncol(TM)) break
      }
      if (z>999) stop("Unable to find a non-singular treatment start design of the required size by random search ")
      TM=TM[rerand,,drop=FALSE]
      TF=TF[rerand,,drop=FALSE]
      if ( replicates==round(replicates) ) {
        gTM=TM[1:nunits,,drop=FALSE]
        gTF=TF[1:nunits,,drop=FALSE]
        gfracEff=1
        break #quit searches
        } else {
        V=chol2inv(chol(crossprod(TM[1:nunits,,drop=FALSE]))) # V is the variance matrix of the the non-singular starting design
        repeat {
          M1VM2=       tcrossprod(tcrossprod(TM[(baseunits+1):nunits,,drop=FALSE],V),TM[(nunits+1):nrow(TM),,drop=FALSE])
          M1VM1=1-diag(tcrossprod(tcrossprod(TM[(baseunits+1):nunits,,drop=FALSE],V),TM[(baseunits+1):nunits,,drop=FALSE]))
          M2VM2=1+diag(tcrossprod(tcrossprod(TM[(nunits+1):nrow(TM),,drop=FALSE],V),TM[(nunits+1):nrow(TM),,drop=FALSE]))
          Z=M1VM2**2 + tcrossprod(M1VM1,M2VM2)
          z=which(Z == max(Z), arr.ind = TRUE)[1,]
          if (  Z[z[1],z[2]]<(1+tol)) break
          V=fractUpDate(V ,TM[baseunits+z[1],], TM[nunits+z[2],] ) # parameters(V,row_swappedout,row_swappedin) 
          TM[ c(baseunits+z[1],nunits+z[2]), ] = TM[ c(nunits+z[2],baseunits+z[1]), ]
          TF[ c(baseunits+z[1],nunits+z[2]), ] = TF[ c(nunits+z[2],baseunits+z[1]), ]
        }
          Dfrac=exp( determinant(crossprod(TM[(1:nunits),,drop=FALSE]),logarithm = TRUE)$modulus[1]/ncol(TM))/nunits
          if (Dfrac>(gDfrac+tol)) {
          gTM=TM[(1:nunits),,drop=FALSE]
          gTF=TF[(1:nunits),,drop=FALSE]
          gDfrac=Dfrac
          gfracEff=(gDfrac/Dfull)
          }
      }
    }
    return(list(TF=gTF,TM=gTM,eff=gfracEff,fraction=replicates))
  }
  # *************************************************************************************************************************
  # Main design program tests, optimizes treatment design then optimizes block design sequentially for each added block factor
  # *************************************************************************************************************************

  Z=factorial(TF,treatments_model)
  TM=scale(Z$TM, center = TRUE, scale = FALSE)[,-1]
  TF=Z$TF
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  treatsModel=data.frame("Treatment model"=treatments_model,"Treatment fraction"=Z$fraction,"D-Efficiency"=Z$eff)

  if (max(sapply(BF, nlevels))>1) {
    Opt=blocksOpt(TF,TM,BF)
    TM=Opt$TM
    TF=Opt$TF
    blocksModel=Opt$Effics
  } else {
    blocksModel=NULL
    weighting=0
  }
  Design=cbind(BF, TF)
  rownames(Design)=NULL

  list(design=Design,treatments_model=treatsModel,blocks_model=blocksModel,weighting=weighting,seed=seed,searches=searches,jumps=jumps)
  }
  