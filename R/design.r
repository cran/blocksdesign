#' @title General block and treatment designs.
#'
#' @description
#' Constructs D-optimal block and treatment designs for feasible combinations of nested or crossed block 
#' factors and feasible linear treatment models.   
#' 
#' @param treatments a single treatment factor or a data frame containing one or more qualitative or 
#' quantitative (numeric) level treatment factors.  
#' 
#' @param blocks a single block factor or a data frame containing one or more qualitative level block
#' factors in the required order of fitting.
#' 
#' @param treatments_model a model formula for the required treatments design where the default 
#' formula assumes a fully crossed factorial treatment model.
#' 
#' @param weighting a weighting factor between 0 and 1 for weighting the second-order interaction effects of
#' crossed blocks designs where the default weighting is 0.5  
#'
#' @param seed an integer initializing the random number generator. The null default gives an arbitrary 
#' random initialization.
#'
#' @param searches the maximum number of local optima searched at each stage of an
#' optimization. The default depends on the design size.
#'
#' @param jumps the number of pairwise random treatment swaps used to escape a local maxima. 
#' The default is a single swap.
#' 
#' @param fix_treatments logical TRUE (default) indicates that the full candidate set is the required
#' treatment set whereas logical FALSE indicates that the required treatment set is obtained by
#' sampling with replacement. Only applicable when \code{treatments} and \code{blocks} are of equal length.
#'
#' @details
#' 
#' The \code{treatments} object is a factor or a data frame containing one or more qualitative or quantitative
#' level factors and is the candidate set from which the final treatment design is selected. 
#' 
#' The \code{blocks} object is a factor or a data frame containing one or more qualitative level block factors where
#' the length of the factors defines the overall size of the required design.
#' The \code{blocks} object must be defined even when it is just a single factor with a single factor level.   
#' 
#' The \code{blocks} object size divided by the \code{treatments} object size is the fractional size of the required 
#' \code{treatments} design. If the fractional treatment design size is nonunity, the \code{treatments} design is
#' always chosen from the candidate set by sampling with replacement. If the fractional treatment design
#' size is unity, the method of choosing the treatment design can be determined by the \code{fix_treatments} parameter
#' where TRUE (the default) means that the whole candidate set is used while FALSE means that the treatment set is
#' chosen by sampling with  replacement.
#' 
#' The design criterion is the ratio of the generalized variance of the
#' full treatment candidate set relative to the generalized variance of the optimized treatment set 
#' for the required treatment design (D-optimality). If the required design is a fractional factorial and the 
#' candidate set is a full factorial, the candidate set will be orthogonal and any design selected from the candidate set will
#' have a relative efficiency less than or equal to 1. For a quantitative level treatment model, however, 
#' a full factorial design may not provide an optimal design and, in that case, the relative efficiency of 
#' the optimized design may well exceed 1.
#'  
#' For unstructured treatment designs, the A-efficiency factor is also shown together with an estimated A-efficiency 
#' upper-bound, where available. 
#'    
#' The \code{design} algorithm fits the blocks design by sequentially adding \code{blocks} factors
#' in the column order of the \code{blocks} data frame. Each block factor is optimized  
#' conditional on all preceding block factors remaining fixed but ignoring all succeeding block factors.
#' This method of sequential optimization allows the blocking factors to be fitted in order of importance with the 
#' largest and most important blocks fitted first and the smaller and less important blocks fitted subsequently. 
#' 
#' For crossed blocks designs, the algorithm applies a differential weighting w to determine the relative importance
#' of the blocks main effects versus the blocks interaction effects. If w = 0 the algorithm fits a simple additive
#' main effects design whereas if w = 1 the algorithm fits a fully crossed blocks design. For intermediate 0 < w < 1, 
#' the block factor interaction effects are downweighted relative to the main effects where 
#' the smaller the value of w, the greater the downweighting. The default weighting is 0.5 and 
#' provided that all block effects are estimable, this weighting will give a design where 
#' main block effects are assumed to be of greater importance than block interaction effects.
#' 
#' For example, a design for 4 replicates of 12 treatments arranged in 4 main rows and 4 main columns with
#' 3 sub-columns nested within each main column (see \code{examples}) is known to have an optimal Trojan 
#' solution with orthogonal main rows, orthogonal main columns and nested sub-columns with A-efficiency 22/31. 
#' The default weighting
#' of 0.5 will find an optimal Trojan design whereas a weighting of w = 0 will find an optimal main column 
#' blocks design with sub-optimal sub-column blocks while a weighting of w = 1 will find an optimal sub-column 
#' blocks design with sub-optimal main column blocks.  
#' 
#' Trojan designs are rare and normally it will not be possible to optimise a design for main 
#' rows, main columns and rows-by-columns interaction effects simultaneously. In that situation, a suitable
#' choice of weighting parameter can help to find a good compromise design that will give good efficiency on the main
#' effects of rows and columns and on the interaction effects of rows and columns simultaneously. 
#' 
#'   
#' Outputs:
#'
#' The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the replication of each individual treatment taken in a standard order. \cr
#'  \item  A data frame showing the randomized allocation of treatments to blocks. \cr
#'  \item  A table showing the fractional size of the treatment design and the D-efficiency factor of
#'   that fraction. \cr
#'  \item  A table showing the blocks sub-model design and the D-efficiency factor of each successively
#'   fitted blocks sub-model. \cr
#' }
#' 
#' @return
#' \item{treatments}{The treatments included in the design and the replication of each individual 
#' treatment}
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
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' DURBAN, M., HACKETT, C., MCNICOL, J., NEWTON, A., THOMAS, W., & CURRIE, I. (2003). The practical use of semi-parametric models
#'  in field trials, Journal of Agric. Biological and Envir. Stats., 8, 48-66.
#'  
#' @examples
#' 
#' ## For optimum results, the number of searches may need to be increased in practice.
#' 
#' ## 4 replicates of 12 treatments with 16 nested blocks of size 3
#' # rectangular lattice see Plan 10.10 Cochran and Cox 1957.
#' treatments = factor(rep(1:12,4))
#' blocks = data.frame(Main = gl(4,12), Sub = gl(16,3))
#' design(treatments,blocks)$blocks_model
#' 
#' ## 4 x 12 design for 4 replicates of 12 treatments with 16 nested column blocks of size 3
#' ## only an intermediate weighting will ensure an optimal Trojan design 
#' treatments = factor(rep(1:12,4))
#' blocks = data.frame(Rows = gl(4,12), Cols = gl(4,3,48), subCols = gl(12,1,48))
#' \donttest{design(treatments,blocks,searches=200)$blocks_model}
#' 
#' ## 4 x 13 Row-and-column design for 4 replicates of 13 treatments 
#' ## Youden design Plan 13.5 Cochran and Cox (1957).
#' treatments = factor(rep(1:13,4))
#' blocks = data.frame(Rows = gl(4,13), Cols = gl(13,1,52))
#' \donttest{design(treatments,blocks,searches = 700)}
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
#' \donttest{design(treatments,blocks,searches=1)$blocks_model
#' ## Finds post-blocked efficiency factors of original design; Durban et al (2003)
#' blockEfficiencies(treatments,blocks)
#' } 
#' 
#' ## differential replication 
#' treatments=factor(c(rep(1:12,2), rep(13,12)))
#' blocks = data.frame(Main = gl(2,18),  Sub = gl(12,3,36))
#' design(treatments,blocks,searches = 5)
#' 
#' ## 48 treatments in 2 replicate blocks of size 48 for a 24 x 4 array 
#' ## with 2 main rows and 3 main columns the cols factor must precede 
#' ## the rows factor otherwise the design will confound one treatment contrast
#' ## in the replicates.rows x columns interactions due to inherent aliasing 
#' treatments=factor(rep(1:48,2))
#' blocks = data.frame(Reps = gl(2,48),Cols = gl(3,8,96),Rows = gl(2,24,96))
#' design(treatments,blocks,searches=5)
#' 
#' ## Factorial treatment designs defined by a factorial model equation.
#' 
#' ## Main effects of five 2-level factors in a half-fraction of 
#' ## a 2/2/2 nested blocks design
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:2), F3 = factor(1:2), 
#' F4 = factor(1:2), F5 = factor(1:2))
#' blocks = data.frame(b1 = gl(2,8),b2 = gl(4,4),b3 = gl(8,2))
#' model="F1 + F2 + F3 + F4 + F5"
#' \donttest{ repeat {
#'  z = design(treatments,blocks,treatments_model=model,searches=5)
#'  if ( z$blocks_model[3,3] == 1) break }
#'  print(z)}
#' 
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:2), F3 = factor(1:2), 
#' F4 = factor(1:2), F5 = factor(1:2))
#' blocks = data.frame(blocks = gl(4,8))
#' model = "(F1 + F2 + F3 + F4 + F5)^2"
#' design(treatments,blocks,treatments_model=model,searches = 10)
#' 
#' # Main effects of five 2-level factors in a half-fraction of 
#' # a 4 x 4 row-and column design.
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:2), F3 = factor(1:2), 
#' F4 = factor(1:2), F5 = factor(1:2))
#' blocks = data.frame(rows = gl(4,4), cols = gl(4,1,16))
#' model = "~ F1 + F2 + F3 + F4 + F5"
#' design(treatments,blocks,treatments_model=model,searches = 50)
#' 
#' # Quadratic regression for three 3-level numeric factor assuming a 10/27 fraction
#' treatments = expand.grid(A = 1:3, B = 1:3, C = 1:3)
#' blocks=data.frame(main=gl(1,10))
#' model = " ~ ( A + B + C)^2 + I(A^2) + I(B^2) + I(C^2)"
#' design(treatments,blocks,treatments_model=model,searches=5) 
#' 
#' # First-order model for 1/3rd fraction of four qualitative 3-level factors in 3  blocks
#' treatments = expand.grid(F1 = factor(1:3), F2 = factor(1:3), F3 = factor(1:3), 
#' F4 = factor(1:3))
#' blocks = data.frame(main = gl(3,9))
#' model = " ~ F1 + F2 + F3 + F4"
#' design(treatments,blocks,treatments_model=model,searches=25)
#' 
#' # Second-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' treatments = expand.grid(F1 = factor(1:3), F2 = factor(1:3), F3 = factor(1:3), 
#' F4 = factor(1:3), F5 = factor(1:3))
#' blocks=data.frame(main=gl(3,27))
#' model = " ~ (F1 + F2 + F3 + F4 + F5)^2"
#' \donttest{design(treatments,blocks,treatments_model=model,searches=500)}
#' 
#' # Second-order model for two qualitative and two quantitative level factors in 4 blocks
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:3), V1 = 1:3, V2 = 1:4)
#' blocks = data.frame(main = gl(4,18))
#' model = " ~ F1 + F2 + poly(V1,2) + poly(V2,2) + (poly(V1,1) + F1 + F2):(poly(V2,1) + F1 + F2)"
#' \donttest{design(treatments,blocks,treatments_model=model,searches=5)}
#'  
#' # Plackett and Burman design for eleven 2-level factors in 12 runs 
#' GF = expand.grid(F1 = factor(1:2,labels=c("a","b")), F2 = factor(1:2,labels=c("a","b")), 
#'                  F3 = factor(1:2,labels=c("a","b")), F4 = factor(1:2,labels=c("a","b")),
#'                  F5 = factor(1:2,labels=c("a","b")), F6 = factor(1:2,labels=c("a","b")),
#'                  F7 = factor(1:2,labels=c("a","b")), F8 = factor(1:2,labels=c("a","b")), 
#'                  F9 = factor(1:2,labels=c("a","b")), F10= factor(1:2,labels=c("a","b")), 
#'                  F11= factor(1:2,labels=c("a","b")) )
#' blocks=data.frame(main=gl(1,12))
#' model = "~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9 + F10 + F11"
#' \donttest{design(GF,blocks,treatments_model=model,searches=25)}
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames 
#' @importFrom lme4 lmer
#' @importFrom plyr count
#' 
  design = function(treatments,blocks=NULL,treatments_model=NULL,fix_treatments=TRUE,
                    weighting=0.5,searches=NULL,seed=NULL,jumps=1) {
  options(contrasts=c('contr.SAS','contr.poly'))
  tol = .Machine$double.eps ^ 0.5
  if (!is.null(seed)) set.seed(seed)
  
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
      if (nrow(TF)>(1+rank))
        s1=sample(pivot[(1+rank):nrow(TF)],1)
      else s1=nrow(TF)
        swaptrts=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
        candidates = (1:nrow(TF))[swaptrts & mainblocks==mainblocks[s1]]
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
    if (ncol(TM)+ncol(BM)>(nrow(TF)-1)) stop(" Block design has too many parameters")
    fullrank=ncol(TM)+ncol(BM) # fullrank not more than (plots-1)
    Q=qr(t(cbind(TM,BM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:500) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      swap=Swaps(TF,pivot,rank,mainblocks)
      TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
      Q=qr(t(cbind(TM,BM)))
      if (Q$rank>rank) {
        TF[c(swap[1],swap[2]),]=TF[c(swap[2],swap[1]),]
        rank=Q$rank
        pivot=Q$pivot
      } else TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
    }
    stop(" Unable to find an initial non-singular design - try different initial seed or simplify the block design")
  }
  # *********************************************************************************************************************
  #  dMat is a matrix where the ith, jth matrix element is the change in the D-max criterion of the block design 
  # due to swapping ith and jth treatments of a sample s of the set of plots 
  # *********************************************************************************************************************
  dMat=function(TM,BM,V,s,dim0) {
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
  DMax=function(TF,TM,BM,V,mainBlocks,dim0) {
    locrelD=1
    mainSets=tabulate(mainBlocks)
    nSamp=pmin(rep(8,nlevels(mainBlocks)), mainSets)
    repeat {
      kmax=1
      for (k in 1: nlevels(mainBlocks)) {
        s=sort(sample((1:nrow(TF))[mainBlocks==levels(mainBlocks)[k]] , nSamp[k])) 
        dMat=dMat(TM,BM,V,s,dim0)
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
      }  else if (sum(nSamp) == nrow(TF)) break
      else nSamp=pmin(mainSets,2*nSamp)
    }
    list(V=V,locrelD=locrelD,TF=TF,TM=TM)
  }
  # ************************************************************************************************************************
  # Optimize the nested blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ************************************************************************************************************************
   blocksOpt=function(TF,TM,BF,searches,nested,treatments_model) {
    dim0 = ncol(TM)
    mainBlocks = data.frame(factor(rep(1,nrow(BF))), matrix(nrow=nrow(BF),ncol=ncol(BF)))
    for (i in 1: ncol(BF)) mainBlocks[,i+1] = droplevels(interaction( mainBlocks[,i],BF[,i]))
    colnames(mainBlocks) = c("mean", unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse=".")})) )
    for (u in 1:ncol(BF)) {
      BM1 = scale(model.matrix(as.formula(paste("~",paste0(colnames(BF)[1:u],collapse="+"))),BF), center = TRUE, scale = FALSE)
      Q = qr(BM1)
      BM1 = BM1[,Q$pivot[1:Q$rank],drop=FALSE] # full rank
      if ((1+dim0+ncol(BM1))>nrow(TF)) stop("Additive block design has too many fixed effects for number of plots ")
      if(u>1 & weighting>0){
        BM2=scale(model.matrix(as.formula(paste("~",paste0("(",paste0(colnames(BF)[1:u],collapse="+"),")^2"))),BF), center = TRUE, scale = FALSE)
        Q=qr(BM2)
        BM2 = BM2[,Q$pivot[1:Q$rank],drop=FALSE] # full rank
      } else BM2=BM1
      deg2=(weighting>0 & ncol(BM2)>ncol(BM1) & (dim0 + ncol(BM2))<nrow(TF))
      if (deg2) BM = BM2 else BM = BM1
      
      BM = qr.Q(qr(BM)) # orthogonal basis for BM 
      NS = NonSingular(TF,TM,BM,mainBlocks[,u],searches)
      TM = NS$TM
      TF = NS$TF  
      globrelD = 0
      relD = 1
      globTM = TM
      globTF = TF
      if (deg2) BM[,(ncol(BM1)+1):ncol(BM2)] = BM[,(ncol(BM1)+1):ncol(BM2)]*weighting
      Info=crossprod(cbind(TM,BM))
      if (deg2) 
        Info[(1+dim0+ncol(BM1)):(dim0+ncol(BM2)),(1+dim0+ncol(BM1)):(dim0+ncol(BM2))] = 
        Info[(1+dim0+ncol(BM1)):(dim0+ncol(BM2)),(1+dim0+ncol(BM1)):(dim0+ncol(BM2))]/(weighting^2)
      V=chol2inv(chol(Info))
      for (r in 1:searches) {
        dmax=DMax(TF,TM,BM,V,mainBlocks[,u],dim0)
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
            s1=sample(seq_len(nrow(TF)),1)
            available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
            z= seq_len(nrow(TF))[mainBlocks[,u]==mainBlocks[s1,u] & mainBlocks[,u+1]!=mainBlocks[s1,u+1] & available]  
            if (length(z)==0 & counter<501) next 
            else if (length(z)==0 & counter>500) break
            if (length(z)>1) s=c(s1,sample(z,1)) else s=c(s1,z)
            testDswap=dMat(TM,BM,V,s,dim0)[2,1]
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
  # Updates variance matrix for swapped rows where mi is swapped out and qj is swapped in
  # ******************************************************************************************************************
  fractUpDate=function(V,mi,qj) {
    f=crossprod(V,qj)
    g=crossprod(V,mi)
    a=as.numeric(crossprod(qj,f))
    b=as.numeric(crossprod(mi,g))
    c=as.numeric(crossprod(mi,f))
    W=g*sqrt(1+a)-f*c/sqrt(1+a)
    U=f*sqrt(1-b+(c*c)/(1+a))
    V=V-(tcrossprod(U)-tcrossprod(W))/((1+a)*(1-b)+c*c)
    return(V)
  }
  # ******************************************************************************************************************
  # non singular starting design expands TF if necessary to contain ceiling(repeats) of the original TF then ensures
  # the first plots of TF is non-singular for the rewuired treatment design matrix TM
  # ******************************************************************************************************************
  nonSingular=function(TF,TM,nunits) {
    singular=FALSE
    for (t in 1:searches) {
      rank=qr(TM[1:nunits,])$rank
      if (identical(ncol(TM),rank)) break
      if (nunits==nrow(TF)) stop("Unable to find non-singular starting design by random search")
      for (z in seq_len(100000)) {
        p1 = sample(1:nunits,1)
        if (nrow(TF)>(nunits+1))
          p2 = sample( (nunits+1):nrow(TF) ,1)
        else p2 = nrow(TF)
        TM[c(p1,p2),] = TM[c(p2,p1),]
        TF[c(p1,p2),] = TF[c(p2,p1),]
        newrank = qr(TM[1:nunits,])$rank
        if (newrank<rank) {
          TM[c(p1,p2),] = TM[c(p2,p1),]
          TF[c(p1,p2),] = TF[c(p2,p1),]
          } else rank=newrank
      if (identical(ncol(TM),rank)) break
      }
      if (identical(ncol(TM),rank)) break
      if (t>=searches) singular=TRUE
    }
    list(TM=TM[1:nunits,,drop=FALSE],TF=TF[1:nunits,,drop=FALSE],singular=singular)
  }
  # *********************************************************************************************************************
  # Fractional factorials for quantitative level factors - allows replacement of swapped rows from candidate set
  # *******************************************************************************************************************
  factorial = function(TF,BF,treatments_model,fix_treatments) {
    if (fix_treatments & nrow(TF)==nrow(BF))   
      return(list(TF = TF[sample(1:nrow(TF)), ,drop=FALSE], Eff = 1 ,fraction = 1))
    TM=model.matrix(as.formula(treatments_model),TF)
    if (ncol(TM)>nrow(BF)) stop("too many treatment parameters to be estimated for the given number of experimental units")
    fraction=nrow(BF)/nrow(TF)
    Dmax = exp( (determinant(crossprod(TM),logarithm = TRUE)$modulus)/ncol(TM))/nrow(TM)
    Dfacteff=1
    gDfrac = 0
    for (i in 1:searches) {
      fullrand = unlist(lapply(1:ceiling(fraction),function(j){sample(seq_len(nrow(TF)))}))
      newTF = TF[fullrand,,drop=FALSE]
      newTM = TM[fullrand,,drop=FALSE]
      W = nonSingular(newTF,newTM,nrow(BF)) 
      if (W$singular==TRUE)  stop("Unable to find non-singular starting design for factorial treatment design by random search")
      newTM = W$TM
      newTF = W$TF
      Dfacteff = NA
      V = chol2inv(chol(crossprod(newTM))) # V is the variance matrix of the the non-singular starting design
      counter=0
      repeat {
        s=sample(1:nrow(TF),min(1024,nrow(TF)))
        M1VM2 =       tcrossprod(tcrossprod(  newTM,V),TM[s,])
        M1VM1 = 1-diag(tcrossprod(tcrossprod( newTM,V),newTM))
        M2VM2 = 1+diag(tcrossprod(tcrossprod( TM[s,],V),TM[s,]))
        Z = M1VM2**2 + tcrossprod(M1VM1,M2VM2)
        z = which(Z == max(Z), arr.ind = TRUE)[1,]
        improved=(Z[z[1],z[2]]>(1+tol))
        if (improved) counter=0
        else counter=counter+1
        if ( counter==9 | ( counter==1 & length(s)==nrow(TF))  )  break
        if (!improved) next
        V = fractUpDate(V,newTM[z[1],],TM[s[z[2]],])   # parameters(V,row_swappedout,row_swappedin) 
        newTM[z[1],] = TM[s[z[2]],]
        newTF[z[1],] = TF[s[z[2]],]
      }
      Dfrac = exp((determinant(crossprod(newTM[1:nrow(BF),,drop=FALSE]),logarithm = TRUE)$modulus)/ncol(newTM))/nrow(BF)
      if (Dfrac>gDfrac) {
        gTM=newTM
        gTF=newTF
        gDfrac=Dfrac
      }
      if (isTRUE(all.equal(gDfrac,Dmax,tolerance = tol))) break
    }
    TF=gTF[1:nrow(BF),,drop=FALSE]
    return(list(TF = TF, Eff = round(gDfrac/Dmax,7) ,fraction = fraction))
  }
  # *************************************************************************************************************************
  # Main design program tests, optimizes treatment design then optimizes block design sequentially for each added block factor
  # *************************************************************************************************************************
  TF=data.frame(treatments)
  if (is.null(blocks)) blocks= gl(1,nrow(TF))
  BF=data.frame(blocks)
  
  if (!all(sapply(BF,is.factor))) stop("blocks must be factors")
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) stop(" max. jumps is 10")

  if (is.null(searches)) 
    if (nrow(BF)<100) 
      searches=50 else if (nrow(BF)<1000) 
        searches=5000%/%nrow(BF) else if (nrow(BF)<5000) 
          searches=5000%/%nrow(BF) else 
            searches=1
  if (!is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  
  if (is.null(treatments_model)) treatments_model=paste("~",paste(colnames(TF),collapse="*"))
  if (!substring(trimws(treatments_model, "l"),1,1)=="~") treatments_model=paste("~",treatments_model)
  
  Z=factorial(TF,BF,treatments_model,fix_treatments)
  TF=data.frame(Z$TF)
  treatsModel=data.frame(cbind("Treatment model" = treatments_model,"D-Efficiency" = Z$Eff))
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  TM = model.matrix(as.formula(treatments_model),TF)
  TM=scale(TM, center = TRUE, scale = FALSE)[,-1,drop=FALSE]
  nested=all(sapply(1:ncol(BF),function(i){nlevels(droplevels(interaction(BF[,c(1:i)])))==nlevels(BF[,i])}))
 # for nested equi-block designs with a single unstructured treatment factor uses blocks() function

  bf = table(BF[,ncol(BF)])
  if (nested  & max(bf)==min(bf) & ncol(TF)==1 & is.factor(TF[,1])) {
    tlevs=levels(TF[,1])
    blevs=c(1,sapply(BF, nlevels))
    blevs=sapply(2:length(blevs), function(i){blevs[i]/blevs[i-1]})
    Z=blocks(rep(1,nlevels(TF[,1])),table(TF[,1]),blevs,searches=searches,seed=seed,jumps=1)
    Treatments=Z$Treatments
    levels(Treatments[,1]) = tlevs
    design=Z$Design
    levels(design[,ncol(design)]) = tlevs
    blocksModel=Z$blocks_model
    weighting=NULL
  } else {
    if (!is.numeric(weighting)) stop("weighting must be a number between 0 and 1")
    if (weighting<0 | weighting>1) stop("weighting must be a number between 0 and 1")
    if (max(sapply(BF, nlevels))>1) {
      Opt=blocksOpt(TF,TM,BF,searches,nested,treatments_model)
      TM=Opt$TM
      TF=Opt$TF
      blocksModel=Opt$Effics
    } else blocksModel=NULL
    design=cbind(BF,TF)
    Treatments=count(TF,colnames(TF))
  }

  rownames(design)=NULL
  rownames(Treatments)=NULL
  list(treatments=Treatments,design=design,treatments_model=treatsModel,blocks_model=blocksModel,
       weighting=weighting,seed=seed,searches=searches,jumps=jumps)
  }
  