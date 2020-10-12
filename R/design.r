#' @title General block and treatment designs.
#'
#' @description
#' Constructs general D-optimal designs for any feasible linear treatment model with any feasible combinations of block 
#' factors.
#' 
#' @param treatments a single treatment factor or data frame for the candidate set 
#' of treatment factor combinations assuming any combination of qualitative or quantitative factor levels.
#' 
#' @param blocks a single block factor or data frame for the required combinations of 
#' block factors in the required order of fitting assuming quantitative block factor levels only.
#' 
#' @param treatments_model a character vector containing one or more nested treatment model formula. 
#' 
#' @param weighting a weighting factor between 0 and 1 for weighting the 2-factor interaction effects of
#' factorial blocks.
#'
#' @param seed an integer initializing the random number generator.
#'
#' @param searches the maximum number of local optima searched at each stage of an
#' optimization.
#'
#' @param jumps the number of pairwise random treatment swaps used to escape a local maxima. 
#' 
#' 
#' @details
#' 
#' \code{treatments} is a factor or data frame containing one or more qualitative
#' or quantitative level treatment factors defining the candidate treatment set. The 
#' required treatment design is selected from the candidate treatment set without replacement.
#' 
#'  i) If the candidate set is the same size as the required design, the full candidate set will be
#' selected. 
#' 
#' ii) If the candidate set is larger than the required design, a subset of the candidate set 
#' will be selected by optimizing the treatment D-optimality criterion. 
#' 
#' iii) If the candidate set is smaller than the required design size, the candidate set will be replicated 
#' until the candidate set equals or exceeds the required design size.
#' 
#' \code{blocks} is a factor or data frame containing one or more qualitative level block factors, added
#' sequentially. If the blocks are fully nested or if the design has additive crossed blocks,
#' each added block factor is optimized by swaps made within the levels of all previously added blocks. If,
#' however, a design has crossed blocks with estimable interaction effects, the relative importance of the main
#' factorial block effects versus the 2-factor block interaction effects is weighted by a parameter w in the range 0 to 1. 
#' 
#' i) If w = 0, the design is a simple additive block effects. 
#' 
#' ii) If w = 1, the design is a simple interaction blocks design for the crossed blocks interaction effects
#' 
#' iii) If 0 < w < 1, the 2-factor block interaction effects are weighted relative to the additive block effects with
#' the importance of the interaction effects assumed to increase as w increases from 0 to 1. 
#' 
#' The length of the \code{blocks} object defines the total number of plots in the design. 
#'  
#' \code{treatments_model} is a character vector containing one or more nested treatments formula where
#' each model formula, taken in order, is optimized with the treatment factors of all previously
#' optimized model formula remaining constant. Sequential model fitting provides improved flexibility
#' for fitting factors or variables of different status or importance (see examples below).
#' 
#' The treatment design criterion for each treatment model is the generalized variance of the information matrix for that design 
#' and the design efficiency is the ratio of the generalized variance of the full candidate treatment model
#' relative to the generalized variance of the optimized design. The efficiency will be less than or equal to 1 for factorial
#' models but may exceed 1 for polynomial models. 
#' 
#' For more details see \code{vignette(package = "blocksdesign")}  
#' 
#' @return
#' \item{Treatments}{The treatments included in the design and the replication of each individual 
#'  treatment taken in de-randomized standard order.}
#' \item{Design}{The design layout showing the randomized allocation of treatments to blocks and plots.}
#' \item{Treatments_model}{The fitted treatment model, the number of model parameters (DF)
#'   and the D-efficiency of each sequentially fitted treatment model}
#' \item{Blocks_model}{The blocks sub-model design and 
#'  the D- and A-efficiency factors of each successively fitted sub-blocks model.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#' 
#' @references
#' 
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' @examples
#' ## For optimum results, the number of searches may need to be increased.
#' 
#' ## 4 replicates of 12 treatments with 16 nested blocks of size 3
#' ## rectangular lattice see Plan 10.10 Cochran and Cox 1957.
#' blocks = data.frame(Main = gl(4,12), Sub = gl(16,3))
#' design(treatments = gl(12,1,48), blocks)
#' 
#' ## 3 replicates of 15 treatments in 3 main blocks with 3 nested blocks and one control treatment
#' blocks=data.frame( Main = gl(3,18,54),Sub = gl(9,6,54))
#' treatments=factor(rep(c(1:15,rep("control",3)),3),levels = c(1:15,"control") )
#' Z=design(treatments,blocks)
#' incid=table(interaction(Z$Design$Main,Z$Design$Sub,lex.order = TRUE),Z$Design$treatments)
#' Z # print design
#' incid # print incidences of treatments in blocks
#' crossprod(incid) # print pairwise concurrences within blocks
#'  
#' ## 4 x 12 design for 4 replicates of 12 treatments with 3 plots in each intersection block
#' ## The optimal design is Trojan with known A-efficiency = 22/31 for the intersection blocks
#' blocks = data.frame(Rows = gl(4,12), Cols = gl(4,3,48))
#' Z=design(treatments =gl(12,1,48),blocks)
#' incid=table(interaction(Z$Design$Rows,Z$Design$Cols,lex.order = TRUE),Z$Design$treatments)
#' Z # print design
#' incid # print incidences of treatments in blocks
#' crossprod(incid) # print pairwise concurrences within blocks
#' ## as above but showing 3 sub-columns nested within each main column
#' blocks = data.frame(Rows = gl(4,12), Cols = gl(4,3,48), subCols = gl(12,1,48))
#' \donttest{Z=design(treatments = gl(12,1,48),blocks,searches=200)
#' Z # print design}
#' 
#' ## 4 x 13 Row-and-column design for 4 replicates of 13 treatments 
#' ## Youden design Plan 13.5 Cochran and Cox (1957).
#' blocks = data.frame(Rows = gl(4,13), Cols = gl(13,1,52))
#' \donttest{Z=design(treatments = gl(13,1,52),blocks,searches = 700)
#' incid=table(Z$Design$Cols,Z$Design$treatments)
#' Z # print design
#' crossprod(incid) # print pairwise concurrences of treatments within column blocks (BIB's)
#' tcrossprod(incid) # print pairwise concurrences of column blocks within treatments (Dual design)}
#' 
#' ## 48 treatments in 2 replicate blocks with 2 nested rows in each replicate and 3 main columns
#' ##  (Reps/Rows) x Cols
#' blocks = data.frame(Reps = gl(2,48), Rows = gl(4,24,96), Cols = gl(3,8,96))
#' design(treatments=gl(48,1,96),blocks,searches=5)
#' 
#' ## 48 treatments in 2 replicate blocks with 2 main columns
#' ## The default weighting gives non-estimable Reps:Cols effects due to inherent aliasing
#' ## Increased weighting gives estimable Reps:Cols effects but non-orthogonal main effects
#' blocks = data.frame(Reps = gl(2,48), Cols = gl(2,24,96))
#' design(treatments=gl(48,1,96),blocks,searches=5)
#' design(treatments=gl(48,1,96),blocks,searches=5,weighting=.9)
#' 
#' ## Factorial treatment designs defined by a single factorial treatment model
#' 
#' ## Main effects of five 2-level factors in a half-fraction in 2/2/2 nested blocks design 
#' ## (may require 100's of repeats to find a fully orthogonal solution - a VERY long wait!)
#' treatments = expand.grid(F1 = gl(2,1), F2 = gl(2,1),F3 = gl(2,1), F4 = gl(2,1), F5 = gl(2,1))
#' blocks = data.frame(b1 = gl(2,8),b2 = gl(4,4),b3 = gl(8,2))
#' model=" ~ F1 + F2 + F3 + F4 + F5"
#' \donttest{repeat {z = design(treatments,blocks,treatments_model=model,searches=50)
#' if ( isTRUE(all.equal(z$Blocks_model[3,3],1) ) ) break }
#' print(z)}
#'  
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' treatments = expand.grid(F1 = gl(2,1), F2 = gl(2,1),F3 = gl(2,1), F4 = gl(2,1), F5 = gl(2,1))
#' blocks = data.frame(blocks = gl(4,8))
#' model = " ~ (F1 + F2 + F3 + F4 + F5)^2"
#' design(treatments,blocks,treatments_model=model,searches = 10)
#' 
#' # Main effects of five 2-level factors in a half-fraction of a 4 x 4 row-and column design.
#' treatments = expand.grid(F1 = gl(2,1), F2 = gl(2,1),F3 = gl(2,1), F4 = gl(2,1), F5 = gl(2,1))
#' blocks = data.frame(rows = gl(4,4), cols = gl(4,1,16))
#' model = "~ F1 + F2 + F3 + F4 + F5"
#' \donttest{repeat {z = design(treatments,blocks,treatments_model=model,searches=50)
#' if ( isTRUE(all.equal(z$Blocks_model[2,3],1) ) ) break }
#'  print(z)}
#' 
#' # Quadratic regression for three 3-level numeric factor assuming a 10/27 fraction
#' treatments = expand.grid(A = 1:3, B = 1:3, C = 1:3)
#' blocks=data.frame(main=gl(1,10))
#' model = " ~ ( A + B + C)^2 + I(A^2) + I(B^2) + I(C^2)"
#' design(treatments,blocks,treatments_model=model,searches=10) 
#' 
#' # Quadratic regression for three 3-level numeric factor crossed with a qualitative 2-level factor
#' treatments = expand.grid(F = factor(1:2), A = 1:3, B = 1:3, C = 1:3)
#' blocks=data.frame(main=gl(1,18))
#' model = " ~ F + A + B + C + F:A + F:B + F:C + A:B + A:C + B:C + I(A^2) + I(B^2) + I(C^2)"
#' design(treatments,blocks,treatments_model=model,searches=5) 
#' 
#' # 1st-order model for 1/3rd fraction of four qualitative 3-level factors in 3  blocks
#' treatments = expand.grid(F1 = gl(3,1), F2 = gl(3,1), F3 = gl(3,1), F4 = gl(3,1))
#' blocks = data.frame(main = gl(3,9))
#' model = " ~ F1 + F2 + F3 + F4"
#'\donttest{ design(treatments,blocks,treatments_model=model,searches=25)}
#' 
#' # 2nd-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' # (may require many repeats to find a fully orthogonal solution)
#' treatments = expand.grid(F1 = gl(3,1), F2 = gl(3,1),F3 = gl(3,1), F4 = gl(3,1), F5 = gl(3,1))
#' blocks=data.frame(main=gl(3,27))
#' model = " ~ (F1 + F2 + F3 + F4 + F5)^2"
#' \donttest{ repeat {z = design(treatments,blocks,treatments_model=model,searches=50)
#' if ( isTRUE(all.equal(z$Blocks_model[1,3],1) ) ) break}
#'  print(z) }
#' 
#' # 2nd-order model for two qualitative and two quantitative level factors in 2 blocks of size 18
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:3), V1 = 1:3, V2 = 1:4)
#' blocks = data.frame(main = gl(2,18))
#' model = " ~ (F1 + F2 + V1 + V2)^2 +  I(V1^2) +  I(V2^2)"
#'\donttest{design(treatments,blocks,treatments_model=model,searches=5)}
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
#' \donttest{D=design(GF,blocks,treatments_model=model,searches=25)
#' round(crossprod(scale(data.matrix(D$Design)[,-1])),6) }
#' 
#' ## Factorial treatment designs defined by sequentially fitted factorial treatment models
#' 
#' ## 4 varieties by 3 levels of N by 3 levels of K assuming degree-2 treatment
#' ## interaction effects and two blocks of 12 plots
#' ## the single stage model gives an unequal split for the replication of the four varieties
#' ## which may be undesirable whereas the two stage model forces an equal split of 6 plots
#' ## per variety. The single stage model appears slightly more efficient but
#' ## in this example global D-optimality does not give the most suitable design structure. 
#' treatments = expand.grid(Variety = factor(1:4), N = 1:3, K = 1:3)
#' blocks=data.frame(main=gl(2,12))
#' treatments_model = " ~  (Variety + N + K)^2  + I(N^2) + I(K^2)"
#' design(treatments,blocks,treatments_model=treatments_model,searches=10) 
#' treatments_model = c(" ~ Variety" ," ~ Variety + (Variety + N + K)^2 + I(N^2) + I(K^2)")
#' design(treatments,blocks,treatments_model=treatments_model,searches=10)
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames 
#' @importFrom plyr count
#' @importFrom plyr match_df
#' 
  design = function(treatments,blocks=NULL,treatments_model=NULL,weighting=0.5,searches=NULL,seed=NULL,jumps=1) {
    tol = .Machine$double.eps ^ 0.5
    options(contrasts=c('contr.SAS','contr.poly'))
    options(warn=0)
    if (!is.null(seed)) set.seed(seed) 
    
    # *********************************************************************************************************************
    #  dMat is a matrix where the ith, jth matrix element is the change in the D-max criterion of the block design 
    # due to swapping ith and jth treatments of a sample s of the set of plots 
    # *********************************************************************************************************************
    dMat=function(TM,BM,V,s) {
      sTM=TM[s,,drop=FALSE]
      sBM=BM[s,,drop=FALSE]
      TMT=crossprod(t(crossprod(t(sTM),V[1:ncol(TM),1:ncol(TM),drop=FALSE])),t(sTM))
      TMB=crossprod(t(crossprod(t(sTM),V[1:ncol(TM),(ncol(TM)+1):ncol(V), drop=FALSE])),t(sBM))
      BMB=crossprod(t(crossprod(t(sBM),V[(ncol(TM)+1):ncol(V),(ncol(TM)+1):ncol(V),drop=FALSE])),t(sBM))
      TT=TMT-diag(TMT)
      BB=BMB-diag(BMB)
      TB=TMB-diag(TMB)
      dMat=(1+TB+t(TB))**2 - (TT + t(TT))*(BB + t(BB))
      return(dMat)
    }
    
    # ******************************************************************************************************************
    # Maximises the matrix dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
    # Sampling is used initially but later a full search is used to ensure steepest ascent optimization.
    # *******************************************************************************************************************
    DMax=function(TF,TM,BM,V,mainBlocks) {
      locrelD=1
      mainSets=tabulate(mainBlocks)
      nSamp=pmin(rep(8,nlevels(mainBlocks)), mainSets)
      repeat {
        kmax=1
        for (k in 1: nlevels(mainBlocks)) {
          s=sort(sample((1:nrow(TF))[mainBlocks==levels(mainBlocks)[k]], nSamp[k])) 
          dMat=dMat(TM,BM,V,s)
          z=which(dMat == max(dMat,na.rm=TRUE), arr.ind = TRUE)[1,]
          if (dMat[z[1],z[2]]>kmax) {
            kmax=dMat[z[1],z[2],drop=FALSE]
            pi=s[z[1]]
            pj=s[z[2]]
          } 
        }
        if (kmax>(1+tol)) {
          locrelD=locrelD*kmax
          V=UpDate(V,(TM[pi,,drop=TRUE]-TM[pj,,drop=TRUE]),(BM[pj,,drop=TRUE]-BM[pi,,drop=TRUE]),ncol(TM)  )
          TM[c(pi,pj),]=TM[c(pj,pi),,drop=FALSE]
          TF[c(pi,pj),]=TF[c(pj,pi),,drop=FALSE]
        }  else if (sum(nSamp) == nrow(TF)) break
        else nSamp=pmin(mainSets,2*nSamp)
      }
      list(V=V,locrelD=locrelD,TF=TF,TM=TM)
    }
 
    # ********************************************************************************************************************
    # Random swaps for blocks optimization
    # ********************************************************************************************************************    
    Swaps=function(TF,BF,pivot,rank,restrict,blocks) {
      candidates=NULL
      counter=0
      while (length(candidates)==0 & counter<500) {
        counter=counter+1
        swapOut=pivot[1:rank]
        swapIn=pivot[(1+rank):nrow(TF)]
        if (length(swapIn)>1) 
          s1=sample(swapIn,1) else s1=swapIn
        rownames(TF)=NULL
        sameTrts=suppressMessages(as.integer(rownames(match_df(TF, TF[s1,,drop=FALSE])))) # vector of row numbers of all replications of TF[s1,]
        candidates = seq_len(nrow(TF))[restrict==restrict[s1] & blocks!=blocks[s1] & !(seq_len(nrow(TF))%in%sameTrts) & (seq_len(nrow(TF))%in%swapOut)]
      }
      if (length(candidates)==0) stop(" 1. Unable to find an initial non-singular starting design - perhaps try a simpler block design?")
      if (length(candidates)>1) s2=sample(candidates,1) else s2=candidates 
      return(c(s1,s2))
    } 
    # **********************************************************************************************************************
    # Initial randomized starting design. If the initial design is rank deficient, 
    # random swaps with positive selection are used to to increase design rank
    # **********************************************************************************************************************
    blocksNonSingular=function(TF,TM,BM,restrict,blocks) {
      Q=qr(t(cbind(TM,BM)))
      rank=Q$rank
      pivot=Q$pivot
      for ( i in 1:500 ) {
        if (rank == (ncol(TM) + ncol(BM))) 
            return(list(TF = TF, TM = TM, fullrank = TRUE))
        swap=Swaps(TF,BF,pivot,rank,restrict,blocks)
        TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),,drop=FALSE]
        Q=qr(t(cbind(TM,BM)))
        if (Q$rank>rank) {
          rank=Q$rank
          pivot=Q$pivot
          TF[c(swap[1],swap[2]),]=TF[c(swap[2],swap[1]),,drop=FALSE]
        } else 
          TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),,drop=FALSE]
      }
      return(list(TF=TF,TM=TM,fullrank=FALSE))
    }
    
    # ************************************************************************************************************************
    # Orthogonal basis for M
    # ************************************************************************************************************************
    orthogM=function(M) {
      QR = qr(M) # qr transformation
      if (QR$rank < min(nrow(M),ncol(M)))
        QR = qr(M[, QR$pivot[1:QR$rank] ,drop=FALSE]) # removes any aliased block effects then finds qr transformation
      M = qr.Q(QR) # orthogonal basis for M where Q'Q=I
      M
    }
    
    # ************************************************************************************************************************
    # Optimize the  blocks assuming a possible set of Main block constraints Initial completely randomized starting design.
    # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
    # Each factor in BF may include levels from preceding factors hence qr pivoting is essential to sweep out factor dependencies  
    # ************************************************************************************************************************
    blocksOpt=function(TF,TM,BF) {
      IBF=data.frame(factor(rep(1,nrow(TF))),lapply(1:ncol(BF), function(i) droplevels(interaction(BF[,c(1:i)]))))
      colnames(IBF)=c("mean",unlist(lapply(1:ncol(BF),function(j){paste0(colnames(BF)[1:j],collapse=".")})))
      D_Effic=rep(0,length=ncol(BF))
      A_Effic=rep(0,length=ncol(BF))
      D_IntEff=rep(0,length=ncol(BF))
      A_IntEff=rep(0,length=ncol(BF))
      Int_levs=rep(0,length=ncol(BF))
      Add_levs=rep(0,length=ncol(BF))
      for (u in 1:ncol(BF)) {
        BM1 = scale(model.matrix(as.formula(paste("~",addfactors[u])),BF), center = TRUE, scale = FALSE)[,-1,drop=FALSE]
        BM1=orthogM(BM1)
        Add_levs[u]=ncol(BM1)
         BM2 = scale(model.matrix(as.formula(paste("~",prodfactors[u])),BF), center = TRUE, scale = FALSE)[,-1,drop=FALSE] 
         BM2=orthogM(BM2)
        Int_levs[u]=ncol(BM2)
        if ((ncol(TM) + ncol(BM2)+1) > nrow(TM)) BM2=NULL
        if (!is.null(BM2)) {
          NS = blocksNonSingular(TF,TM,BM2,IBF[,u,drop=TRUE],BF[,u,drop=TRUE])
          if (NS$fullrank==FALSE) BM2=NULL
        }
        if (is.null(BM2))
          NS = blocksNonSingular(TF,TM,BM1,IBF[,u,drop=TRUE],BF[,u,drop=TRUE])
        
        if (NS$fullrank==FALSE)  stop("Cannot find a non-singular starting block design - perhaps try a simpler design? ")
        TM = NS$TM
        TF = NS$TF 
    
        if ( !is.null(BM2) & weighting>0) BM=BM2 else BM=BM1
        Info=crossprod(cbind(TM,BM))
        # down weights two factor interactions if present
        if (!is.null(BM2) && ncol(BM2)>ncol(BM1) && weighting>0)
          Info[(1+ncol(BM1)+ncol(TM)):ncol(Info),(1+ncol(TM)+ncol(BM1)):ncol(Info)] = 
          Info[(1+ncol(BM1)+ncol(TM)):ncol(Info),(1+ncol(TM)+ncol(BM1)):ncol(Info),drop=FALSE]/(weighting^2)
        globrelD = 0
        relD = 1
        globTM = TM
        globTF = TF
        V=chol2inv(chol(Info))
        for (r in 1:searches) {
          dmax=DMax(TF,TM,BM,V,IBF[,u,drop=TRUE])
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
              rownames(TF)=NULL
              unavailable=suppressMessages(as.integer(rownames(match_df(TF,TF[s1,,drop=FALSE]))))
              z= seq_len(nrow(TF))[IBF[,u]==IBF[s1,u] & BF[,u]!=BF[s1,u] & !(seq_len(nrow(TF))%in%unavailable)]
              if (length(z)==0 & counter<501) next 
              else if (length(z)==0 & counter>500) break
              if (length(z)>1) s=c(s1,sample(z,1)) else s=c(s1,z)
              testDswap=dMat(TM,BM,V,s)[2,1,drop=FALSE]
              if (testDswap<tol & counter<501) next
              else if (testDswap<tol & counter>500) break
              Dswap=testDswap
              break
            }
            if (counter>500) break
            relD=relD*Dswap 
            V=UpDate(V,(TM[s[1],,drop=TRUE]-TM[s[2],,drop=TRUE]),(BM[s[2],,drop=TRUE]-BM[s[1],,drop=TRUE]),ncol(TM))
            TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),,drop=FALSE] 
            TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),,drop=FALSE]
          } #jumps
        } #searches
        TF=globTF 
        TM=globTM
        E1=eigen(diag(ncol(TM))-tcrossprod(crossprod(TM,BM1)),symmetric=TRUE,only.values = TRUE)
        if (all(E1$values>tol)) { 
          D_Effic[u]=prod(E1$values)**(1/ncol(TM))
          A_Effic[u]=ncol(TM)/sum(1/E1$values)
        }
        if (!is.null(BM2)) {
          E2=eigen(diag(ncol(TM))-tcrossprod(crossprod(TM,BM2)),symmetric=TRUE,only.values = TRUE)
          if (all(E2$values>tol)) { 
            D_IntEff[u]=prod(E2$values)**(1/ncol(TM))
            A_IntEff[u]=ncol(TM)/sum(1/E2$values)
          }
        } else {
          D_IntEff[u]=0
          A_IntEff[u]=0
        }
      } # list length
      Effics1=data.frame(First_order_model=addfactors,effects=Add_levs,"D-Efficiency"= round(D_Effic,5),"A-Efficiency"= round(A_Effic,5))
      Effics2=data.frame(Second_order_model=prodfactors,effects=Int_levs,"D-Efficiency"= round(D_IntEff,5),"A-Efficiency"= round(A_IntEff,5))
      Effics=cbind(Effics1,Effics2)
      list(TF=TF,TM=TM,Effics=Effics)
    } 
 
   # *******************************************************************************************************************
   # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
   # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
   # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
   # *******************************************************************************************************************
   UpDate=function(V,t,b,v) {
     VTTt=crossprod(V[1:v,1:v,drop=FALSE],t)
     VBBb=crossprod(V[(1+v):ncol(V),(1+v):ncol(V),drop=FALSE],b)
     VTBt=crossprod(V[1:v,(1+v):ncol(V),drop=FALSE],t)
     VTBb=crossprod(V[(1+v):ncol(V),1:v,drop=FALSE],b)
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
     div=(1+a)*(1-b)+c*c
     V[1:v,1:v] = V[1:v,1:v, drop=FALSE]-(tcrossprod(f1)*(1-b) - tcrossprod(g1)*(1+a) + (tcrossprod(g1,f1)+tcrossprod(f1,g1))*c)/div
     V[(1+v):ncol(V),(1+v):ncol(V)] = V[(1+v):ncol(V),(1+v):ncol(V),drop=FALSE]- 
       (tcrossprod(f2)*(1-b) - tcrossprod(g2)*(1+a) + (tcrossprod(g2,f2)+tcrossprod(f2,g2))*c)/div
     V[1:v,(1+v):ncol(V)] = V[1:v,(1+v):ncol(V),drop=FALSE]-
       (tcrossprod(f1,f2)*(1-b) - tcrossprod(g1,g2)*(1+a) + (tcrossprod(g1,f2)+tcrossprod(f1,g2))*c)/div
     V[(1+v):ncol(V), 1:v]= t(V[1:v,(1+v):ncol(V),drop=FALSE])
     return(V)
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
  # non singular starting design for the first nunits of TF given the design matrix TM
  # ******************************************************************************************************************
  nonSingular=function(TF,mod,fixed) {
    TM=model.matrix(as.formula(treatments_model[mod]),TF)
    rank=qr(TM[1:nrow(BF),])$rank
      for (z in seq_len(1000)) {
        if (identical(ncol(TM),rank)) break
        k=sample(nlevels(fixed),1)
        available=which(fixed==levels(fixed)[k])
        swapout=available[available<=nrow(BF)]
        swapin= available[available>nrow(BF)]
        if (length(swapin)==0 | length(swapout)==0) next 
        if (length(swapin)>1)  p2 = sample(swapin, 1) else  p2 = swapin
        if (length(swapout)>1) p1 = sample(swapout, 1) else p1 = swapout 
        TM[c(p1,p2),] = TM[c(p2,p1),]
        TF[c(p1,p2),] = TF[c(p2,p1),]
        newrank = qr(TM[1:nrow(BF),])$rank
        if (newrank<rank) {
          TM[c(p1,p2),] = TM[c(p2,p1),]
          TF[c(p1,p2),] = TF[c(p2,p1),]
        } else rank=newrank
      }
      if (!identical(ncol(TM),rank)) stop ("Cannot find a starting design")
      return(TF) 
  }
  # *********************************************************************************************************************
  # Fractional factorials for quantitative level factors - allows replacement of swapped rows from candidate set
  # *******************************************************************************************************************
  factorial = function(TF,treatments_model) {
    fullrand = unlist(lapply(1:ceiling(nrow(BF)/nrow(TF)),function(j){sample(seq_len(nrow(TF)))}))
    TF = TF[fullrand,,drop=FALSE]
    for (mod in 1:length(treatments_model)) {
      TM=model.matrix(as.formula(treatments_model[mod]),TF)
      if (ncol(TM)>nrow(BF)) stop("too many treatment parameters to be estimated for the given number of experimental units")
      Dmax = exp((determinant(crossprod(TM),logarithm = TRUE)$modulus)/ncol(TM))/nrow(TM)
      if (mod==1) fixed = factor(rep(1,nrow(TF)))
      if (mod>1) {
        included=sapply(1:ncol(TF),function(j){grepl(colnames(TF)[j],treatments_model[mod-1],fixed=TRUE)}) 
        fixed = interaction( data.frame(lapply(TF[,included,drop=FALSE], factor)), drop=TRUE) }
      if (nlevels(fixed)>nrow(TF)) stop("treatments_model has too many restrictions - is the model fully restricted by previous terms?")
      gDfrac = 0
      for (i in 1:searches) {
        for (x in 1:nlevels(fixed)) 
          TF[seq(nrow(TF))[fixed==levels(fixed)[x]],] = TF[sample(seq(nrow(TF))[fixed==levels(fixed)[x]]),,drop=FALSE]
        TF= nonSingular(TF,mod,fixed) 
        TM=model.matrix(as.formula(treatments_model[mod]),TF)
        V = chol2inv(chol(crossprod(TM[1:nrow(BF),,drop=FALSE]))) # V is the variance matrix of the the non-singular starting design
        counter=0
        locrelD=1
        repeat {
          kmax=1
          counter=counter+1
          for (k in 1: nlevels(fixed)) {
            available = which(fixed==levels(fixed)[k])
            swapout = available[available<=nrow(BF)]
            swapin = available[available> nrow(BF)]
            if ( (length(swapin)==0 | length(swapout)==0) & k<nlevels(fixed)  ) next 
            if ( (length(swapin)==0 | length(swapout)==0) & k==nlevels(fixed) ) break 
            if (length(swapin)>1)  swap_in = sample(swapin, min(1024,length(swapin)))  else swap_in = swapin
            if (length(swapout)>1) swap_out=sample(swapout, min(1024,length(swapout))) else swap_out = swapout
            M1VM2 =        tcrossprod(tcrossprod(TM[swap_out,,drop=FALSE],V),TM[swap_in,,drop=FALSE])
            M1VM1 = 1-diag(tcrossprod(tcrossprod(TM[swap_out,,drop=FALSE],V),TM[swap_out,,drop=FALSE]))
            M2VM2 = 1+diag(tcrossprod(tcrossprod(TM[swap_in,,drop=FALSE],V),TM[swap_in,,drop=FALSE]))
            Z = M1VM2**2 + tcrossprod(M1VM1,M2VM2)
            z = which(Z == max(Z), arr.ind = TRUE)[1,]
            if (Z[z[1],z[2]] >= kmax) {
              kmax=Z[z[1],z[2]]
              pi=swap_out[z[1]]
              pj=swap_in[z[2]]
            } 
          }
          if (kmax>(1+tol)) {
            counter=0
            locrelD=locrelD*kmax
            V = fractUpDate(V,TM[pi,],TM[pj,])   # parameters(V,row_swappedout,row_swappedin) 
            TM[c(pi,pj),]=TM[c(pj,pi),]
            TF[c(pi,pj),]=TF[c(pj,pi),]
          }
          if (counter==5)  break
          if (counter==1 & length(swapin)<=1024 & length(swapout)<=1024) break
        }
        Dfrac = exp((determinant(crossprod(TM[1:nrow(BF),,drop=FALSE]),logarithm = TRUE)$modulus)/ncol(TM))/nrow(BF)
        if (Dfrac>gDfrac) {
          gTF=TF
          gDfrac=Dfrac
        }
        if (isTRUE(is.factor(TF[,mod])) & isTRUE(all.equal(gDfrac,Dmax,tolerance = tol))) break
      }
    }
    return(gTF[1:nrow(BF),,drop=FALSE])
  }
  # *********************************************************************************************************************
  #  Tests if all factors are fully nested AND if all nested blocks are of equal size for each level of nesting
  # *******************************************************************************************************************
  equinested=function(BF) {
    for (i in 1:(ncol(BF))) { 
      if ( !isTRUE(nlevels(interaction( BF[,c(1:i)],drop=TRUE))==nlevels(BF[,i,drop=TRUE]))) return(FALSE)
      if ( length(unique(table(BF[,i]))) > 1) return(FALSE)
    }
    return(TRUE)
  }

  # *********************************************************************************************************************
  #  treatment information
  # *******************************************************************************************************************
  information=function(TF,treatments_model) {
  info=sapply(1:length(treatments_model),function(i) {
    TM = model.matrix(as.formula(treatments_model[i]),TF)
    exp(determinant(crossprod(TM),logarithm = TRUE)$modulus/ncol(TM))/nrow(TM)})
  return(info)
  }
  # *************************************************************************************************************************
  # Main design program tests, optimizes treatment design then optimizes block design sequentially for each added block factor
  # *************************************************************************************************************************
  TF=data.frame(treatments)
  if (is.null(blocks)) blocks= gl(1,nrow(TF))
  BF=data.frame(blocks)
  if (!all(sapply(BF,is.factor))) stop("blocks must be factors")
  if (jumps<1 | jumps%%1!=0 | jumps>10) stop("min. jumps is 1 and max. jumps is 10")
  if (is.null(searches)) searches=1+5000%/%nrow(BF)
  if (searches<1 | searches%%1!=0) stop(" searches parameter is invalid")
  if (is.null(treatments_model)) treatments_model=paste("~",paste(colnames(TF),collapse="*"))
  TF=do.call(rbind,lapply(1:ceiling(nrow(BF)/nrow(TF)),function(j){TF})) # ensure TF is at least as long as BF
  infobase=information(TF,treatments_model)
  TF=factorial(TF,treatments_model)
  addfactors  =unlist(lapply(1:ncol(BF),function(j){ paste0("(",paste0(colnames(BF)[1:j],collapse="+"),")"  )  }) )
  prodfactors =unlist(lapply(1:ncol(BF),function(j){ paste0("(",paste0(colnames(BF)[1:j],collapse="+"),")^2")  }) )
   # nested equi-block designs with a single unstructured treatment factor optimized by using blocks() function
  if (equinested(BF) & ncol(TF)==1 & is.factor(TF[,1])) {
    blevs=c(1,sapply(BF, nlevels))
    blevs=sapply(2:length(blevs), function(i){blevs[i]/blevs[i-1]})
    Z=blocks(rep(1,nlevels(TF[,1])),table(TF[,1,drop=TRUE]),blevs,searches=searches,seed=seed,jumps=1)
    Treatments=Z$Treatments
    Design=Z$Design
    levels(Treatments[,1]) = levels(TF[,1])
    levels(Design[,ncol(Design)]) = levels(TF[,1])
    blocksModel=Z$Blocks_model[,-1]
    blocksModel=data.frame(Model=addfactors,blocksModel)
    weighting=NULL
  } else {
    if (!is.numeric(weighting)) stop("weighting must be a number between 0 and 1")
    if (weighting<0 | weighting>1) stop("weighting must be a number between 0 and 1")
    if (max(sapply(BF, nlevels))>1) {
      TM = scale(model.matrix(as.formula(treatments_model[length(treatments_model)]),TF),center=TRUE,scale=FALSE)[,-1,drop=FALSE] # centred model.matrix
      TM = qr.Q(qr(TM)) # orthogonal basis for TM 
      Opt=blocksOpt(TF,TM,BF)
      TM=Opt$TM
      TF=Opt$TF
      blocksModel=Opt$Effics
    } else blocksModel=NULL
    Design=cbind(BF,TF)
    Treatments=count(TF)
  }
  rownames(Design)=NULL
  rownames(Treatments)=NULL
  
  # Table for Treatments_model DF and efficiency
  EfficsTrts=round(information(TF,treatments_model)/infobase,6)
  DF=sapply(1:length(treatments_model),function(i) {ncol(model.matrix(as.formula(treatments_model[i]),TF)) -1})
  treatsModel=data.frame(cbind("Treatment model" = treatments_model,"Model DF" = DF ,"D-Efficiency" = EfficsTrts))
  
  list(Treatments=Treatments,Design=Design,Treatments_model=treatsModel,Blocks_model=blocksModel,
       weighting=weighting,seed=seed,searches=searches,jumps=jumps)
  }
  