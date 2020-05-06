#' @title General block and treatment designs.
#'
#' @description
#' Constructs D-optimal block and treatment designs for any feasible combination of block 
#' factors and any feasible linear treatment model.   
#' 
#' @param treatments a single treatment factor or a data frame containing any combination of 
#' treatment factors with any combination of treatment vectors.
#' 
#' @param blocks a single block factor or a data frame containing one or more qualitative level block
#' factors in the required order of fitting.
#' 
#' @param treatments_model a treatment model formula for the required treatment design.
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
#' @details
#' 
#' \code{treatments} is a factor or a data frame containing one or more qualitative
#' or quantitative level treatment vectors defining a set of candidate treatments.
#' If the size of the candidate set is different from the size of the block design, 
#' or if the \code{treatments_model} is structured, the treatment design is optimized by
#' selection with replacement, otherwise the full fixed candidate set is used for the treatment design. 
#' 
#' \code{blocks} is a factor or a data frame containing one or more qualitative level block factors. The design is optimized
#' for each block factor added sequentially from left to right with the optimization constrained to be constant for all
#' previously added block factors. The length of the \code{blocks} object defines the total number of plots in the design. 
#'  
#' \code{treatments_model} is either a single formula or a compound formula split by the \code{|} operator.
#' The left hand side of each \code{|}, assuming all remaining \code{|} are replaced by 
#' \code{+}, is a partial model formula. Partial model formula define partial design matrices which are fitted
#' and optimized sequentially from left to right. Sequential model fitting provides improved flexibility
#' for fitting factors or variables of different status or importance (see examples below).
#' 
#' The treatment design criterion for each partial model is the generalized variance of the treatment design for that
#' model (D-optimality) and the design efficiency is the ratio of the generalized variance of the full candidate treatment model
#' for the required treatment model relative to the generalized variance of the optimized design. The efficiency is necessarily
#' less than or equal to 1 for purely factorial models but may exceed 1 for polynomial models. 
#' 
#' For crossed blocks, a weighting factor w is used to determine the relative importance
#' of the block main effects versus the block 2-factor interaction effects. If w = 0 the algorithm fits a simple additive
#' main effects model whereas if w = 1 the algorithm fits both main effects and 2-factor interaction effects.
#' For intermediate 0 < w < 1, the 2-factor interaction effects are downweighted relative to the main effects where
#' the smaller the value of w the greater the downweighting. The default weighting is 0.5 and 
#' provided that all block effects are estimable, this weighting gives a compromise design where all 
#' main effects and 2-factor interaction effects are included but where the main effects are given
#' greater realative importance than the 2-factor interaction effects.
#' See \code{vignette(package = "blocksdesign")} for more details. 
#' 
#' For more details see \code{vignette(package = "blocksdesign")}  
#' 
#' @return
#' \item{Treatments}{The treatments included in the design and the replication of each individual 
#' treatment taken in de-randomized standard order.}
#' \item{Design}{The design layout showing the randomized allocation of treatments to blocks and plots.}
#' \item{Treatments_model}{The fitted treatment model, the number of model parameters (DF)
#'   and the D-efficiency of each sequentially fitted treatment model}
#' \item{Blocks_model}{The blocks sub-model design and 
#' the D- and A-efficiency factors of each successively fitted sub-blocks model.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#' 
#' @references
#' 
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons.
#' 
#' @examples
#' 
#' ## For optimum results, the number of searches may need to be increased.
#' 
#' ## 4 replicates of 12 treatments with 16 nested blocks of size 3
#' # rectangular lattice see Plan 10.10 Cochran and Cox 1957.
#' blocks = data.frame(Main = gl(4,12), Sub = gl(16,3))
#' design(treatments = factor(rep(1:12,4)),blocks)
#' 
#'  ## 6 replicates of 5 treatments with 2 super blocks of size 15 and 10 nested blocks of size 3
#'  blocks=data.frame(Main = gl(2,15,30), Sub = gl(10,3,30))
#'  design(treatments=gl(5,1,30),blocks)
#' 
#' ## 4 x 12 design for 4 replicates of 12 treatments with 3 plots in each intersection block
#' ## The optimal design is Trojan with known A-efficiency = 22/31 for the intersection blocks
#' blocks = data.frame(Rows = gl(4,12), Cols = gl(4,3,48))
#' design(treatments = factor(1:12),blocks)
#' 
#' ## 4 x 12 design for 4 replicates of 12 treatments with 3 sub-column blocks nested 
#' ## as above but showing 3 sub-columns nested within each main column
#' blocks = data.frame(Rows = gl(4,12), Cols = gl(4,3,48), subCols = gl(12,1,48))
#' \donttest{design(treatments = factor(1:12),blocks,searches=200)}
#' 
#' ## 4 x 13 Row-and-column design for 4 replicates of 13 treatments 
#' ## Youden design Plan 13.5 Cochran and Cox (1957).
#' blocks = data.frame(Rows = gl(4,13), Cols = gl(13,1,52))
#' \donttest{design(treatments = factor(1:13),blocks,searches = 700)}
#' 
#' ## differential replication 
#' treatments=factor(c(rep(1:12,2),rep(13,12)))
#' blocks = data.frame(Main = gl(2,18),  Sub = gl(12,3,36))
#' design(treatments,blocks,searches = 5)
#' 
#' ## 48 treatments in 2 replicate blocks with 2 nested rows in each replicate and 3 main columns
#' ##  (Reps/Rows) x Cols
#' blocks = data.frame(Reps = gl(2,48), Rows = gl(4,24,96), Cols = gl(3,8,96))
#' design(treatments=factor(1:48),blocks,searches=5)
#' 
#' ## 48 treatments in 2 replicate blocks with 2 main columns
#' ## The default weighting gives non-estimable Reps:Cols effects due to inherent aliasing
#' ## Increased weighting gives estimable Reps:Cols effects but non-orthogonal main effects
#' blocks = data.frame(Reps = gl(2,48), Cols = gl(2,24,96))
#' design(treatments=factor(1:48),blocks,searches=5)
#' design(treatments=factor(1:48),blocks,searches=5,weighting=.9)
#' 
#' ## Factorial treatment designs defined by a single factorial treatment model
#' 
#' ## Main effects of five 2-level factors in a half-fraction in 2/2/2 nested blocks design 
#' ## (may require 100's of repeats to find a fully orthogonal solution - a VERY long wait!)
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:2),
#'  F3 = factor(1:2), F4 = factor(1:2), F5 = factor(1:2))
#' blocks = data.frame(b1 = gl(2,8),b2 = gl(4,4),b3 = gl(8,2))
#' model=" ~ F1 + F2 + F3 + F4 + F5"
#' \donttest{repeat {z = design(treatments,blocks,treatments_model=model,searches=50)
#' if ( isTRUE(all.equal(z$Blocks_model[3,3],1) ) ) break }
#'  print(z)}
#'  
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:2), F3 = factor(1:2), 
#' F4 = factor(1:2), F5 = factor(1:2))
#' blocks = data.frame(blocks = gl(4,8))
#' model = " ~ (F1 + F2 + F3 + F4 + F5)^2"
#' design(treatments,blocks,treatments_model=model,searches = 10)
#' 
#' # Main effects of five 2-level factors in a half-fraction of 
#' # a 4 x 4 row-and column design.
#' treatments = expand.grid(F1 = factor(1:2), F2 = factor(1:2), F3 = factor(1:2), 
#' F4 = factor(1:2), F5 = factor(1:2))
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
#' treatments = expand.grid(F1 = factor(1:3), F2 = factor(1:3), F3 = factor(1:3), 
#' F4 = factor(1:3))
#' blocks = data.frame(main = gl(3,9))
#' model = " ~ F1 + F2 + F3 + F4"
#'\donttest{ design(treatments,blocks,treatments_model=model,searches=25)}
#' 
#' # 2nd-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' # (may require many repeats to find a fully orthogonal solution)
#' treatments = expand.grid(F1 = factor(1:3), F2 = factor(1:3), F3 = factor(1:3), 
#' F4 = factor(1:3), F5 = factor(1:3))
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
#'\donttest{ design(treatments,blocks,treatments_model=model,searches=5)}
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
#' ## Factorial treatment designs defined by sequentially fitted factorial treatment models
#' 
#' ## 2 varieties x 3 levels of N x 3 levels of K assuming 1st-order interactions and 12 plots
#' ## the single stage model gives an unequal 7 + 5 split for the two varieties
#' ## whereas the two stage model forces an equal 6 + 6 split
#' ## NB the two stage model is slightly less efficient than the single stage model 
#' treatments = expand.grid(Variety = factor(rep(1:2)), N = 1:3, K = 1:3)
#' blocks=data.frame(main=gl(1,12))
#' treatments_model = " ~  (Variety + N + K)^2  + I(N^2) + I(K^2)"
#' \donttest{design(treatments,blocks,treatments_model=treatments_model,searches=10)} 
#' treatments_model = " ~ Variety | (Variety + N + K)^2 + I(N^2) + I(K^2)"
#' \donttest{design(treatments,blocks,treatments_model=treatments_model,searches=10)}
#' 
#' ## A 6 x 6 row-and-column design with linear row by linear column interaction.
#' ## Crossed blocks with interactions fitted in the treatments model and additive 
#' ## treatments fitted inthe blocks model as a dual design
#' ##  see vignette(package = "blocksdesign") for further discussion
#' ## may require many separate attempts to get the best overall design efficiency 
#' LS_grid   = expand.grid(rows=factor(1:6), cols=factor(1:6))
#' blocks = data.frame(varieties=factor(rep(1:6,6)))
#' lin_rows = as.numeric(levels(LS_grid$rows))[LS_grid$rows]
#' lin_cols = as.numeric(levels(LS_grid$cols))[LS_grid$cols]
#' latin_sq = "~  rows | cols + lin_rows:lin_cols "
#' \donttest{design(LS_grid,blocks,latin_sq,searches=2000)} 
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames 
#' @importFrom lme4 lmer
#' @importFrom plyr count
#' @importFrom plyr match_df
#' 
  design = function(treatments,blocks=NULL,treatments_model=NULL,weighting=0.5,searches=NULL,seed=NULL,jumps=1) {
    tol = .Machine$double.eps ^ 0.5
    options(contrasts=c('contr.SAS','contr.poly'))
    if (!is.null(seed)) set.seed(seed) 
    
    # *********************************************************************************************************************
    #  Full rank
    # *********************************************************************************************************************
    fullrankModel=function(TF,model_formula) {
      TM=model.matrix(as.formula(model_formula),TF)
      Q = qr(TM)
      TM = TM[,Q$pivot[1:Q$rank],drop=FALSE] # full rank
      return(TM)
    }
    
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
    DMax=function(TF,TM,BM,V,mainBlocks) {
      locrelD=1
      mainSets=tabulate(mainBlocks)
      nSamp=pmin(rep(8,nlevels(mainBlocks)), mainSets)
      repeat {
        kmax=1
        for (k in 1: nlevels(mainBlocks)) {
          s=sort(sample((1:nrow(TF))[mainBlocks==levels(mainBlocks)[k]] , nSamp[k])) 
          dMat=dMat(TM,BM,V,s)
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
    # ********************************************************************************************************************
    # Random swaps
    # ********************************************************************************************************************    
    Swaps=function(TF,BF,pivot,rank,restrict,blocks) {
      candidates=NULL
      counter=0
      while (length(candidates)==0 & counter<500) {
        counter=counter+1
        swapOut=pivot[1:rank]
        swapIn=pivot[(1+rank):nrow(TF)]
        if (length(swapIn)>1) s1=sample(swapIn,1) else s1=swapIn
        rownames(TF)=NULL
        sameTrts=suppressMessages(as.integer(rownames(match_df(TF, TF[s1,,drop=FALSE]))))
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
        TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
        Q=qr(t(cbind(TM,BM)))
        if (Q$rank>rank) {
          rank=Q$rank
          pivot=Q$pivot
          TF[c(swap[1],swap[2]),]=TF[c(swap[2],swap[1]),]
        } else 
          TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
      }
      return(list(TF=TF,TM=TM,fullrank=FALSE))
    }
    # ************************************************************************************************************************
    # Optimize the  blocks assuming a possible set of Main block constraints Initial completely randomized starting design.
    # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
    # Each factor in BF may include levels from preceding factors hence qr pivoting is essential to sweep out factor dependencies  
    # ************************************************************************************************************************
    blocksOpt=function(TF,TM,BF,searches) {
      IBF=data.frame(factor(rep(1,nrow(TF))),lapply(1:ncol(BF), function(i) droplevels(interaction(BF[,c(1:i)]))))
      colnames(IBF)=c("mean",unlist(lapply(1:ncol(BF),function(j){paste0(colnames(BF)[1:j],collapse=".")})))
      D_Effic=rep(0,length=ncol(BF))
      A_Effic=rep(0,length=ncol(BF))
      D_IntEff=rep(0,length=ncol(BF))
      A_IntEff=rep(0,length=ncol(BF))
      Int_levs=rep(0,length=ncol(BF))
      Add_levs=rep(0,length=ncol(BF))
      
      orthogM=function(M) {
        QR = qr(M) # qr transformation
        if (QR$rank < min(nrow(M),ncol(M)))
          QR = qr(M[, QR$pivot[1:QR$rank] ,drop=FALSE]) # removes any aliased block effects then finds qr transformation
        M = qr.Q(QR) # orthogonal basis for M where Q'Q=I
        M
      }
      
      for (u in 1:ncol(BF)) {
        BM1 = scale(model.matrix(as.formula(paste("~",paste0(colnames(BF)[1:u],collapse="+"))),BF), center = TRUE, scale = FALSE)[,-1,drop=FALSE]
        BM1=orthogM(BM1)
        Add_levs[u]=ncol(BM1)
        
        BM2=scale(model.matrix(as.formula(paste("~",paste0("(",paste0(colnames(BF)[1:u],collapse="+"),")^2"))),BF),center=TRUE,scale = FALSE)[,-1,drop=FALSE]
        BM2=orthogM(BM2)
        
        Int_levs[u]=ncol(BM2)
        if ((ncol(TM) + ncol(BM2)+1) > nrow(TM)) BM2=NULL
        
        if (!is.null(BM2)) {
          NS = blocksNonSingular(TF,TM,BM2,IBF[,u],BF[,u])
          if (NS$fullrank==FALSE) BM2=NULL
        }
        
        if (is.null(BM2))
          NS = blocksNonSingular(TF,TM,BM1,IBF[,u],BF[,u])
        
        if (NS$fullrank==FALSE)  stop("Cannot find a non-singular starting block design - perhaps try a simpler design? ")
        TM = NS$TM
        TF = NS$TF 
        
        if ( !is.null(BM2) & weighting>0)BM=BM2 else BM=BM1
        
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
          
          dmax=DMax(TF,TM,BM,V,IBF[,u])
          
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
              testDswap=dMat(TM,BM,V,s)[2,1]
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
  nonSingular=function(TF,TM,nunits,fixed) {
    row.names(TF)=NULL
    singular=FALSE
    for (t in 1:searches) {
      rank=qr(TM[1:nunits,])$rank
      if (identical(ncol(TM),rank)) break
      if (nunits==nrow(TF)) stop("Unable to find non-singular starting design by random search")
      for (z in seq_len(10000)) {
        if (nlevels(fixed)>1) k=sample(nlevels(fixed),1) else k=1
        available=which(fixed==levels(fixed)[k])
        swapout=available[available<(nunits+1)]
        swapin=available[available>nunits]
        if (length(swapin)==0 | length(swapout)==0) next 
        if (length(swapin)>1) p2 = sample(swapin, 1) 
        else p2 = swapin
        if (length(swapout)>1) p1 = sample(swapout, 1)
        else p1 = swapout 
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
      if (t>=searches) TF=NULL
    }
    return(list(TF=TF,TM=TM))
  }
  # *********************************************************************************************************************
  # Fractional factorials for quantitative level factors - allows replacement of swapped rows from candidate set
  # *******************************************************************************************************************
  factorial = function(TF,nunits,treatments_model) {
    fullrand = unlist(lapply(1:ceiling(nunits/nrow(TF)),function(j){sample(seq_len(nrow(TF)))}))
    TF = TF[fullrand,,drop=FALSE]
    if (is.null(treatments_model) & isTRUE(all.equal(nunits,nrow(TF)))) {
      TM = fullrankModel(TF,paste("~",paste(colnames(TF),collapse="*")))
      return(list(TF = TF, TM = TM, Eff = 1, DF = ncol(TM)))
    } else if (is.null(treatments_model)) 
      treatments_model = paste("~",paste(colnames(TF),collapse="*"))
    
    Effics=rep(0,length(treatments_model))
    DF=rep(0,length(treatments_model))
    
    for (mod in 1:length(treatments_model)) {
      if (mod>1) { 
        included=sapply(1:ncol(TF),function(j){grepl(colnames(TF)[j],treatments_model[mod-1],fixed=TRUE)}) 
        fixed = interaction( data.frame(lapply(TF[,included,drop=FALSE], factor))  ,drop=TRUE)
      } else fixed = factor(rep(1,nrow(TF))) 
      gDfrac = 0
      
      if (nlevels(fixed)>=nrow(TF)) stop("treatments_model has too many restrictions - is a model term fully restricted by previous terms?")
      for (i in 1:searches) {
        for (x in 1:nlevels(fixed)) 
          TF[seq(nrow(TF))[fixed==levels(fixed)[x]],]=TF[sample(seq(nrow(TF))[fixed==levels(fixed)[x]]),,drop=FALSE]
        TM = fullrankModel(TF,treatments_model[mod])
        if (ncol(TM)>nunits) stop("too many treatment parameters to be estimated for the given number of experimental units")
        W = nonSingular(TF,TM,nunits,fixed) 
        TF = W$TF
        TM = W$TM
        if (is.null(TF)) stop("Unable to find non-singular starting design for factorial treatment design by random search")
        if (i==1) Dmax = exp((determinant(crossprod(TM),logarithm = TRUE)$modulus)/ncol(TM))/nrow(TM)
        V = chol2inv(chol(crossprod(TM[1:nunits,,drop=FALSE]))) # V is the variance matrix of the the non-singular starting design
        counter=0
        locrelD=1
        repeat {
          kmax=1
          counter=counter+1
          for (k in 1: nlevels(fixed)) {
            available=which(fixed==levels(fixed)[k])
            swapout=available[available<(nunits+1)]
            swapin =available[available>nunits]
            if ( (length(swapin)==0 | length(swapout)==0)) 
             if (k<nlevels(fixed)) next else break
            
            if (length(swapin)>1) swap_in = sample(swapin, min(1024,length(swapin)) ) 
            else swap_in = swapin
            if (length(swapout)>1) swap_out = sample(swapout, min(1024,length(swapout)) ) 
            else swap_out = swapout
            
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
          if (counter==8 | (counter==1 & length(swapin)<=1024 & length(swapout)<=1024))  break
        }
        Dfrac = exp((determinant(crossprod(TM[1:nunits,,drop=FALSE]),logarithm = TRUE)$modulus)/ncol(TM))/nunits
        if (Dfrac>gDfrac) {
          gTM=TM
          gTF=TF
          gDfrac=Dfrac
        }
        if (isTRUE(is.factor(TF[,mod])) & isTRUE(all.equal(gDfrac,Dmax,tolerance = tol))) break
      }
      Effics[mod]=gDfrac/Dmax
      DF[mod]=ncol(gTM)-1
    }
   list(TF = gTF[1:nunits,,drop=FALSE], TM = gTM[1:nunits,,drop=FALSE], Eff = round(Effics,7),DF = DF)
  }

  # *************************************************************************************************************************
  # Main design program tests, optimizes treatment design then optimizes block design sequentially for each added block factor
  # *************************************************************************************************************************
  TF=data.frame(treatments)
  if (is.null(blocks)) blocks= gl(1,nrow(TF))
  BF=data.frame(blocks)
  if (!all(sapply(BF,is.factor))) stop("blocks must be factors")
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) stop(" max. jumps is 10")
  nunits=nrow(BF)
  if (is.null(searches)) 
    if (nunits<100) 
      searches=25 else if (nunits<1000) 
        searches=2500%/%nunits else if (nunits<5000) 
          searches=2 else 
            searches=1
  if (!is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")

  null_model=is.null(treatments_model)
  if (null_model) treatments_model=paste("~",paste(colnames(TF),collapse="*"))
  model=strsplit(treatments_model,split='|', fixed=TRUE)
  treatments_model = sapply( 1:length(model[[1]]),function(j){ paste0( model[[1]][1:j] ,collapse="+")  })
  Z=factorial(TF,nunits,treatments_model)
  TF=data.frame(Z$TF)
  TM=scale(Z$TM, center = TRUE, scale = FALSE)[,-1,drop=FALSE]
  TM = qr.Q(qr(TM)) # orthogonal basis for TM 
  prodfactors=unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse="*")}) )
  addfactors =unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse="+")}) )
  treatsModel=data.frame(cbind("Treatment model" = treatments_model,"Model DF" = Z$DF ,"D-Efficiency" = Z$Eff))
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  equinested =  all(sapply( 1:ncol(BF),function(i) {
    T  = table( interaction( BF[,c(1:i)],drop=TRUE) , BF[,i]  )
    all(nrow(T)==ncol(T) & T[!diag(nrow(T))] == 0 & abs( max(diag(T)) - min(diag(T)) ) < tol ) 
  }))
  
 # nested equi-block designs with a single unstructured treatment factor uses blocks() function
  if (equinested & ncol(TF)==1 & is.factor(TF[,1])) {
    tlevs=levels(TF[,1])
    blevs=c(1,sapply(BF, nlevels))
    blevs=sapply(2:length(blevs), function(i){blevs[i]/blevs[i-1]})
    Z=blocks(rep(1,nlevels(TF[,1])),table(TF[,1]),blevs,searches=searches,seed=seed,jumps=1)
    Treatments=Z$Treatments
    levels(Treatments[,1]) = tlevs
    Design=Z$Design
    levels(Design[,ncol(Design)]) = tlevs
    blocksModel=Z$Blocks_model[,-1]
    blocksModel=data.frame(Model=addfactors,blocksModel)
    weighting=NULL
  } else {
    if (!is.numeric(weighting)) stop("weighting must be a number between 0 and 1")
    if (weighting<0 | weighting>1) stop("weighting must be a number between 0 and 1")
    if (max(sapply(BF, nlevels))>1) {
      Opt=blocksOpt(TF,TM,BF,searches)
      TM=Opt$TM
      TF=Opt$TF
      blocksModel=Opt$Effics
    } else blocksModel=NULL
    Design=cbind(BF,TF)
    Treatments=count(TF)
  }
  rownames(Design)=NULL
  rownames(Treatments)=NULL
  list(Treatments=Treatments,Design=Design,Treatments_model=treatsModel,Blocks_model=blocksModel,
       weighting=weighting,seed=seed,searches=searches,jumps=jumps)
  }
  