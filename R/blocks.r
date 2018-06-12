#' @title Block designs for unstructured treatment sets
#'
#' @description
#'
#' Constructs randomized nested block designs for unstructured treatment sets with any
#' feasible depth of nesting.
#'
#' @details
#'
#' Constructs randomized nested block designs with arbitrary depth of nesting for arbitrary unstructured treatment sets.
#' 
#' The \code{treatments} parameter is a set of numbers that partitions the total number of treatments into equally 
#' replicated treatment sets while the \code{replicates} parameter is a matching set of numbers that defines the replication of each equally 
#' replicated treatment set.
#' 
#' The \code{blocks} parameter, if any, defines the number of blocks for each level of nesting from the highest to the lowest. The
#' first number, if any, is the number of nested row blocks in the first-level of nesting, the second number, if any, is the number of nested row blocks in
#' the second-level of nesting and so on down to any required feasible depth of nesting.
#'
#' Block sizes are as nearly equal as possible and will never differ by more than a single plot for any particular block classification. 
#'
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number of single unreplicated treatments 
#' to give augmented blocks that never differ in size by more than a single plot. General crossed block designs are more complex and currently 
#'
#' Square lattice designs are resolvable incomplete block designs for r replicates of p*p treatments arranged in blocks of size p where
#' r < p+2 for prime or prime power p or r < 4 for general p. Square lattice designs are constructed algebraically from Latin squares or MOLS.
#'
#' Lattice designs based on prime-power MOLS require the \code{\link[crossdes]{MOLS}} package.
#'
#' All other designs are constructed numerically by optimizing a D-optimality criterion.
#'
#' Outputs:
#'
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
#'  \item  A table showing the replication number of each treatment in the design. \cr
#'  \item  A table showing the block levels and the achieved D-efficiency and A-efficiency factor for each nested level together
#'   with A-efficiency upper bounds, where available. \cr
#'  \item  A plan showing the allocation of treatments to blocks in the bottom level of the design.\cr
#' }
#'
#' @param treatments  a partition of the total required number of treatments into equally replicated treatment sets.
#'
#' @param replicates  a set of treatment replication numbers with one replication number for each partitioned treatment set.
#'
#' @param blocks the number of blocks nested in each preceding block for each level of nesting from the top-level block downwards. The top-level block is a
#' single super-block which need not be defined explicitly. 
#' 
#' @param seed  an integer initializing the random number generator. The default is a random seed.
#'
#' @param searches  the maximum number of local optima searched for a design optimization. The default number decreases
#' as the design size increases.
#'
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#'
#' @return
#' \item{Treatments}{A table showing the replication number of each treatment in the design.}
#' \item{Design}{Data frame giving the optimized block and treatment design in plot order.}
#' \item{Plan}{Data frame showing a plan view of the treatment design in the bottom level of the design.}
#' \item{BlocksEfficiency}{The D-efficiencies and the A-efficiencies of the blocks in each nested level of the design together with A-efficiency upper-bounds, where available.}
#' \item{seed}{Numerical seed used for random number generator.}
#' \item{searches}{Maximum number of searches used for each level.}
#' \item{jumps}{Number of random treatment swaps used to escape a local maxima.}
#'
#' @references
#'
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. https://CRAN.R-project.org/package=crossdes
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
#' ## Unstructured treatments partitioned into equally replicated treatment sets
#'
#' # 3 treatments x 2 replicates + 2 treatments x 4 replicates in two complete randomized blocks
#' blocks(treatments=c(3,2),replicates=c(2,4),searches=10)
#'
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block
#' blocks(treatments=50,replicates=4,blocks=c(4,5))
#'
#' # as above but with 20 additional single replicate treatments, one single treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),blocks=c(4,5))
#'
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,blocks=4)
#'
#' # 128 treatments x 2 replicates with two main blocks and 3 levels of nesting
#' blocks(128,2,c(2,2,2,2))
#' 
#' #' # 64 treatments x 4 replicates with 4 main blocks nested blocks of size 8 (lattice square)
#' blocks(64,4,c(4,8)) 
#' 
#' # 100 treatments x 4 replicates with 4 main blocks nested blocks of size 10 (lattice square)
#' blocks(100,4,c(4,10)) 
#' 
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames 
#' @importFrom crossdes MOLS
#'
 blocks = function(treatments,replicates,blocks=NULL,searches=NULL,seed=NULL,jumps=1) {
   options(contrasts=c('contr.SAS','contr.poly'))
   tol=.Machine$double.eps^0.5
   if (missing(treatments)|is.null(treatments)) stop(" Treatments missing or not defined ")
   if (missing(replicates)|is.null(replicates)) stop(" Replicates missing or not defined ")
   if (is.null(seed)) seed=sample(10000,1)
   if (is.na(seed) | !is.finite(seed) | is.nan(seed) | seed%%1!=0 | seed<0 ) stop(" seed parameter invalid  ")
   set.seed(seed)
   if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) stop(" number of jumps parameter is invalid (max is 10) ")
   if (any(replicates%%1!=0)|any(replicates<1)) stop(" replication numbers must be integers")
   if (anyNA(treatments)|any(is.nan(treatments))|any(!is.finite(treatments))|any(treatments%%1!=0)|any(treatments<1)) stop(" treatments parameter invalid")
   if (length(replicates)!=length(treatments)) stop("treatments and replicates parameters must both be the same length")

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
        mols=MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])
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
    if (!is.null(TF)) TF=factor(TF)
    return(TF)
  }
 
     # ******************************************************************************************************************************************************** 
     # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
     # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
     # ********************************************************************************************************************************************************
     DMax=function(VTT,VBB,VTB,TF,MF,fBF,BF,TM) {  
       locrelD=1
       mainSizes=tabulate(MF)
       nSamp=pmin(rep(8,nlevels(MF)), mainSizes)
       repeat {
         kmax=1
         for (k in 1: nlevels(MF)) {
           s=sort(sample(seq_len(length(TF))[MF==levels(MF)[k]], nSamp[k])) 
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
       list(VTT=VTT,VBB=VBB,VTB=VTB,TF=TF,locrelD=locrelD,TM=TM)
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
  Optimise=function(TF,MF,BF,VTT,VBB,VTB,TM,BM) {
    fBF=interaction(MF:BF)
    globrelD=0
    relD=1
    globTF=TF
    blocksizes=tabulate(fBF)
    if (regReps & max(blocksizes)==min(blocksizes) )
      blocksEffBound=upper_bounds(length(TF),nlevels(TF),nlevels(fBF)) 
    else blocksEffBound=1
    for (r in 1:searches) {
      dmax=DMax(VTT,VBB,VTB,TF,MF,fBF,BF,TM)
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
          if (!is.na(blocksEffBound)) {
            reff=EstEffics(globTF,fBF)[2]
            if (isTRUE(all.equal(blocksEffBound,reff))) return(globTF)
          }  
        }
      }  
      if (r<searches) { 
        for (iswap in 1:jumps) {
          counter=0
          repeat {
            counter=counter+1
            s1=sample(seq_len(nunits),1)
            z= seq_len(nunits)[MF==MF[s1] & fBF!=fBF[s1] & TF!=TF[s1] ] 
            if (length(z)==0) next
            if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z)
            v11=VTT[c(TF[s[1]],TF[s[2]]),c(TF[s[1]],TF[s[2]])]
            v22=VBB[c(fBF[s[1]],fBF[s[2]]),c(fBF[s[1]],fBF[s[2]])]
            v12=VTB[c(TF[s[1]],TF[s[2]]),c(fBF[s[1]],fBF[s[2]])]
            Dswap=(1-2*sum(diag(v12))+sum(v12))**2-(2*sum(diag(v22))-sum(v22))*(2*sum(diag(v11))-sum(v11))
            if (Dswap>.1| counter>1000) break
          }
          if (counter>1000) return(globTF) # finish with no non-singular swaps
          relD=relD*Dswap 
          up=UpDate(VTT,VBB,VTB,TF[s[1]],TF[s[2]], fBF[s[1]], fBF[s[2]])
          VTT=up$VTT
          VBB=up$VBB
          VTB=up$VTB
          TF[c(s[1],s[2])]=TF[c(s[2],s[1])]  
          TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
        } #end of jumps
      }
    }
    return(globTF)
  }
  # ******************************************************************************************************************************************************** 
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank) {
    candidates=NULL
    while (is.null(candidates)) {
      if (rank<(nunits-1)) s1=sample(pivot[(1+rank):nunits],1) else s1=pivot[nunits]
      candidates = (1:nunits)[ MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1] ]
    }
    if ( length(candidates)>1 ) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  }
  # ********************************************************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************
  NonSingular=function(TF,MF,BF,TM,BM) {
    fullrank=ncol(cbind(TM,BM))
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:100) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      s=Swaps(TF,MF,BF,pivot,rank)
      tindex=1:length(TF)
      tindex[c(s[1],s[2])]=tindex[c(s[2],s[1])]
      Q=qr(t(cbind(BM,TM[tindex,])))
      if (Q$rank>rank) {
        TF=TF[tindex]
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
  blocksOpt=function(TF,MF,BF) {
      TM=model.matrix(as.formula("~Treatments"),data.frame("Treatments"=TF))[,-1,drop=FALSE] # drops mean contrast
      TM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(TM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)}))
      if (nlevels(MF)==1) BM=model.matrix(as.formula(~BF))[,-1,drop=FALSE] else BM=model.matrix(as.formula(~MF+MF:BF))[,-c(1:nlevels(MF)),drop=FALSE]
      BM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(BM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)})) #centres sub-blocks within main blocks
      BM=BM[ ,as.numeric(matrix(1:(nlevels(BF)*nlevels(MF)),nrow=nlevels(BF),ncol=nlevels(MF),byrow=TRUE)[-nlevels(BF),]) ,drop=FALSE] # reorder within main blocks
      
      if ( (ncol(cbind(TM,BM))+1) > nrow(TM)  ) stop( paste("Too many parameters: plots df = ",nrow(TM)-1," model:df = ",ncol(cbind(TM,BM)) )) 
      nonsing=NonSingular(TF,MF,BF,TM,BM)
      TF=nonsing$TF
      TM=nonsing$TM
      V=chol2inv(chol(crossprod(cbind(TM,BM))))
      VTT = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
      VBB = V[ (ncol(TM)+1):ncol(V) , (ncol(TM)+1):ncol(V),drop=FALSE]
      VTB = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
      VTT=rbind ( cbind( VTT, rep(0,ncol(TM) )) , rep(0,(1+ncol(TM)) ))
      
      VBB=rbind(  cbind( VBB, matrix(0, nrow=ncol(BM),ncol=nlevels(MF)) ),  matrix(0, nrow=nlevels(MF),ncol=(ncol(BM)+nlevels(MF))))
      VTB=rbind(  cbind( VTB, matrix(0, nrow=ncol(TM),ncol=nlevels(MF)) ),  rep(0,(ncol(BM)+nlevels(MF)  )))
      
      reorder=c(rbind( matrix(seq_len(ncol(BM)), nrow=(nlevels(BF)-1), ncol=nlevels(MF)), seq(ncol(BM)+1, nlevels(BF)*nlevels(MF) )))
      VBB=VBB[reorder,reorder] 
      VTB=VTB[,reorder] 
      TF=Optimise(TF,MF,BF,VTT,VBB,VTB,TM,BM)
    return(TF)
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
    names =unlist(lapply(1:strata, function(j) {paste0("Level_",j)}))
    blocklevs=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocklevs,effics,bounds))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    return(efficiencies)
  }
 
  # ********************************************************************************************************************************************************
  # Main design program which tests input variables, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************
  hcf=HCF(replicates) 
  if (is.null(blocks)) blocks=hcf
  if (anyNA(blocks)|any(is.nan(blocks))|any(!is.finite(blocks))|any(blocks%%1!=0)|any(blocks<1)|is.null(blocks)) stop(" blocks invalid")
  fulltrts=treatments
  treatments=factor(unlist(lapply(1:hcf,function(i){sample(rep(1:sum(treatments), rep(replicates/hcf,treatments)))})))
  nunits=length(treatments)
  strata=length(blocks)
    # omit any single replicate treatments and find hcf
  if ( min(replicates)==1 & max(replicates)>1) {
    replications=rep(replicates,fulltrts)
    newlevs=(1:sum(fulltrts))[replications>1]
    hcf=HCF(replicates[replicates>1])
    treatments=factor(unlist(lapply(1:hcf,function(i){sample(rep(1:sum(fulltrts[replicates>1]), rep(replicates[replicates>1]/hcf,fulltrts[replicates>1])))})))
    levels(treatments) = newlevs
    nunits=length(treatments)
  }
  blocksizes=nunits
    for (i in 1:length(blocks))
      blocksizes=unlist(lapply(1:length(blocksizes), function(j) {
        nestblocksizes=rep(blocksizes[j]%/%blocks[i],blocks[i])
        blocksizes=nestblocksizes+rep(c(1,0),c(blocksizes[j]-sum(nestblocksizes), blocks[i]-blocksizes[j]+sum(nestblocksizes)))
      }))
  
  regBlocks=isTRUE(all.equal(max(blocksizes), min(blocksizes)))
  if (is.null(searches)) 
    if (nunits<1000) searches=10000%/%nunits else if (nunits<5000) searches=5000%/%nunits else searches=1
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  regReps=isTRUE(all.equal(max(replicates), min(replicates)))
  if (prod(blocks)*2>nunits) stop("Too many blocks for the available plots  - each block must contain at least two plots")
 # nested factor level data frames
  nestblocks=c(1,blocks)
  cumblocks=c(1,cumprod(blocks))
  nblkdesign=as.data.frame(lapply(1:(length(blocks)+1),function(i) {gl(nestblocks[i],cumblocks[length(blocks)+1]/cumblocks[i],
                                                                       labels=unlist(lapply(1:nestblocks[i], function(j) {paste0("Blocks_",j)})))}))
  fblkdesign=as.data.frame(lapply(1:(length(blocks)+1),function(i) {gl(cumblocks[i],cumblocks[length(blocks)+1]/cumblocks[i],
                                                                       labels=unlist(lapply(1:cumblocks[i], function(j) {paste0("Blocks_",j)})))}))
  colnames(nblkdesign)=unlist(lapply(1:ncol(fblkdesign), function(j) {paste0("Level_",j-1)}))
  colnames(fblkdesign)=unlist(lapply(1:ncol(fblkdesign), function(j) {paste0("Level_",j-1)}))
  nblkDesign=nblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  fblkDesign=fblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  
  # lattice squares
  v=sqrt(nlevels(treatments))  # dimension of a lattice square
  k=nunits/cumblocks[strata+1]  # average block size
  orthoMain=(regReps & identical(replicates[1],blocks[1]))
  Lattice=(regReps & regBlocks & orthoMain & identical(v,floor(v)) & identical(k,v) & length(blocks)==2)
  if (Lattice) TF=lattice(v,replicates[1]) else TF=NULL
  attempts=0
  CRB= ( length(blocks)==1 & hcf%%blocks[1]==0)
  while (is.null(TF) & attempts<10) {
    attempts=attempts+1
    TF=treatments
    if (!CRB)
      for ( i in 1:strata) 
        if (hcf%%prod(blocks[1:i])!=0) TF=blocksOpt(TF,fblkDesign[,i],nblkDesign[,i+1]) 
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  # add back single rep treatments for nested blocks only
  if (min(replicates)==1 & max(replicates)>1) {
    replications=rep(replicates,fulltrts)
    nunits=sum(replicates*fulltrts)
    fullblocksizes=nunits
    for (i in 1:length(blocks))
      fullblocksizes=unlist(lapply(1:length(fullblocksizes), function(j) {
        nestblocksizes=rep(fullblocksizes[j]%/%blocks[i],blocks[i])
        fullblocksizes=nestblocksizes+rep(c(1,0),c(fullblocksizes[j]-sum(nestblocksizes), blocks[i]-fullblocksizes[j]+sum(nestblocksizes)))
      }))
    TrtsInBlocks= split( levels(TF)[TF] , rep(1:length(blocksizes),blocksizes))
    singleTF=split( sample( (seq_len(sum(fulltrts)))[replications==1] ), rep(1:length(fullblocksizes),(fullblocksizes-blocksizes)))
    for (i in names(singleTF)) 
      TrtsInBlocks[[i]]=(append( TrtsInBlocks[[i]],singleTF[[i]]))
    TF=unlist(TrtsInBlocks)
    blocksizes=fullblocksizes
  }
    fblkDesign=cbind(fblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE],Plots=factor(1:nunits))
    fblkDesign=data.frame(lapply(1:ncol(fblkDesign), function(r){sample(nlevels(fblkDesign[,r]))[fblkDesign[,r]]})) # Randomize labels - gives numeric columns
    fblkDesign=cbind(fblkDesign,TF)
    fblkDesign=fblkDesign[do.call(order, fblkDesign), ] # re-order block and plot labels - randomizes treatments within nested blocks
    fblkDesign[]=lapply(fblkDesign, as.factor)
    blocksizes=table(fblkDesign[,ncol(fblkDesign)-2])[unique(fblkDesign[,ncol(fblkDesign)-2])] # randomized block sizes
    fDesign=data.frame(fblkdesign[ rep(1:length(blocksizes),blocksizes),-1,drop=FALSE],factor(1:nunits),fblkDesign[,ncol(fblkDesign)])
    Design =data.frame(nblkdesign[rep(1:length(blocksizes),blocksizes),-1,drop=FALSE],factor(1:nunits),fblkDesign[,ncol(fblkDesign)])
    Efficiencies=BlockEfficiencies(fDesign) 
    colnames(Design)=c( colnames(nblkdesign)[-1],"Plots","Treatments")
    V = split(Design[,ncol(Design)],fDesign[,(ncol(fDesign)-2)])
    V = lapply(V, function(x){ length(x) =max(blocksizes); x })
    Plan = data.frame(nblkdesign[,-1 ,drop=FALSE],rep("",length(V)),matrix(unlist(V),nrow=length(V),byrow=TRUE))
    colnames(Plan)=c(colnames(nblkdesign[,-1 ,drop=FALSE]),"Blocks.Plots:", c(1:max(blocksizes)))
  row.names(Plan)=NULL
  row.names(Design)=NULL
  row.names(Efficiencies)=NULL
  # treatment replications
  TreatmentsTable=as.data.frame(table(Design[,ncol(Design)]))
  TreatmentsTable=TreatmentsTable[order( as.numeric(levels(TreatmentsTable[,1])) [TreatmentsTable[,1]]),]
  TreatmentsTable[]=lapply(TreatmentsTable, as.factor)
  colnames(TreatmentsTable)=c("Treatments","Replicates")
  row.names(TreatmentsTable)=NULL
  list(Treatments=TreatmentsTable,BlocksEfficiency=Efficiencies,Plan=Plan,Design=Design,seed=seed,searches=searches,jumps=jumps)
}
