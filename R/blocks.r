#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested and crossed block designs for unstructured treatment sets where treatments can have any arbitrary levels of replication
#' and blocks can have any arbitrary feasible depth of nesting.
#' 
#' @details
#' 
#' The \code{blocksdesign} package constructs nested or crossed block designs with arbitrary number of nested strata and arbitrary nesting or crossing
#' in each blocks stratum.
#' 
#' The \code{treatments} and \code{replicates} vectors give a partition of the total number of treatments into sets of equally replicated treatments where the
#' \code{treatments} vector provides the required treatments partition and the \code{replicates} vector provides the required replication for each treatment set. 
#' The sum of the cross-products of the \code{treatments} and \code{replicates} vectors is the total number of plots and treatments are numbered consecutively 
#' according to the ordering of the consecutive treatment sets. 
#'   
#' The \code{rows} vector, if specified, defines the nested row blocks in each nested stratum taken in order from the highest to the lowest.
#' The first number, if any, is the number of rows in the blocks of the top stratum, the second, if any, is the number of rows in the nested blocks of
#'  the second stratum, the third, if any, is the number of rows in the nested blocks of the third stratum and so on for all the required strata in the design. 
#'  
#' The \code{columns} vector, if specified, defines the nested column blocks in each nested stratum taken in order from the highest to the lowest. 
#' The first number, if any, is the number of columns in the blocks of the top stratum, the second, if any, is the number of columns in the nested blocks of
#'  the second stratum, the third, if any, is the number of columns in the nested blocks of the third stratum and so on for all the required strata in the design. 
#'   
#' The \code{rows} and \code{columns} vectors, if defined, must be of equal length and if a simple nested blocks design is required in 
#' any particular stratum then the number of columns in that stratum should be set to unity.
#' 
#' If both the \code{rows} vector and the \code{columns} vector are null, the default block design will be a single set of orthogonal
#' main blocks equal in number to the highest common factor of the replication numbers. If the \code{rows} vector is defined but the \code{columns} vector
#' is null, the design will comprise simple nested blocks in each stratum defined by the \code{rows} vector.
#' 
#' The block sizes in any particular rows or columns stratum are always as nearly equal in size as possible and will never differ by more
#' than a single unit. The row blocks and column blocks must always comprise at least two plots and this restriction will constrain the 
#' allowed number of row and column blocks.
#' 
#' Unreplicated treatments with a single replication are allowed but unreplicated treatments are added heuristically after the replicated
#' part of the design has been optimized and to ensure that the blocks containing the unreplicated treatments do not differ in size by more than
#'  a single unit, unreplicated treatments should only be used in designs where the replicated treatment blocks 
#'  are all equal in size in each stratum of the design.     
#' 
#'  NB. For two-replicate 2 x 2 row-and-column block designs, the algorithm will always find a semi-Latin squarere with complete
#'  replicate rows and complete replicate columns but will always fully confound one treatment contrast with the row-by-column interaction contrast.
#'  For this reason, the algorithm is unable to fit any nested sub-block design within the row-by-column intersections of any two replicate 2 x 2 semi-Latin square. 
#'  A better design would be a simple nested block design with four main blocks and sub-blocks nested within the main blocks.
#'  
#'  Lattice designs where v is a prime-power require the \code{\link[crossdes]{MOLS}} package.
#' 
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order.  \cr
#'  \item  A table showing the replication number of each treatment in the design . \cr
#'  \item  A table showing the block levels and the achieved D- and A-efficiency factors for each blocks stratum together with A-efficiency upper bounds, where available. \cr
#'  \item  Plans showing the allocation of treatments to blocks or to rows and columns in the bottom stratum of the design. 
#'  Each plan layout is indexed by the factor levels of any higher level blocks containing the bottom level plans.    \cr
#' } 
#' 
#' @param treatments  a vector giving a partition of the total required number of treatments into sets of equally replicated treatments.
#' 
#' @param replicates  a vector giving the replication number of each equally replicated treatment set in the \code{treatments} vector.
#' 
#' @param rows  a vector of factor levels for the row blocks in each succesive stratum of the blocks design taken in order from the highest to the lowest. 
#' The default is a single set of main blocks equal to the hcf of the replication numbers.
#' 
#' @param columns  a vector of factor levels forthe column blocks in each succesive stratum of the blocks design taken in order from the highest to the lowest. 
#' The \code{rows} and the \code{columns} vectors, if both present, must be of equal length. The default is the null vector.  
#' 
#' @param seed  an integer initializing the random number generator. The default is a random seed.
#' 
#' @param searches  the maximum number of local optima searched for a design optimization. The default is 1 plus the floor of 2000 divided by the number of model parameters.
#' 
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#' 
#' @return  
#' \item{Design}{Data frame giving the optimized block and treatment factors in plot order}
#' \item{Plan}{Data frame giving a plan view of the treatments design in the bottom stratum of the design classified by rows and columns}
#' \item{Efficiencies}{The achieved A- and D-efficiencies for each stratum of the design together with an A-efficiency upper-bound, where available}
#' \item{seed}{Numerical seed for random number generator}
#' \item{searches}{Maximum number of searches in each stratum}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima}
#' 
#' @references
#' 
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. http://CRAN.R-project.org/package=crossdes
#' 
#' @examples
#' 
#' # 3 treatments x 2 replicates, 2 treatments x 4 replicates and 4 treatments x 3 replicates  
#' # the hcf of the replication numbers is 1 therefore the default design is completely randomized 
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 4 treatments x 4 replicates with 2 main rows each containing two complete replicates  
#' blocks(treatments=4,replicates=4,rows=2)
#' 
#' # 50 treatments x 4 replicates with 4 main rows and 5 nested sub-rows in each main block 
#' blocks(treatments=50,replicates=4,rows=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments 
#' # giving exactly one single replicate treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),rows=c(4,5))
#' 
#' # 64 treatments x 2 replicates with 2 main rows and five succesively nested 2-level factors
#' blocks(treatments=64,replicates=2,rows=c(2,2,2,2,2,2))
#' 
#' # 6 replicates of 6 treatments in 4 rows of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,rows=4)
#' 
#' # 4 replicates of 13 treatments arranged in a 13 x 4 Youden rectangle 
#' blocks(treatments=13,replicates=4,rows=13,columns=4)
#' 
#' # 64 treatments x 2 replicates with nested 8 x 8 row-and-column designs in two main blocks 
#' blocks(treatments=64,replicates=2,rows=c(2,8),columns=c(1,8)) 
#' 
#' # 2**9 treatments x 2 replicates in 2**9 blocks giving a fully saturated block design 
#' # (requires a considerable time to run!)
#' \dontrun{ d=blocks(2**9,2,rep(2,9)) }
#'          
#' @export
#' @importFrom stats anova lm
#' 
blocks = function(treatments,replicates,rows=HCF(replicates),columns=NULL,searches=(1+2000%/%(sum(treatments)+prod(rows))),seed=sample(10000,1),jumps=1) { 
  
  # ******************************************************************************************************************************************************** 
  # Finds the highest common factor (hcf) of a set of numbers omitting any zero values (Euclidean algorithm)
  # ********************************************************************************************************************************************************
  HCF=function(replicates)  {
    replicates=sort(replicates[replicates>0])
    for (i in  seq_len(length(replicates)))
      while (!isTRUE(all.equal(replicates[i]%%replicates[1],0))) replicates[c(1,i)] = c(replicates[i]%%replicates[1], replicates[1])
      replicates[1]
  }   
  # ******************************************************************************************************************************************************** 
  # Tests a given number for primality and returns TRUE or FALSE
  # ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v < 4) return(TRUE) 
    if ( isTRUE(all.equal(v %% 2,0)) ||  isTRUE(all.equal(v %% 3,0)) ) return(FALSE) 
    if (v<25) return(TRUE)
    for(i in  6*seq_len(length(floor((sqrt(v)+1)/6)))        )
      if ( isTRUE(all.equal(v %% (i-1) , 0)) ||   isTRUE(all.equal(v %% (i+1) , 0)) ) return(FALSE) 
    return(TRUE)
  } 
  # ******************************************************************************************************************************************************** 
  # Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only  
  # ********************************************************************************************************************************************************
  Contrasts=function(MF,NF) {
    NM=matrix(0,nrow=length(NF),ncol=nlevels(NF))
    NM[cbind(seq_len(length(NF)),NF)]=1 # factor indicator matrix  
    for (i in seq_len(nlevels(MF))) 
      NM[MF==i,]=scale(NM[MF==i,] , center = TRUE, scale = FALSE)
    NM
  }
  # ******************************************************************************************************************************************************** 
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb   
  # ********************************************************************************************************************************************************
  UpDate=function(MTT,MBB,MTB,ti,tj,bi,bj) {  
    mtt=MTT[ti,ti]+MTT[tj,tj]-2*MTT[tj,ti]
    mbb=MBB[bi,bi]+MBB[bj,bj]-2*MBB[bi,bj]
    mtb=1-MTB[ti,bi]+MTB[tj,bi]+MTB[ti,bj]-MTB[tj,bj]  
    TBbij=MTB[,bi]-MTB[,bj]
    TBtij=MTB[ti,]-MTB[tj,]
    TTtij=MTT[,ti]-MTT[,tj]
    BBbij=MBB[bi,]-MBB[bj,]
    Z1 = (TBbij-TTtij)/sqrt(2*mtb+mtt+mbb)   
    Z2 = (BBbij-TBtij)/sqrt(2*mtb+mtt+mbb)
    W1 = ( sqrt(2*mtb+mtt+mbb)*(TTtij+TBbij) - (mbb-mtt)*Z1) /(2*sqrt(mtb**2-mtt*mbb))
    W2 = ( sqrt(2*mtb+mtt+mbb)*(TBtij+BBbij) - (mbb-mtt)*Z2) /(2*sqrt(mtb**2-mtt*mbb))
    MTT = MTT - tcrossprod(Z1) + tcrossprod(W1)
    MBB = MBB - tcrossprod(Z2) + tcrossprod(W2)
    MTB = MTB - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }  
  # ******************************************************************************************************************************************************** 
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  optEffics=function(TF,BF,ntrts,k) { 
    if (ntrts<=k) 
      e=eigen( (diag(ntrts)-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)] else    
        e=c(rep(1,(ntrts-k)), 
            eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])  
      round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6)
  }
  # ******************************************************************************************************************************************************** 
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
  # ********************************************************************************************************************************************************
  DMax=function(MTT,MBB,MTB,TF,BF,Restrict) {   
    relD=1
    mainSizes=tabulate(Restrict)
    nSamp=pmin(rep(8,nlevels(Restrict)),mainSizes)
    repeat {
      improved=FALSE
      for (k in seq_len(nlevels(Restrict))) {
        S=sort(sample(seq_len(length(TF))[Restrict==k], nSamp[k])) 
        TB=MTB[TF[S],BF[S],drop=FALSE]-tcrossprod(MTB[cbind(TF[S],BF[S])],rep(1,nSamp[k]))
        dMat=(TB+t(TB)+1)**2-
          (2*MTT[TF[S],TF[S],drop=FALSE]-tcrossprod(MTT[cbind(TF[S],TF[S])]+rep(1,nSamp[k]) ) + tcrossprod(MTT[cbind(TF[S],TF[S])]) + 1)*
          (2*MBB[BF[S],BF[S],drop=FALSE]-tcrossprod(MBB[cbind(BF[S],BF[S])]+rep(1,nSamp[k]) ) + tcrossprod(MBB[cbind(BF[S],BF[S])]) + 1)
        sampn=which.max(dMat) 
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        if ( !isTRUE(all.equal(dMat[i,j],1)) && dMat[i,j]>1) {
          improved=TRUE
          relD=relD*dMat[i,j]
          up=UpDate(MTT,MBB,MTB,TF[S[i]],TF[S[j]],BF[S[i]],BF[S[j]])
          MTT=up$MTT
          MBB=up$MBB
          MTB=up$MTB
          TF[c(S[i],S[j])]=TF[c(S[j],S[i])]
        }
      } 
      if (improved) next
      if (sum(nSamp) < min(length(TF),512)) nSamp=pmin(mainSizes,2*nSamp) else break
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,TF=TF,relD=relD)
  }  
  # ******************************************************************************************************************************************************** 
  #  Number of searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,BF,MTT,MBB,MTB,Restrict)  {
    globrelD=0
    relD=1
    globTF=TF
    treps=tabulate(TF)
    breps=tabulate(BF)     
    if (identical(max(treps),min(treps)) && identical(max(breps),min(breps))  )
      bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) else bound=NA
    for (r in seq_len(searches)) {
      dmax=DMax(MTT,MBB,MTB,TF,BF,Restrict)
      if ( !isTRUE(all.equal(dmax$relD,1)) && dmax$relD>1) {
        relD=relD*dmax$relD
        TF=dmax$TF
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
        if (!isTRUE(all.equal(relD,globrelD)) && relD>globrelD) {
          globTF=TF
          globrelD=relD
          if ( !is.na(bound) && isTRUE(all.equal(bound,optEffics(globTF,BF,nlevels(TF),nlevels(BF))[2]))) break
        }
      }
      if (r==searches) break
      for (iswap in seq_len(jumps)) {
        counter=0
        repeat {  
          counter=counter+1
          s1=sample(seq_len(length(TF)),1)
          z=seq_len(length(TF))[Restrict==Restrict[s1] & BF!=BF[s1] & TF!=TF[s1]]
          if (length(z)==0) next
          if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z[1])
          dswap = (1+MTB[TF[s[1]],BF[s[2]]]+MTB[TF[s[2]],BF[s[1]]]-MTB[TF[s[1]],BF[s[1]]]-MTB[TF[s[2]],BF[s[2]]])**2-
            (2*MTT[TF[s[1]],TF[s[2]]]-MTT[TF[s[1]],TF[s[1]]]-MTT[TF[s[2]],TF[s[2]]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]])  
          if (dswap>.1 | counter>1000) break
        }
        if (counter>1000) return(globTF) # no non-singular swaps
        relD=relD*dswap
        up=UpDate(MTT,MBB,MTB,TF[s[1]],TF[s[2]], BF[s[1]], BF[s[2]])
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TF[c(s[1],s[2])]=TF[c(s[2],s[1])]  
      } 
    }
    globTF
  } 
  # ******************************************************************************************************************************************************** 
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank,n,Restrict) {
    candidates=NULL
    while (isTRUE(all.equal(length(candidates),0))) {
      if (rank<(n-1)) s1=sample(pivot[(1+rank):n],1) else s1=pivot[n]
      candidates = seq_len(n)[MF==MF[s1] & Restrict==Restrict[s1] & BF!=BF[s1] & TF!=TF[s1]]
    }
    if ( length(candidates)>1 )
      s2=sample(candidates,1) else s2=candidates[1] 
      s=c(s1,s2)
  }
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  NonSingular=function(TF,MF,BF,Restrict) { 
    BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
    BM[cbind(1:length(BF),BF)]=1
    fullrank=nlevels(TF)+nlevels(BF)-1
    TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
    TM[cbind(seq_len(length(TF)),TF)]=1
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    times=0
    n=length(TF)
    while (rank<fullrank & times<1000) {
      times=times+1
      s=Swaps(TF,MF,BF,pivot,rank,n,Restrict)
      rindex=seq_len(length(TF))
      rindex[c(s[1],s[2])]=rindex[c(s[2],s[1])]
      newQ=qr(t(cbind(BM,TM[rindex,])))
      if (isTRUE(all.equal(newQ$rank,rank)) || newQ$rank>rank) { 
        TF=TF[rindex]
        TM=TM[rindex,]
        rank=newQ$rank
        pivot=newQ$pivot
      } 
    }
    if (times>999) TF=NULL
    return(TF)
  }  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  colsOpt=function(TF,Main,Columns,Blocks) { 
    TF=NonSingular(TF,Main,Columns,Blocks)
    if (is.null(TF)) return(TF)
    blevels=nlevels(Columns)%/%nlevels(Main)
    BM=Contrasts(Main,Columns)[, rep(c(rep(TRUE,(blevels-1)),FALSE),nlevels(Main)),drop=FALSE]
    TM=Contrasts(Main,TF)[,-nlevels(TF),drop=FALSE] 
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    MBB=matrix(0,nrow=nlevels(Columns),ncol=nlevels(Columns))
    MTT=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))  
    MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(Columns))
    indicb=seq_len(nlevels(Columns)-nlevels(Main))
    indict=seq_len(nlevels(TF)-1)
    indicv=seq_len(ncol(V)-nlevels(Columns)+nlevels(Main)) + nlevels(Columns)-nlevels(Main)
    MBB[indicb,indicb]=V[indicb,indicb,drop=FALSE]
    MTT[indict,indict]=V[indicv,indicv,drop=FALSE]
    MTB[indict,indicb]=V[indicv,indicb,drop=FALSE]
    perm=order(order(seq_len(nlevels(Columns))%%blevels ==0  ))  
    MTB=MTB[,perm]
    MBB=MBB[perm,perm] 
    TF=Optimise(TF,Columns,MTT,MBB,MTB,Blocks)
    return(TF)
  }  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ******************************************************************************************************************************************************** 
  rowsOpt=function(TF,Main,Blocks) { 
    nblocks=nlevels(Blocks)
    mblocks=nlevels(Main)
    v=sqrt(ntrts)  
    regrep=isTRUE(all.equal(max(replicates),min(replicates)))
    regular=isTRUE(all.equal(nunits%%nblocks,0)) && isTRUE(all.equal(hcf%%mblocks,0))
    sqrLattice  = regular && isTRUE(all.equal(v,floor(v))) && isTRUE(all.equal(nunits,v*nblocks))
    w=sqrt(nblocks)
    s=ntrts/w
    trojanic = regular && isTRUE(all.equal(replicates[1],w)) && isTRUE(all.equal(s,floor(s))) && (s<w)
    if (sqrLattice  && replicates[1]<4) {
      t=c(rep(0:(v-1),each=v),rep(0:(v-1),v)+v)
      if (replicates[1]>2)
        for (i in 0: (v-1)) 
          t=c(t,(rep(0:(v-1))+i)%%v + 2*v)
        TF=  as.factor(rep(seq_len(v*v),replicates[1])[order(t)])
    } else if (sqrLattice  &&  replicates[1]<(v+2)  && isPrime(v) ) {
      t=c(rep(0:(v-1),each=v),rep(0:(v-1),v)+v)
      for (z in seq_len(replicates[1]-2))
        for (j in 0: (v-1)) 
          t=c(t,(rep(0:(v-1))*z +j)%%v + v*(z+1) )
        TF=as.factor(rep(seq_len(v*v),replicates[1])[order(t)])
    } else if (sqrLattice  && replicates[1]<(v+2)  &&  ntrts%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==ntrts)
      mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])			
      TF=c(seq_len(ntrts),seq_len(ntrts)[order(rep(0:(v-1),v))])
      for (i in seq_len(replicates[1]-2))
        TF=c(TF, seq_len(ntrts)[order(as.numeric(mols[,,i]))]) 
      TF=as.factor(TF)
    } else if (sqrLattice  && v==10  && replicates[1]==4) {
      TF=as.factor(c(
        seq_len(100),seq_len(100)[order(rep(0:9,10))],
        seq_len(100)[order(c(
          1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
          5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8))],
        seq_len(100)[order(c(
          1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
          6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6))]
      )) 
    } else if (trojanic  && isPrime(w) ) {
      for (z in 1:(nunits/nblocks))
        for (j in 0: (w-1)) 
          for (k in 0: (w-1)) 
            TF[ (k + j*w)*nunits/nblocks + z]=(j+k*z)%%w + (z-1)*w +1
          TF=as.factor(TF)
    } else if (trojanic &&  nblocks%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==nblocks)
      mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])		
      TF=NULL
      for (i in  seq_len(w))
        for (j in seq_len(w))
          for (k in seq_len(nunits/nblocks))
            TF=c(TF,mols[i,j,k]+(k-1)*w)
      TF=as.factor(TF)
    } else {
      TF=NonSingular(TF,Main,Blocks,rep(1,length(TF)))
      if (is.null(TF)) return(TF)
      blevels=nlevels(Blocks)%/%nlevels(Main)
      BM=Contrasts(Main,Blocks)[, rep(c(rep(TRUE,(blevels-1)),FALSE),nlevels(Main)),drop=FALSE]
      TM=Contrasts(Main,TF)[,-nlevels(TF),drop=FALSE] 
      V=chol2inv(chol(crossprod(cbind(BM,TM))))
      MBB=matrix(0,nrow=nlevels(Blocks),ncol=nlevels(Blocks))
      MTT=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))  
      MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(Blocks))
      indicb=seq_len(nlevels(Blocks)-nlevels(Main))
      indict=seq_len(nlevels(TF)-1)
      indicv=seq_len(ncol(V)-nlevels(Blocks)+nlevels(Main)) + nlevels(Blocks)-nlevels(Main)
      MBB[indicb,indicb]=V[indicb,indicb,drop=FALSE]
      MTT[indict,indict]=V[indicv,indicv,drop=FALSE]
      MTB[indict,indicb]=V[indicv,indicb,drop=FALSE]
      perm=order(order(seq(nlevels(Blocks))%%blevels ==0 ))  
      MTB=MTB[,perm]
      MBB=MBB[perm,perm] 
      TF=Optimise(TF,Blocks,MTT,MBB,MTB,Main)
    }
    return(TF)
  }  
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for the rows in each stratum of a design 
  # ********************************************************************************************************************************************************     
  A_Efficiencies=function(Design) {
    levels=as.numeric(rbind(rows,columns))
    ntrts=nlevels(Design[,ncol(Design)])
    effics=matrix(1,nrow=2*strata,ncol=2)
    bounds=rep(NA,2*strata)
    nblocks=unlist(lapply(1:(2*strata), function(i) {nlevels(Design[,i])})) 
    
    if (isTRUE(all.equal(max(replicates),min(replicates))) ) 
      for (i in seq_len(2*strata))  
        if ( isTRUE(all.equal(nunits%%nblocks[i],0))) 
          bounds[i]=upper_bounds(nunits,ntrts,nblocks[i]) else if (hcf%%nblocks[i]==0) bounds[i]=1
      
    for (i in seq_len(2*strata))    
      if (ntrts>1 && nblocks[i]>1)
        effics[i,]=optEffics(Design$Treatments,Design[,i],ntrts,nblocks[i])  

    names=colnames(Design)
    efficiencies=data.frame(cbind( names[1:(2*strata)],levels,nblocks,effics, bounds)) 
    colnames(efficiencies)=c("Stratum","Levels","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies[,'Blocks'] = as.factor(efficiencies[,'Blocks'])
    # omit any null strata
      efficiencies=efficiencies[(levels>1),,drop=FALSE]
    efficiencies
  }
  
  # ******************************************************************************************************************************************************** 
  # Carries out some input validation
  # ********************************************************************************************************************************************************     
  testInputs=function(treatments,replicates,rows,columns,seed) {  
    if (missing(treatments) | missing(replicates)) return(" Treatments or replicates not defined ")   
    if (is.null(treatments) | is.null(replicates)) return(" Treatments or replicates list is empty ")   
    if (anyNA(treatments) | anyNA(replicates)) return(" NA values not allowed")
    if (any(!is.finite(treatments)) | any(!is.finite(replicates)) | any(is.nan(treatments)) | any(is.nan(replicates))) 
      return(" Treatments and replicates can contain only finite integers ")
    if ( length(treatments)!=length(replicates)) 
      return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
    if (any(treatments<1)) return("Treatments must be non-negative integers")
    if (any(replicates<1)) return("Replicates must be non-negative integers")  
    if (anyNA(rows) ) return(" NA rows values not allowed") 
    if (any(!is.finite(rows)) | any(is.nan(rows)) ) return(" rows can contain only finite integers ")
    if (min(rows)<1) return (" rows must be at least one ")
    if (anyNA(columns) ) return(" NA columns values not allowed") 
    if (any(!is.finite(columns)) | any(is.nan(columns)) ) return(" columns can contain only finite integers ")
    if (min(columns)<1) return (" columns must be at least one ")
    cumcols=cumprod(columns)
    cumrows=cumprod(rows)
    cumblocks=c(1,cumprod(rows*columns))
    strata=length(rows)
    plots=sum(treatments[replicates>1]*replicates[replicates>1])
    if (cumblocks[strata]*rows[strata]*2>plots) return("Too many rows for the available plots  - every row must contain at least two replicated plots")
    if (cumblocks[strata]*columns[strata]*2>plots) return("Too many columns for the available plots  - every column must contain at least two replicated plots")
    if (cumblocks[strata+1]>plots) return("Too many rows for the available plots  - every row-by-column intersection must contain at least one replicated plot")
    if (!is.null(searches)) {
      if (anyNA(searches) ) return(" NA searches values not allowed") 
      if ( any(!is.finite(searches)) | any(is.nan(searches))) return(" Searches must be a finite integer ") 
      if (searches<1)  return(" Repeats must be at least one ")   
    }  
    if (!is.null(jumps)) {
      if (anyNA(jumps) ) return(" NA jumps values not allowed") 
      if ( !all(is.finite(jumps)) | !all(!is.nan(jumps))) return(" jumps must be a finite integer ") 
      if (jumps<1)  return(" Random jumps must be at least one ")   
    }    
    if (!is.null(seed)) {
      if (anyNA(seed) ) return(" NA seed values not allowed") 
      if (any(!is.finite(seed)) | any(is.nan(seed))) return(" Seed must be a finite integer ") 
      if (seed<1)  return(" Seed must be at least one ")   
    } 
    if ( isTRUE( sum(treatments) < 2 ) )
      return(paste("The number of treatments must be at least two "))  
    
    if ( isTRUE( sum(treatments*replicates) < (prod(rows) + sum(treatments)-1) ) )
      return(paste("The total number of plots is",  sum(treatments*replicates) , 
                   "whereas the total required number of model parameters is", prod(rows) + sum(treatments),", which is not feasible. "))  
    
    if (max(replicates)==2 && length(rows)>1 )
      for (i in seq_len(length(rows)-1)) 
        if (rows[i]==2 && columns[i]==2) return( paste(" Warning: this design tries to fit a sub-block design within the blocks of a 2-replicate
        semi-Latin square. However, because the row-by-column interaction of a 2-replicate semi-Latin square is always confounded with a treatment contrast, 
        the nested sub-blocks design will always be singular. Instead, try replacing the crossed 2 x 2 row-and column block design by an uncrossed design with 4 rows and 
        a single column."))
    return(TRUE)
  }
  
  # ******************************************************************************************************************************************************** 
  # Finds row and column sizes in each stratum of a design 
  # ********************************************************************************************************************************************************     
  Sizes=function(blocksizes,stratum) {
    nblocks=length(blocksizes)
    newblocksizes=NULL
    for (j in seq_len(nblocks)) {
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
      newblocksizes=c(newblocksizes , unlist(rowcolsizes))
    }
  newblocksizes
  }
 
  # ******************************************************************************************************************************************************** 
  # Main body of rows design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************     
  if (isTRUE(all.equal(length(rows),0))) rows=1
  if (isTRUE(all.equal(length(columns),0))) columns=rep(1,length(rows))
  if (!isTRUE(all.equal(length(columns),length(rows)))) stop("The number of row rows strata and the number of column rows strata must be equal")
  testout=testInputs(treatments,replicates,rows,columns,seed) 
  if (!isTRUE(testout)) stop(testout)
  set.seed(seed)
  strata=length(rows)
  fulltreatments=treatments
  fullreplicates=replicates
  treatments=treatments[replicates>1]
  replicates=replicates[replicates>1]
  nunits=sum(treatments*replicates)
  ntrts=sum(treatments)
  hcf=HCF(replicates)
  cumblocks=c(1,cumprod(rows*columns))
  if (prod(columns)>1) 
  stratumnames=unlist(lapply(1:strata, function(i) { c(paste("Rows",i), paste("Columns", i) )})) else
    stratumnames=unlist(lapply(1:strata, function(i) { c(paste("Blocks",i), paste("Columns", i) )}))
  blocksizes=nunits
  for (i in seq_len(strata))
    blocksizes=Sizes(blocksizes,i) 
  # Design factors
  Blocks=data.frame(do.call(cbind,lapply(1:strata,function(r){ rep(1:cumblocks[r],each=(cumblocks[strata+1]/cumblocks[r]))})))
  fRows=do.call(cbind,lapply(1:strata,function(i) {(Blocks[,i]-1)*rows[i]+rep(rep(1:rows[i],each=cumblocks[strata+1]/rows[i]/cumblocks[i]),cumblocks[i])}))
  fCols=do.call(cbind,lapply(1:strata,function(i) {(Blocks[,i]-1)*columns[i]+rep(rep(1:columns[i],each=cumblocks[strata+1]/rows[i]/cumblocks[i]/columns[i]),rows[i]*cumblocks[i])}))  
  fDesign=data.frame(cbind(fRows,fCols))[, order(c( seq(1,2*strata, by=2) ,seq(2,2*strata, by=2)))]
  fDesign[]=lapply(fDesign, as.factor) 
  #permBlocks will randomise blocks in each stratum
  tempDesign=data.frame( do.call(cbind,lapply(1:(2*strata), function(r){ sample(nlevels(fDesign[,r]))[fDesign[,r]] })) , seq_len(nrow(fDesign)   ))
  permBlocks=tempDesign[ do.call(order, tempDesign), ][,ncol(tempDesign)]
  # names for plans
  if (strata>1 & max(columns)>1)  
    rcnames=unlist(lapply(1:cumblocks[strata], function(i) {c(paste(fDesign[((i-1)*rows[strata]*columns[strata]+1),],collapse= " "))})) 
  else if (strata==1 & max(columns)>1) rcnames="Rows x Columns"
  else if (strata>1 & max(columns)==1)  
    rcnames=unlist(lapply(seq_len(cumblocks[strata]), function(i) {c(paste(fDesign[((i-1)*rows[strata]+1),seq(1, 2*(strata-1),by = 2)],collapse= " "))})) 
  else  rcnames=" "
  Design    = as.data.frame(fDesign[rep(seq_len(nrow(fDesign)),  blocksizes ),])  
  Blocks    = as.data.frame(Blocks[rep(seq_len(nrow(Blocks)),  blocksizes ),]) 
  Design[]  = lapply(Design, as.factor) 
  Blocks[]  = lapply(Blocks, as.factor) 
 
  for ( z in seq_len(10)) {
    TF=rep(rep( seq_len(ntrts)  ,rep(replicates/hcf,treatments)), hcf)
    rand=sample(nunits)
    TF=as.factor( TF[rand][order(rep(seq_len(hcf),each=nunits/hcf)[rand])] )
    for ( i in seq_len(strata)) {
    if (rows[i]>1 && !all(hcf%%cumprod(rows[1:i])==0)) 
      TF=rowsOpt(TF,Blocks[,i],Design[,2*i-1])
    if (is.null(TF)) break
    if (columns[i]>1 ) 
      TF=colsOpt(TF,Blocks[,i],Design[,2*i],Design[,2*i-1])
    if (is.null(TF)) break
    }
    if (!is.null(TF)) break
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  
  #add back single rep treatments
  if ( min(fullreplicates)==1 && max(fullreplicates)>1 ) {
    addTF=((sum(treatments)+1) :sum(fulltreatments))
    TF=as.factor( c(TF, addTF ) )
    newblocksizes=length(addTF)
    for ( i in seq_len(strata)) 
      newblocksizes=Sizes(newblocksizes,i)
    reptrts=rep(fullreplicates==1,fulltreatments)
    levels(TF)= seq_len(sum(fulltreatments*fullreplicates))[order(reptrts)] 
    TF=as.numeric(levels(TF))[TF]
    nunits=sum(fulltreatments*fullreplicates)
    TF=as.factor(TF[order(c(rep( seq_len(length(blocksizes)),blocksizes),rep(   seq_len(length(blocksizes)),newblocksizes)))])
    blocksizes=blocksizes+newblocksizes
    treatments=fulltreatments
    replicates=fullreplicates
    ntrts=sum(treatments)
  }
  # Randomize
  D=as.data.frame(cbind(rep(1:length(blocksizes),blocksizes),sample(seq_len(nunits)),TF))
  D[]=lapply(D, as.factor)
  D[,1] = factor(D[,1],levels(D[,1])[permBlocks])
  TF=D[ do.call(order, D),][,ncol(D)]
  blocksizes=blocksizes[permBlocks]
   Design =as.data.frame(cbind(  fDesign[rep(seq_len(nrow(fDesign)),  blocksizes ),] ,TF) )  
  colnames(Design)=c(stratumnames,"Treatments")
  rownames(Design)=NULL
  #Plan
  V=split(Design[,(2*strata+1)],rep(1:cumblocks[strata+1],blocksizes))
  baseblocks=unlist(lapply(seq_len(length(V)), function(r){ paste( "  ", paste(format(V[[r]], width=nchar(sum(treatments)) ), collapse = " "))}))
  Plan=vector(mode="list",length=cumblocks[strata+1]/rows[strata]/columns[strata]  )
  if (columns[strata]>1) {
  colhead="Column 1"
  for (i in 2: columns[strata])
    colhead=c(colhead,paste("Column",i))
  rowhead="Row 1"
  if (rows[strata]>1)
  for (i in 2: rows[strata])
    rowhead=c(rowhead,paste("Row",i))
  } else {
    colhead="treatments"
    rowhead="Blocks 1"
    if (rows[strata]>1)
      for (i in 2: rows[strata])
        rowhead=c(rowhead,paste("Blocks",i))
  }
  for (i in seq_len(cumblocks[strata])) {
    Plan[[i]]=as.data.frame(matrix(baseblocks[((i-1)*rows[strata]*columns[strata]+1):(i*rows[strata]*columns[strata])],nrow=rows[strata],ncol=columns[strata],byrow=TRUE))
    colnames(Plan[[i]])=colhead
    rownames(Plan[[i]])=rowhead
  }
  names(Plan)=rcnames 
  # efficiencies
  Efficiencies=A_Efficiencies(Design)
  row.names(Efficiencies)=NULL
  # treatment replications
  Treatments=data.frame(table(Design[,"Treatments"]))
  Treatments[]=lapply(Treatments, as.factor) 
  colnames(Treatments)=c("Treatments","Replicates")
  
  list(Treatments=Treatments,Efficiencies=Efficiencies,Plan=Plan,Design=Design,Seed=seed,Searches=searches,Jumps=jumps) 
} 