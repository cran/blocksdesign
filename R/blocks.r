#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested block designs for unstructured treatment sets where treatments can have any arbitrary levels of replication
#' and blocks can have any arbitrary feasible depth of nesting.
#' 
#' @details
#' 
#' The \code{treatments} and \code{replicates} parameters partition the total number of treatments into sets of equally replicated treatments where 
#' \code{treatments} contains the set sizes and \code{replicates} contains the set replication numbers. 
#' The sum of the set sizes is the total number of treatments and the sum of the cross-products of the set sizes and the replication numbers
#' is the total number of plots. Treatments are numbered consecutively according to the ordering of the consecutive treatment sets. 
#' 
#' The \code{blocklevels} parameter contains the nested blocks for each individual stratum taken in order from the highest to the lowest.
#' The first number is the number of main blocks, the second, if any, is the number of sub-blocks nested in each main block, the third, if any, 
#' is the number of sub-sub-blocks nested in each sub-block,and so on for all the reqired strata. If left blank, the default block design is a 
#'  maximal set of orthogonal main blocks, where the maximal number of of orthogonal main blocks is the highest common factor of the replication numbers. 
#'
#' The block sizes in each blocks stratum are always as equal as possible and never differ by more than a single unit. If the number of nested blocks 
#' in a particular stratum exactly divdes the number of units, the block sizes in that stratum will be exactly equal; otherwise the block sizes
#' will differ by at most one unit.
#' 
#'  Special square and rectangular lattice designs are constructed algebraically and include all designs that can be constructed from 
#'  a single latin square, from a complete sets of prime or prime-power orthogonal latin squares or from a pair of orthogonal 10 x 10 Latin squares. 
#'  All other non-orthogonal block designs are constructed by a D-optimality swapping algorithm that makes improving swaps between 
#'  blocks within the same stratum until a local optima is atttained. The swapping algorithm works from the top stratum downwards and
#'  is always constrained to make improving swaps within the levels of any existing blocks. The whole process is repeated according to the 
#'  number of searches defined by the search parameter and the design returned will be the design with the best overall stratum efficiencies in each stratum 
#'  taken in top-down order.
#'  
#'  Lattice designs where v is a prime-power require the \code{\link[crossdes]{MOLS}} package.
#' 
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A design matrix showing the allocation of treatments to blocks with successive nested blocks factors arranged in successive columns in standard block order.  \cr
#'  \item  A design matrix as above but with the last (bottom) blocks factor shown arranged horizontally to give a plan view. \cr
#'  \item  A set of incidence matrices, one for each blocks stratum, showing the number of times each treatment occurs in each block for each stratum. \cr
#'  \item  A table showing the achieved D- and A-efficiency factors for each nested blocks stratum together with an A-efficiency upper bound, where available. \cr
#'  \item  A table showing a skeleton analysis of degrees of freedom for the combined block and treatment design. \cr
#' } 
#' 
#' Very occasionally, the algorithm may fail to converge due to a near-singular design with a large number of single plot blocks.
#' In that case, it may be best to build a simpler design with larger blocks and then to add the extra block constraints by hand using ad hoc or heuristic methods.     
#' 
#' @param treatments treatment numbers giving a partition of the total required number of treatments into sets of equally replicated treatments.
#' 
#' @param replicates replication numbers giving the replictaion for each set of equally replicated treatments defined by the \code{treatments} partition.
#' 
#' @param blocklevels factor levels that define the number of nested blocks in each succesive blocks stratum taken in order from the highest to the lowest. 
#' The default is the hcf of the replication numbers.
#' 
#' @param seed integer initializing the random number generator. The default is a random seed.
#' 
#' @param searches maximum number of local optima searched for a design optimization. The default is 1 plus the floor of 2000 divided by the number of model parameters.
#' 
#' @param jumps number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#' 
#' @return  
#' \item{Design}{Data frame giving the optimized block and treatment factors in plot order}
#' \item{Plan}{Data frame giving a plan view of the treatments in the bottom stratum of the design}
#' \item{AOV}{Data frame giving a skeleton analysis of variance of the degrees of freedom of the design}
#' \item{Incidences}{Blocks-by-treatments incidence matrices for each stratum of the design}
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
#' # 4 treatments x 4 replicates with 2 main blocks each containing two complete replicates  
#' blocks(treatments=4,replicates=4,blocklevel=2)
#' 
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block 
#' blocks(treatments=50,replicates=4,blocklevels=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments 
#' # giving exactly one single replicate treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),blocklevels=c(4,5))
#' 
#' # 64 treatments x 2 replicates with 2 main blocks and five succesively nested 2-level factors
#' blocks(treatments=64,replicates=2,blocklevels=c(2,2,2,2,2,2))
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,blocklevels=4)
#' 
#' # concurrence matrix of balanced incomplete block design 
#' crossprod(blocks(13,4,13,searches=100)$Incidences[[1]])
#' 
#' # concurrence matrix for 13 treatments x 4 replicates and 13 treatments with one rep in 13 blocks 
#' crossprod(blocks(c(13,13),c(4,1),13)$Incidences[[1]])
#' 
#' # 2**10 treatments x 2 replicates in 2**10 blocks giving a fully saturated blocks design 
#' # (requires a considerable time to run!)
#' \dontrun{ d=blocks(1024,2,rep(2,10)) }
#'          
#' @export
#' @importFrom stats anova lm
#' 
blocks = function( treatments, replicates, blocklevels=HCF(replicates),searches=(1+2000%/%(sum(treatments)+prod(blocklevels))),seed=sample(10000,1),jumps=1) { 
 
# ******************************************************************************************************************************************************** 
#  Generates a vector of block sizes for a particular stratum where all blocks are as equal as possible and never differ by more than a single unit
# ********************************************************************************************************************************************************
  Sizes=function(sizes,blocklevels) { 
    for  (i in 1:length(blocklevels)) {    
      newsizes=NULL
      for (z in 1: length(sizes)) 
        newsizes=c(newsizes, rep(sizes[z] %/% blocklevels[i], blocklevels[i]) + 
                     c( rep(1, sizes[z] %% blocklevels[i]), rep(0,(blocklevels[i]-sizes[z] %% blocklevels[i])))) 
      sizes=newsizes
    }   
    sizes 
  } 
# ******************************************************************************************************************************************************** 
# Finds the highest common factor (hcf) of a set of numbers omitting any zero values (Euclidean algorithm)
# ********************************************************************************************************************************************************
  HCF=function(replevs)  {
    replevs=sort(replevs[replevs>0])
    for (i in 1: length(replevs))
      while (!isTRUE(all.equal(replevs[i]%%replevs[1],0))) replevs[c(1,i)] = c(replevs[i]%%replevs[1], replevs[1])
      replevs[1]
  }   
# ******************************************************************************************************************************************************** 
# Tests a given number for primality and returns TRUE or FALSE
# ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v < 4) return(TRUE) 
    if ( isTRUE(all.equal(v %% 2,0)) ||  isTRUE(all.equal(v %% 3,0)) ) return(FALSE) 
    if (v<25) return(TRUE)
    for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
      if ( isTRUE(all.equal(v %% (i-1) , 0)) ||   isTRUE(all.equal(v %% (i+1) , 0)) ) return(FALSE) 
    return(TRUE)
  } 
# ******************************************************************************************************************************************************** 
# Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only  
# ********************************************************************************************************************************************************
 Contrasts=function(MF,NF) {
    NM=matrix(0,nrow=length(NF),ncol=nlevels(NF))
    NM[cbind(1:length(NF),NF)]=1 # factor indicator matrix  
    for (i in 1:nlevels(MF)) 
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
  optEffics=function(TF,BF) { 
    r=nlevels(TF)
    k=nlevels(BF)
    if (r<=k) 
      e=eigen( (diag(r)-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(r-1)] else    
      e=c(rep(1,(r-k)),
          eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])  
    round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6)
  }
# ******************************************************************************************************************************************************** 
# Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
# Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
# ********************************************************************************************************************************************************
DMax=function(MTT,MBB,MTB,TF,MF,BF) {   
  relD=1
  mainSizes=tabulate(MF)
  nSamp=pmin(rep(8,nlevels(MF)),mainSizes)
  repeat {
    improved=FALSE
    for (k in 1:nlevels(MF)) {
      S=sort(sample((1:length(TF))[MF==k],nSamp[k])) 
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
  Optimise=function(TF,BF,MF,MTT,MBB,MTB,searches,jumps)  {
    globrelD=0
    relD=1
    globTF=TF
    treps=tabulate(TF)
    breps=tabulate(BF)
    if (identical(max(treps),min(treps)) && identical(max(breps),min(breps))  )
      bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) else bound=NA
    for (r in 1 : searches) {
      dmax=DMax(MTT,MBB,MTB,TF,MF,BF) 
      if ( !isTRUE(all.equal(dmax$relD,1)) && dmax$relD>1) {
        relD=relD*dmax$relD
        TF=dmax$TF
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
        if (!isTRUE(all.equal(relD,globrelD)) && relD>globrelD) {
          globTF=TF
          globrelD=relD
          if ( !is.na(bound) && isTRUE(all.equal(bound,optEffics(globTF,BF)[2]))) break
        }
      }
      if (r==searches) break
      for (iswap in 1 : jumps) {
        repeat {     
          s1=sample(1:length(TF),1)
          z=(1:length(TF))[MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1]]
          if (length(z)==0) next
          if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z[1])
          dswap = (1+MTB[TF[s[1]],BF[s[2]]]+MTB[TF[s[2]],BF[s[1]]]-MTB[TF[s[1]],BF[s[1]]]-MTB[TF[s[2]],BF[s[2]]])**2-
            (2*MTT[TF[s[1]],TF[s[2]]]-MTT[TF[s[1]],TF[s[1]]]-MTT[TF[s[2]],TF[s[2]]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]])  
          if (dswap>.1) break
        }
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
  Swaps=function(TF,MF,BF,pivot,rank) {
    n=length(TF)
    candidates=NULL
    while (isTRUE(all.equal(length(candidates),0))) {
      if (rank<(n-1)) s1=sample(pivot[(1+rank):n],1) else s1=pivot[n]
      candidates = rep(1:n)[MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1]]
    }
    if ( length(candidates)>1 )
        s2=sample(candidates,1) else s2=candidates[1] 
    s=c(s1,s2)
  }
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  NonSingular=function(TF,MF,BF) { 
      BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
      BM[cbind(1:length(BF),BF)]=1
      fullrank=nlevels(TF)+nlevels(BF)-1
      TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
      TM[cbind(1:length(TF),TF)]=1
      Q=qr(t(cbind(BM,TM)))
      rank=Q$rank
      pivot=Q$pivot
      times=0
      while (rank<fullrank & times<1000) {
        times=times+1
        s=Swaps(TF,MF,BF,pivot,rank)
        rindex=(1:length(TF))
        rindex[c(s[1],s[2])]=rindex[c(s[2],s[1])]
        newQ=qr(t(cbind(BM,TM[rindex,])))
        if (isTRUE(all.equal(newQ$rank,rank)) || newQ$rank>rank) { 
          TF=TF[rindex]
          TM=TM[rindex,]
          rank=newQ$rank
          pivot=newQ$pivot
        } 
      }
    if (times>999) TF=as.factor(NULL)
    TF
  }  
  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  GenOpt=function(TF,Design,searches,jumps,stratum,blocklevels,hcf) { 
    Design=cbind(as.factor(rep(1,nrow(Design))),Design)
    MF=Design[,stratum]
    BF=Design[,stratum+1]
    rand=sample(1:length(TF))
    TF=TF[rand][order(MF[rand])]
    if (!isTRUE(all.equal(hcf %% prod(blocklevels[1:stratum]),0))) 
    TF=NonSingular(TF,MF,BF)
    if (length(TF)==0) return(TF)
    blevels=nlevels(BF)%/%nlevels(MF)
    BM=Contrasts(MF,BF)[, rep(c(rep(TRUE,(blevels-1)),FALSE),nlevels(MF)),drop=FALSE]
    TM=Contrasts(MF,TF)[,-nlevels(TF),drop=FALSE] 
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    MBB=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
    MTT=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))  
    MTB=matrix(0,nrow=nlevels(TF),ncol=nlevels(BF))
    MBB[1:(nlevels(BF)-nlevels(MF)), 1:(nlevels(BF)-nlevels(MF))]=V[1:(nlevels(BF)-nlevels(MF)),1:(nlevels(BF)-nlevels(MF)),drop=FALSE]
    MTT[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[(nlevels(BF)-nlevels(MF)+1):ncol(V),(nlevels(BF)-nlevels(MF)+1):ncol(V), drop=FALSE]
    MTB[1:(nlevels(TF)-1),  1:(nlevels(BF)-nlevels(MF)) ]=V[(nlevels(BF)-nlevels(MF)+1):ncol(V),1:(nlevels(BF)-nlevels(MF)),drop=FALSE]
    perm=order(order( (1:nlevels(BF))%%blevels ==0  ))  
    MTB=MTB[,perm]
    MBB=MBB[perm,perm] 
    TF=Optimise(TF,BF,MF,MTT,MBB,MTB,searches,jumps)
    TF
  }  
# ******************************************************************************************************************************************************** 
# Generates an initial orthogonal design then builds algebraic lattice blocks or calls the general block design algorithm as appropriate
# ********************************************************************************************************************************************************     
    optTF=function(Design,treatlevs,replevs,blocklevels,searches,jumps) {
    nunits=sum(treatlevs*replevs)
    ntrts=sum(treatlevs)
    hcf=HCF(replevs)
    MB=rep(1:hcf,each=(nunits%/%hcf))
      repeat{
      TF=rep( rep(1:ntrts,rep(replevs%/%hcf,treatlevs)), hcf)
      rand=sample(nunits)
      TF=as.factor( TF[rand][order(MB[rand])] )
      firstPass=TRUE
      v=sqrt(ntrts)  
      regrep=isTRUE(all.equal(max(replevs),min(replevs)))
      for (stratum in 1 : length(blocklevels)) { 
        nblocks=prod(blocklevels[1:stratum])
        if (isTRUE(all.equal(hcf%%nblocks,0))) next
        regular=firstPass &&  regrep && isTRUE(all.equal(nunits%%nblocks,0))
        firstPass=FALSE
        sqrLattice  = regular && isTRUE(all.equal(v,floor(v))) && isTRUE(all.equal(nunits,v*nblocks))
        w=sqrt(nblocks)
        s=ntrts/w
        rectLattice = regular && isTRUE(all.equal(replevs[1],w)) && isTRUE(all.equal(s,floor(s))) && (s<w)
        if (sqrLattice  && replevs[1]<4) {
            t=c(rep(0:(v-1),each=v),rep(0:(v-1),v)+v)
            if (replevs[1]>2)
            for (i in 0: (v-1)) 
              t=c(t,(rep(0:(v-1))+i)%%v + 2*v)
            TF=as.factor(rep(1:(v*v),replevs[1])[order(t)])
        } else if (sqrLattice  &&  replevs[1]<(v+2)  && isPrime(v) ) {
            t=c(rep(0:(v-1),each=v),rep(0:(v-1),v)+v)
            for (z in 1: (replevs[1]-2))
              for (j in 0: (v-1)) 
                t=c(t,(rep(0:(v-1))*z +j)%%v + v*(z+1) )
              TF=as.factor(rep(1:(v*v),replevs[1])[order(t)])
        } else if (sqrLattice  && replevs[1]<(v+2)  &&  ntrts%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
              index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==ntrts)
              mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])			
              TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])
              for (i in 1: (replevs[1]-2))
                TF=c(TF, rep(1:ntrts)[order(as.numeric(mols[,,i]) ) ]) 
              TF=as.factor(TF)
        } else if (sqrLattice  && v==10  && replevs[1]==4) {
          TF=as.factor(c(
            rep(1:100), rep(1:100)[order(rep(0:9,10))],
            rep(1:100)[order(c(
              1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
              5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8))],
            rep(1:100)[order(c(
              1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
              6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6))]
          )) 
        } else if (rectLattice  && isPrime(w) ) {
              for (z in 1:(nunits/nblocks))
                for (j in 0: (w-1)) 
                  for (k in 0: (w-1)) 
                    TF[ (k + j*w)*nunits/nblocks + z]=(j+k*z)%%w + (z-1)*w +1
                  TF=as.factor(TF)
        } else if (rectLattice &&  nblocks%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
              index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==nblocks)
              mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])		
              TF=NULL
              for (i in 1: w)
                for (j in 1: w)
                  for (k in 1: (nunits/nblocks))
                    TF=c(TF,mols[i,j,k]+(k-1)*w)
              TF=as.factor(TF)
        } else TF=GenOpt(TF,Design,searches,jumps,stratum,blocklevels,hcf)
      }
      if (length(TF)>0) break
    }
  levels(TF)=sample(1:ntrts) 
  TF
  }
# ******************************************************************************************************************************************************** 
# Finds efficiency factors for the blocks in each stratum of a design 
# ********************************************************************************************************************************************************     
  A_Efficiencies=function(Design,treatments,replicates)  {
    hcf=HCF(replicates)
    strata=ncol(Design)-2
    nunits=nrow(Design)
    nblocks=as.numeric(sapply(Design,nlevels)[1:strata])
    effics=matrix(1,nrow=strata,ncol=2)
    bounds=rep(NA,strata) 
    for (i in 1:strata) { 
      if (isTRUE(all.equal(max(replicates),min(replicates))) ) {
        if ( isTRUE(all.equal(nunits%%nblocks[i],0))) 
          bounds[i]=upper_bounds(nunits,sum(treatments),nblocks[i]) else if (hcf%%nblocks[i]==0) 
          bounds[i]=1
      }
      if ( sum(treatments)>1 && nblocks[i]>1)
        effics[i,]=optEffics(Design$Treatments,Design[,i])  
    }
    efficiencies=data.frame(cbind( names(Design)[1:strata] ,nblocks, effics, bounds))  
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies[,'Blocks'] = as.factor(efficiencies[,'Blocks'])
    efficiencies
  }
# ******************************************************************************************************************************************************** 
# Carries out some input validation
# ********************************************************************************************************************************************************     
 testInputs=function(treatments,replicates,blocklevels,searches,seed,jumps) {  
   if (missing(treatments) | missing(replicates) )  
     return(" Treatments or replicates not defined ")   
   if (is.null(treatments) | is.null(replicates))  
     return(" Treatments or replicates list is empty ")   
   if (anyNA(treatments) | anyNA(replicates) ) 
     return(" NA values not allowed")
   if (any(!is.finite(treatments)) | any(!is.finite(replicates)) | any(is.nan(treatments)) | any(is.nan(replicates))) 
     return(" Treatments and replicates can contain only finite integers ")
   if ( length(treatments)!=length(replicates) ) 
     return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))

  if (any(treatments<1)) 
     return("Treatments must be non-negative integers")
  if (any(replicates<1)) 
    return("Replicates must be non-negative integers")  
  
   if (!is.null(blocklevels)) {
     if (anyNA(blocklevels) ) return(" NA blocklevels values not allowed") 
     if (any(!is.finite(blocklevels)) | any(is.nan(blocklevels)) ) return(" Blocklevels can contain only finite integers ")
     if (min(blocklevels)<1) return (" Blocklevels must be at least one ")
   }
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
   
   if ( isTRUE( sum(treatments*replicates) < (prod(blocklevels) + sum(treatments)-1) ) )
     return(paste("The total number of plots is",  sum(treatments*replicates) , 
                  "whereas the total required number of model parameters is", prod(blocklevels) + sum(treatments),", which is not feasible. "))  
   return(TRUE)
 }
 
 # ******************************************************************************************************************************************************** 
 # Main body of blocks design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
 # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
 # ********************************************************************************************************************************************************     
 testout=testInputs(treatments,replicates,blocklevels,searches,seed,jumps) 
 if (!isTRUE(testout)) stop(testout)
 set.seed(seed)
 blocklevels=blocklevels[blocklevels>1]
 if (isTRUE(all.equal(length(blocklevels),0))) blocklevels=1
 strata=length(blocklevels)
 cumblocklevs=cumprod(blocklevels)
 if (max(replicates)==1) blocksizes=sum(treatments) else
 blocksizes=Sizes(sum(treatments[replicates>1]*replicates[replicates>1])  , blocklevels)
 stratumnames="Main" 
 if (!isTRUE(all.equal(strata,1)))  
   stratumnames=c( stratumnames,sapply(2:strata, function(i) { paste0("Sub_",(i-1)) }))
 facMat= matrix(nrow=prod(blocklevels),ncol=strata)
 factlevs <- function(r){ gl(cumblocklevs[r],prod(blocklevels)/cumblocklevs[r]) }
 facMat=do.call(cbind,lapply(1:strata,factlevs))
 Design=data.frame(facMat[rep(1:length(blocksizes),blocksizes),])
 Design[]=lapply(Design, as.factor) 

 if (!isTRUE(all.equal(max(replicates),1)) && !isTRUE(all.equal(sum(treatments[replicates>1]),1)) )
      TF=optTF(Design,treatments[replicates>1],replicates[replicates>1],blocklevels,searches,jumps) else
       TF=sample(rep(treatments[replicates>1],replicates[replicates>1])) 
    
 #add back single rep treatments
 if ( min(replicates)==1 && max(replicates)>1 ) {
   addTF=((sum(treatments[replicates>1])+1) :sum(treatments))
   if (length(addTF)>1) addTF=sample(addTF)
   TF=as.factor(  c(TF, addTF )  ) 
   reptrts=NULL
   for (i in 1 : length(replicates))
     if (replicates[i]>1)
       reptrts=c(reptrts,rep(FALSE,treatments[i])) else 
         reptrts=c(reptrts,rep(TRUE,treatments[i]))
   levels(TF)= (1:sum(treatments*replicates))[order(reptrts)] 
   TF=as.numeric(levels(TF))[TF]
   newblocksizes=Sizes(sum(treatments*replicates),blocklevels)
   Design=data.frame(facMat[rep(1:length(newblocksizes),newblocksizes),])
   TF=TF[order(c(rep(1:length(blocksizes),blocksizes),rep(1:length(blocksizes),(newblocksizes-blocksizes))))]
   blocksizes=newblocksizes
 }
  Design[,c("Plots","Treatments")]  = c(rep(1:nrow(Design)) ,TF)
  Design[]=lapply(Design, as.factor)
  # randomization
  Design=data.frame(do.call(cbind,lapply(1:(strata+1), function(r){ sample(nlevels(Design[,r]))[Design[,r]] })) ,Design[,strata+2])
  Design=Design[ do.call(order, Design), ]
  # blocksizes for re-labelled block factor levels
  t=tabulate(Design[,strata])
  blocksizes=t[unique(Design[,strata])]
  # re-labelled block factor levels for re-ordered blocks
  Design[,c(1:(ncol(Design)-1))]=cbind(facMat[rep(1:length(blocksizes),blocksizes),],rep(1:nrow(Design)))
  colnames(Design)=c(stratumnames,"Plots","Treatments")  
  row.names(Design)=c(1:nrow(Design))
  Design[,ncol(Design)-1]=unlist(lapply(1:length(blocksizes),function(r){rep(1:blocksizes[r])}))
  Design[]=lapply(Design, as.factor)
  # incidences
  Incidences=vector(mode = "list", length =strata )
  for (i in 1:strata)
    Incidences[[i]]=table( Design[,i] ,Design[,strata+2])  
  names(Incidences)=stratumnames
  # efficiencies
  Efficiencies=A_Efficiencies(Design,treatments,replicates)
  # aov
    DF=blocklevels-1
    if (strata>1) 
      DF[2:strata]=DF[2:strata]*cumblocklevs[1:(strata-1)]
    DF=c(DF,sum(treatments)-1)
    DF=c(DF,nrow(Design)-sum(DF)-1)
    AOV=data.frame(DF)
    colnames(AOV)="DF"
    rownames(AOV)=c(stratumnames,"Treatments","Residual")
  # treatment replications
  Treatments=data.frame(table(Design[,"Treatments"]))
  Treatments[]=lapply(Treatments, as.factor) 
  colnames(Treatments)=c("Treatments","Replicates")
  # blocksizes
  BlockSizes=data.frame(Design[Design["Plots"]==1,1:strata,drop=FALSE] ,blocksizes)
  BlockSizes[]=lapply(BlockSizes, as.factor) 
  row.names(BlockSizes)=c(1:nrow(BlockSizes))
  colnames(BlockSizes)=c(stratumnames," Sizes ")  
  # convert block levels to nested levels
  Design[,c(1:strata)] = do.call(cbind,lapply(1:strata, function(r){ (as.numeric(Design[,r])-1)%%blocklevels[r]+1 }))
  Design[]=lapply(Design, as.factor)
  # plan
  Plots=t(sapply(split(Design[,"Treatments"],rep(1:length(blocksizes),blocksizes)), '[', 1:max(blocksizes) ))
  Plots[is.na(Plots)] = ""
  Plan=data.frame(Design[Design["Plots"]==1,rep(1:strata),drop=FALSE] ,rep("",nrow(Plots)) ,Plots)
  colnames(Plan)=c(stratumnames,"Plots:",rep(1:max(blocksizes)))
  Plan[]=lapply(Plan,as.factor) 
  row.names(Plan)=c(1:nrow(Plan))
 list(Treatments=Treatments,BlockSizes=BlockSizes,Efficiencies=Efficiencies,Design=Design,Plan=Plan,AOV=AOV,Incidences=Incidences,Seed=seed,Searches=searches,Jumps=jumps) 
} 