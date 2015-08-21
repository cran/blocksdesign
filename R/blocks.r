#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested block designs for unstructured treatments where treatments can have arbitrary levels of replication
#' and blocks can have arbitrary depth of nesting.
#' 
#' @details
#' 
#' The \code{treatments} and \code{replicates} parameters partition the treatments design into sets of equally replicated treatments where 
#' \code{treatments} contains the set sizes and \code{replicates} contains the set replication numbers. 
#'  The sum of the set sizes is the total number of treatments and the sum of the cross-products of the set sizes and the replication numbers
#'  is the total number of plots. Treatments are numbered consecutively according to the treatments in the consecutive treatment sets. 
#' 
#' The \code{blocklevels} parameter contains the numbers of nested blocks in the individual strata taken in order from highest to lowest. The
#' first number is the number of main blocks, the second number, if any, is the number of sub-blocks nested in each main block, the third
#' number, if any, is the number of sub-sub-blocks nested in each sub-block,and so on for all the reqired strata.
#' The default number of blocks is the maximum possible number of orthogonal main blocks (the highest common factor of the replication numbers). 
#' Block sizes are always as equal as possible and never differ by more than a single unit in any particular stratum of the design. 
#' 
#'  Block designs for k replicates of v**2 treatments with k complete main blocks of size v**2 and v incomplete blocks of size v in each main block 
#'  where k < (v+2) if v is prime or prime-power, k < 5 if v = 10, or k < 4 generally, are lattice block designs and are constructed algebraically. 
#'  All other non-orthogonal block designs are constructed by a D-optimality swapping algorithm that makes improving swaps between 
#'  blocks until a local optima is atttained. The swapping algorithm works from the top stratum downwards and
#'  is always constrained to make improving swaps within the levels of any existing blocks. The whole process will be repeated according to the 
#'  number of searches defined by the search parameter and the design returned will be the design with the best overall stratum efficiencies in top-down 
#'  order.
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
#' In that case, it may be best to build a simpler design with larger blocks and then to add the extra block constraints by hand using ad hoc or heuristic 
#' methods.     
#' 
#' @param treatments numbers that provide a partition of the total required number of treatments into sets of equally replicated treatments.
#' 
#' @param replicates replication numbers for each equally replicated treatments set defined by the \code{treatments} list.
#' 
#' @param blocklevels factor levels that define the number of nested blocks in each succesive blocks stratum taken in order from the highest to the lowest. 
#' The default is the hcf of the replication numbers.
#' 
#' @param seed integer initializing the random number generator. The default is a random seed.
#' 
#' @param searches maximum number of local optima searched for a design optimization. The default is the maximum of 1 or (100 - total model terms). 
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
#' 
#' @importFrom stats anova lm
#' 
blocks = function(treatments, replicates, blocklevels=HCF(replicates), searches=max(1,100-sum(treatments)-prod(blocklevels)),seed=sample(10000,1),jumps=1) { 
  
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
  # Finds the highest common factor (hcf) of a set of numbers omitting any zero values
  # ********************************************************************************************************************************************************
  HCF=function(replevs)  {
    replevs=replevs[replevs>0]
    if (length(replevs)==0) return(1)
    replevs=sort(replevs)
    v=(c(replevs[1],NULL))
    if (length(replevs)>1) 
      for (i in 2: length(replevs)) {
        v[2]=replevs[i] 
        while (v[2]%%v[1] != 0) v = c(v[2]%%v[1], v[1]) 
      }
    v[1]
  }   
# ******************************************************************************************************************************************************** 
# Tests a given number for primality and returns TRUE or FALSE
# ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v <= 3)  { 
      return(TRUE)
    } else if (v %% 2 == 0 | v %% 3 == 0) {
      return(FALSE) 
    } else if (v<25) {
      return(TRUE)
    } else {
      for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
        if( v %% (i-1) == 0 | v %% (i+1) == 0) return(FALSE) 
    }
    return(TRUE)
  }      
# ******************************************************************************************************************************************************** 
# Contrasts for factor NF centered within the levels of factor MF thus ensuring that NF information is estimated within the levels of factor MF only  
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
# ********************************************************************************************************************************************************
  UpDate=function(MTT,MBB,MTB,ti,tj,bi,bj) {  
    mtt=MTT[ti,ti]+MTT[tj,tj]-2*MTT[tj,ti]
    mbb=MBB[bi,bi]+MBB[bj,bj]-2*MBB[bi,bj]
    mtb=MTB[ti,bi]-MTB[tj,bi]-MTB[ti,bj]+MTB[tj,bj] 
    f = sqrt(2+mtt+mbb-2*mtb)
    m = f/sqrt(1-2*mtb-mtt*mbb+mtb*mtb)/2 
    Z1 = (MTB[,bi]-MTB[,bj]-MTT[,ti]+MTT[,tj])/f     
    Z2 = (MBB[bi,]-MBB[bj,]-MTB[ti,]+MTB[tj,])/f 
    W1 = (MTT[,ti]-MTT[,tj]+MTB[,bi]-MTB[,bj] - Z1*(mbb-mtt)/f)*m
    W2 = (MTB[ti,]-MTB[tj,]+MBB[bi,]-MBB[bj,] - Z2*(mbb-mtt)/f)*m
    MTT = MTT - tcrossprod(Z1) + tcrossprod(W1)
    MBB = MBB - tcrossprod(Z2) + tcrossprod(W2)
    MTB = MTB - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }   
# ******************************************************************************************************************************************************** 
# Calculates A-efficiency for treatment factor TF assuming block factor BF
# ********************************************************************************************************************************************************
  optEffics=function(TF,BF) { 
    ntrts=nlevels(TF)
    nblks=nlevels(BF)
    if (ntrts<=nblks) {
      e=eigen( (diag(ntrts)-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)]     
   } else {     
      e=c(rep(1,(ntrts-nblks)),
          eigen((diag(nblks)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(nblks-1)])  
   }
    c(round(exp(sum(log(e))/(ntrts-1)),6),round(1/mean(1/e),6))
  }
# ******************************************************************************************************************************************************** 
# Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
# Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
# ********************************************************************************************************************************************************
D_Max=function(MTT,MBB,MTB,TF,MF,BF) {   
  relD=1
  mainSizes=tabulate(MF)
  nSamp=pmin(rep(8,nlevels(MF)),mainSizes)
  repeat {
    improved=FALSE
    for (k in 1:nlevels(MF)) {
      S=sort(sample( (1:length(TF))[MF==k]  ,nSamp[k])) 
      TT=2*MTT[TF[S],TF[S],drop=FALSE]-tcrossprod(MTT[cbind(TF[S],TF[S])]+rep(1,nSamp[k]) ) + tcrossprod(MTT[cbind(TF[S],TF[S])]) + 1
      BB=2*MBB[BF[S],BF[S],drop=FALSE]-tcrossprod(MBB[cbind(BF[S],BF[S])]+rep(1,nSamp[k]) ) + tcrossprod(MBB[cbind(BF[S],BF[S])]) + 1
      TB=MTB[TF[S],BF[S],drop=FALSE]-tcrossprod(MTB[cbind(TF[S],BF[S])],rep(1,nSamp[k]))
      dMat=(TB+t(TB)+1)**2-TT*BB
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
    if (sum(nSamp) < min(length(TF),512))
      nSamp=pmin(mainSizes,2*nSamp)
     else 
      break
  }  
  list(MTT=MTT,MBB=MBB,MTB=MTB,TF=TF,relD=relD)
}  

# ******************************************************************************************************************************************************** 
#  Number of searches for an optimization with selected number of searches and selected number of junps to escape from local optima
# ********************************************************************************************************************************************************
  Optimise=function(TF,BF,MF,MTT,MBB,MTB,searches,jumps)  {
    newrelD=1
    globrelD=0
    globTF=TF
    treps=tabulate(TF)
    breps=tabulate(BF)
    if (identical(max(treps),min(treps)) && identical(max(breps),min(breps))  )
        bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) 
    else
      bound=NA
    for (r in 1 : searches) {
      dmax=D_Max(MTT,MBB,MTB,TF,MF,BF) 
      
      if (!isTRUE(all.equal(dmax$relD,1)) &&  dmax$relD>1) {
        newrelD=newrelD*dmax$relD 
        TF=dmax$TF
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
      
        if (!isTRUE(all.equal(newrelD,globrelD)) &&  newrelD>globrelD) {
          globTF=TF
          globrelD=newrelD
          if ( !is.na(bound) &&  isTRUE( all.equal(bound,  optEffics(globTF,BF)[2])) ) break
        }
      }
      
      if (r==searches) break
      for (iswap in 1 : jumps) {
        dswap=0
        while(isTRUE(all.equal(dswap,0)) | dswap<0) {
          s1=sample(1:length(TF),1)
          z=(1:length(TF))[MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1]]
          if (length(z)==0) next
           s=c(s1,sample(z,1)) 
          dswap = (1+MTB[TF[s[1]],BF[s[2]]]+MTB[TF[s[2]],BF[s[1]]]-MTB[TF[s[1]],BF[s[1]]]-MTB[TF[s[2]],BF[s[2]]])**2-
            (2*MTT[TF[s[1]],TF[s[2]]]-MTT[TF[s[1]],TF[s[1]]]-MTT[TF[s[2]],TF[s[2]]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]])  
        }
        
        newrelD=newrelD*dswap
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
     repeat { 
      s1=sample(pivot[(1+rank):length(pivot)],1)
      s2=(1:length(TF))[MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1]]
      if (length(s2)>0) break
     }
    s=c(s1,s2)
  }
  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  NonSingular=function(TF,MF,BF,cycles) { 
      BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
      BM[cbind(1:length(BF),BF)]=1
      fullrank=nlevels(TF)+nlevels(BF)-1
      TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
      TM[cbind(1:length(TF),TF)]=1
      Q=qr(t(cbind(BM,TM)))
      rank=Q$rank
      pivot=Q$pivot
      searches=0
      while (rank<fullrank & searches<(cycles*5)) {
        searches=searches+1
        s=Swaps(TF,MF,BF,pivot,rank)
        TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]
        newQ=qr(t(cbind(BM,TM)))
          if ( newQ$rank>=rank) { 
            TF[c(s[1],s[2])]=TF[c(s[2],s[1])]
            rank=newQ$rank
            pivot=newQ$pivot
          } else 
            TM[c(s[2],s[1]),]=TM[c(s[1],s[2]),]
        }
      if (searches>=(cycles*5)) TF=NULL 
      TF
  }  
  
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************    
  GenOpt=function(TF,Design,searches,jumps,stratum,blocklevels,hcf,cycles) { 
    MF=Design[,stratum]
    BF=Design[,stratum+1]
    rand=sample(1:length(TF))
    TF=TF[rand][order(MF[rand])]
    if ( !identical( hcf %% prod(blocklevels[1:stratum]), 0)) 
    TF=NonSingular(TF,MF,BF,cycles)
    if (!is.null(TF)) {
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
    }
    TF
  }  
 
# ******************************************************************************************************************************************************** 
# Generates an initial orthogonal design then builds algebraic lattice blocks or calls the general block design algorithm GenOpt as appropriate
# ********************************************************************************************************************************************************     
    optTF=function(Design,treatlevs,replevs,blocklevels,searches,jumps,cycles) {
    nunits=sum(treatlevs*replevs)
    ntrts=sum(treatlevs)
    hcf=HCF(replevs)
    TF=rep( rep(1:ntrts,rep(replevs%/%hcf,treatlevs)), hcf)
    rand=sample(nunits)
    TF=as.factor(TF[rand][order(rep(1:hcf,each=(nunits%/%hcf))[rand])])
      firstNest=TRUE
      for (i in 1 : length(blocklevels)) { 
       if ( identical( hcf %% prod(blocklevels[1:i]), 0)) next
        v=sqrt(ntrts)
        nblocks=prod(blocklevels[1:i])
        reglat=(identical(max(replevs),min(replevs)) && firstNest && identical(nunits/nblocks , v ) && identical(v, floor(v)) && identical(nunits %% nblocks, 0)  ) 
        firstNest=FALSE
        if (  reglat &&  replevs[1]<4   ) {
          TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])
          if (replevs[1]>2) {
            set=NULL
            for (j in 0: (v-1)) 
              for (k in 0: (v-1)) 
                set=c(set, (j+k)%%v )
            TF=c(TF, rep(1:ntrts)[order(set)])
            TF=as.factor(TF)
            levels(TF)=sample(1:ntrts)
          }
        } else if ( reglat &&  replevs[1]<(v+2)  && isPrime(v) ) { 
          TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])
          for (z in 1: (replevs[1]-2)) {
            set=NULL
            for (j in 0: (v-1)) 
              for (k in 0: (v-1)) 
                set=c(set,(j+k*z)%%v)
            TF=c(TF, rep(1:ntrts)[order(set)])
            TF=as.factor(TF)
          levels(TF)=sample(1:ntrts)
          }
        } else if (reglat  &&  replevs[1]<(v+2)  &&  ntrts%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
          index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==ntrts)
          mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])			
          TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])
          for (i in 1: (replevs[1]-2))
            TF=c(TF, rep(1:ntrts)[order(    as.numeric(mols[,,i]) ) ]) 
          TF=as.factor(TF)
          levels(TF)=sample(1:ntrts)
        } else if (  reglat && v==10  && replevs[1]<5  ) {
          TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])  
          if (replevs[1]>2)
            TF=c(TF, rep(1:ntrts)[order(c(
            1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
            5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8))]) 
          if (replevs[1]==4) 
          TF=c(TF, rep(1:ntrts)[order(c(
            1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
            6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6))]) 
          TF=as.factor(TF)
          levels(TF)=sample(1:ntrts)
        } else {
          TF=GenOpt(TF,Design,searches,jumps,i,blocklevels,hcf,cycles)
        }
        if (is.null(TF)) break
      }   
   TF 
    }
    
# ******************************************************************************************************************************************************** 
# Single replicate treatments are omitted before blocking and are then replaced after optimization. Single replicate treatments are added back one block 
# at a time while maintaining block size equality or near-equality with not more than one unit difference in each block for each stratum of the design.
# ********************************************************************************************************************************************************     
  fullDesign=function(Design,facMat,treatments,replicates,oldblocksizes,blocklevels) {
    strata=ncol(Design)-3
    TF=Design[,ncol(Design)]
    nunits=sum(treatments*replicates)
    ntrts=sum(treatments)
    TF=as.factor(c( as.numeric(TF), ((nlevels(TF)+1):ntrts)[sample(ntrts-nlevels(TF))]  ) ) 
    trtlabs=NULL  
    extlabs=NULL
    index=0
    for (i in 1 : length(treatments)) {
      if (replicates[i]>1) {
        trtlabs=c(trtlabs,  (index+1):(index+treatments[i]) )
     } else {
       extlabs=c(extlabs,  (index+1):(index+treatments[i]) )
     }
        index=index+treatments[i] 
    }    
    levels(TF)=c(trtlabs,extlabs)
    TF=as.numeric(levels(TF))[TF]
    newblocksizes=Sizes(nunits,blocklevels)
    BF=c( rep( 1:length(oldblocksizes),oldblocksizes),  rep( 1:length(oldblocksizes),(newblocksizes-oldblocksizes) ) )
    Design=facMat[rep(1:length(newblocksizes),newblocksizes),]
    Design=as.data.frame(cbind(rep(1,nunits), Design, rep(1:nunits),TF[order(BF)]))
    Design[]=lapply(Design, factor) 
    Design
  } 
  
# ******************************************************************************************************************************************************** 
# Finds A-efficiency factors for the blocks in each stratum of a design 
# ********************************************************************************************************************************************************     
  A_Efficiencies=function(Design)  {
    strata=ncol(Design)-2
    treps=tabulate(Design$Treatments)
    effics=matrix(1,nrow=strata,ncol=2)
    bounds=rep(NA,strata) 
    blocks=rep(0,strata)  
    for (i in 1:strata) { 
      blocks[i]=nlevels(Design[,i])   
    breps=tabulate(Design[,i])
      
    if ( all(treps==treps[1]) && all(breps==breps[1]) )
     bounds[i]=upper_bounds(nrow(Design),nlevels(Design$Treatments),blocks[i])  
      
    if (nlevels(Design$Treatments)>1 && nlevels(Design[,i])>1)
      effics[i,]=optEffics(Design$Treatments,Design[,i])  
    }
    efficiencies=as.data.frame(cbind(names(Design)[1:strata]    ,blocks, effics, bounds))  
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies[, 'Blocks'] = as.factor(efficiencies[, 'Blocks'])
    efficiencies
  }
# ******************************************************************************************************************************************************** 
# Randomizes nested blocks for each nested stratum of a design. Nested blocks are randomised within the blocks of the next higher blocks stratum.
# Plots are randomised within the bottom stratum of the design.
# ********************************************************************************************************************************************************     
  randBlocks=function(Design,facMat) {
    for (r in 1 : (ncol(Design)-1) )
      Design[,r]=as.numeric( sample(nlevels( Design[,r]) ))[Design[,r]]  
    Design=Design[ do.call(order, Design), ] 
    blocksizes=tabulate(as.numeric( order(unique(Design[,ncol(Design)-2])))[Design[,ncol(Design)-2]])
    Design[,1 : (ncol(Design)-2)] = facMat[rep(1:length(blocksizes),blocksizes),1 : (ncol(Design)-2)]
    Design[,(ncol(Design)-1)]=rep(1:nrow(Design))
    Design[]=lapply(Design, factor) 
    Design
  }
# ******************************************************************************************************************************************************** 
#  Plan output for printing
# ********************************************************************************************************************************************************     
  Plan=function(Design)  {
    nblocks=nlevels(Design[,strata])
    trts=as.numeric(levels(Design[,ncol(Design)]))[Design[,ncol(Design)]]
    strata=ncol(Design)-2
    bSizes=c(0,tabulate(Design[,strata]))
    cumSizes=bSizes
    for (i in 1:nblocks)
      cumSizes[i+1]=cumSizes[i]+cumSizes[i+1]    
    for (i in 1:nblocks)
      Design[i,]=Design[1+sum(bSizes[1:i]) ,]
    Design=Design[1:nblocks,1:strata, drop = FALSE]
    fullsize=max(bSizes)
    fulltrts=rep(NA, fullsize*nblocks)
    for (i in 1:nblocks) 
      fulltrts[ ((i-1)*fullsize+1) : ((i-1)*fullsize+bSizes[i+1]) ] = trts[  (cumSizes[i]+1 ) : cumSizes[i+1]   ]   
    plan=matrix(fulltrts,nrow=nblocks,ncol=max(bSizes),byrow=TRUE)
    stratumnames=colnames(Design)
    Design=cbind(Design,rep(" ",nblocks),plan)
    Design[is.na(Design)]  = " "   
    colnames(Design)=c(stratumnames , "Sub_plots", 1:ncol(plan) )
    Design
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
   if (!all(is.finite(treatments)) | !all(is.finite(replicates)) | !all(!is.nan(treatments)) | !all(!is.nan(replicates))) 
     return(" Treatments and replicates can contain only finite integers ")
   if ( length(treatments)!=length(replicates) ) 
     return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
  if (!all(treatments>=0)) 
     return("Treatments must be non-negative integers")
  if (!all(replicates>=0)) 
    return("Replicates must be non-negative integers")  
  if (  sum(treatments*replicates) <=0 ) 
    return("Design cannot be fitted : number of plots must be greater than zero")  
  
   if (!is.null(blocklevels)) {
     if (anyNA(blocklevels) ) return(" NA blocklevels values not allowed") 
     if (!all(is.finite(blocklevels)) | !all(!is.nan(blocklevels)) ) return(" Blocklevels can contain only finite integers ")
     if (min(blocklevels)<1) return (" Blocklevels must be at least one ")
   }
   if (!is.null(searches)) {
     if (anyNA(searches) ) return(" NA searches values not allowed") 
     if ( !all(is.finite(searches)) | !all(!is.nan(searches))) return(" Searches must be a finite integer ") 
     if (searches<1)  return(" Repeats must be at least one ")   
   }  
  if (!is.null(jumps)) {
    if (anyNA(jumps) ) return(" NA jumps values not allowed") 
    if ( !all(is.finite(jumps)) | !all(!is.nan(jumps))) return(" jumps must be a finite integer ") 
    if (jumps<1)  return(" Random jumps must be at least one ")   
  }    
   if (!is.null(seed)) {
     if (anyNA(seed) ) return(" NA seed values not allowed") 
     if ( !all(is.finite(seed)) | !all(!is.nan(seed))) return(" Seed must be a finite integer ") 
     if (seed<1)  return(" Seed must be at least one ")   
   } 
   if (  sum(treatments*replicates) < (prod(blocklevels) + sum(treatments)-1) ) 
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
 if (is.null(seed)) seed=sample(1:100000,1)
 set.seed(seed) 
 if (is.null(jumps)) jumps=1
 if (blocklevels[1]>1)
   stratumnames="Main" else
     stratumnames=NULL
 
 if (length(blocklevels)>1)
   for (i in 2:length(blocklevels))
     if (blocklevels[i]>1)
       stratumnames=c(stratumnames,paste("Sub",(i-1), sep="_"))  
 
 sets=treatments*replicates
 treatments=treatments[sets>0]
 replicates=replicates[sets>0]
 
 if (!all(blocklevels==1))
   blocklevels=blocklevels[blocklevels>1]
 else {
   blocklevels=1
   stratumnames="Main"
 }
 if (max(replicates)==1) {
   nunits=sum(treatments) 
   strata=1
   Design=as.data.frame(cbind(rep(1,each=nunits), rep(1:nunits), sample(1:nunits)    ))
   colnames(Design)=c(stratumnames,"Plots","Treatments") 
   Design[]=lapply(Design, factor) 
   
   facMat= matrix(nrow=prod(blocklevels),ncol=strata)
   for (r in 1 : strata) 
     facMat[,r]=gl(prod(blocklevels[1:r]),prod(blocklevels)/prod(blocklevels[1:r])  )  
 } else {  
   treatlevs=treatments[replicates>1]
   replevs = replicates[replicates>1]
   
   nunits=sum(treatlevs*replevs) 
   if (is.null(searches)) 
     searches=1+2000%/%(sum(treatments)+prod(blocklevels))
   
   strata=length(blocklevels)
   blocksizes=Sizes(nunits,blocklevels)
   
   facMat= matrix(nrow=prod(blocklevels),ncol=strata)
   for (r in 1 : strata) 
     facMat[,r]=gl(prod(blocklevels[1:r]),prod(blocklevels)/prod(blocklevels[1:r])  )  
   
   Design=facMat[rep(1:length(blocksizes),blocksizes),]
   Design=as.data.frame(cbind(rep(1,nunits), Design, rep(1:nunits)))
   Design[]=lapply(Design, factor) 
   TF=NULL
   cycles=1
   if (sum(treatlevs)==1) 
     TF=as.factor(rep(1,replevs))
   else
     while (is.null(TF) && cycles<100) {
       TF=optTF(Design,treatlevs,replevs,blocklevels,searches,jumps,cycles)
       cycles=cycles+1
     }
   
   if (cycles>=100) stop("Cannot find a non-singular starting design for every blocks stratum - please try a simpler design structure")  
   Design=cbind(Design,TF)  
   
   # add back single replicate treatments here 
   if (!all(replicates>1) )
     Design= fullDesign(Design,facMat,treatments,replicates,blocksizes,blocklevels) 
   Design=Design[,-1] 
   
   # randomization
   Design=randBlocks(Design,facMat)
   
   colnames(Design)=c(stratumnames,"Plots","Treatments")   
   rownames(Design) = NULL 
 }
 Incidences=vector(mode = "list", length =strata )
 for (i in 1:strata)
   Incidences[[i]]=table( Design[,i] ,Design[,strata+2])  
 names(Incidences)=stratumnames
 BlockSizes=as.data.frame(cbind(facMat,tabulate(Design[,strata])))
 BlockSizes[]=lapply(BlockSizes, factor) 
 colnames(BlockSizes)=c(stratumnames," Sizes ")  
 
 Treatments=as.data.frame(table(Design[,"Treatments"]))
 Treatments[]=lapply(Treatments, factor) 
 colnames(Treatments)=c("Treatments","Replicates")
 
 dfDesign = Design[-(ncol(Design)-1)]
 factlevels=unlist(lapply( dfDesign , nlevels))
 if ( max(factlevels)>1 )
   dfDesign=dfDesign[, unlist(lapply( dfDesign , nlevels))>1, drop=FALSE]  else {
     dfDesign=as.data.frame(cbind( c(1,2)  ))
     colnames(dfDesign)="Treatments"
     dfDesign[]=lapply(dfDesign, factor)  
   }
 Dummy=suppressWarnings(anova(lm(rnorm(nrow(dfDesign)) ~ ., data = dfDesign))) 
 if (max(factlevels)==1) Dummy[1,1]=0
 Dummy=Dummy[,1,drop=FALSE]
 list(Treatments=Treatments,BlockSizes=BlockSizes,Efficiencies=A_Efficiencies(Design),Design=Design,Plan=Plan(Design),AOV=Dummy,Incidences=Incidences,Seed=seed,Searches=searches,Jumps=jumps) 
} 
 
 
 