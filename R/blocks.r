#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested block designs for unstructured treatments with arbitrary replication, 
#' not necessarily all equal, and arbitrary depth of nesting.
#' 
#' @details
#' 
#' The \code{blocks} function optimizes nested blocks designs where treatments can have any arbitrary level of replication, not necessarily all equal,
#' and blocks can be either a simple main blocks design or a nested blocks design with any feasible depth of nesting. 
#' 
#' The \code{treatments} and \code{replicates} arguments define the treatment structure of the design. \code{treatments} is a
#' set of cardinal numbers that define sets of equally replicated treatments and \code{replicates} is a matching set of replication numbers for those sets.
#' The sum of the cardinal numbers is the total number of treatments and the sum of the cross-products of the cardinal numbers and the replication 
#' numbers is the total number of units. Treatments are numbered consecutively according to the ordering of the treatment sets but
#'  different sets with the same replication can be used if arbitrary numbering is required. 
#'  Single replicate treatments sets are permitted provided there are at least some replicated treatment sets in the design.
#' 
#' The \code{blocklevels} argument is a vector of integers that defines the blocks structure of the design. The length of the vector is the total number of 
#' strata while the numeric elements are the block levels of the successive strata. 
#' The first number is the number of main blocks, the second, if any, is the number of sub-blocks nested in each main block, the third, if any,
#' is the number of sub-sub-blocks nested in each sub-block and so on. The default is a main blocks design with the maximum possible number of
#' orthogonal main blocks. Block sizes are always as equal as possible and will never differ by more than a single unit 
#' in any one stratum. 
#' 
#' The \code{searches} argument is the integer number of searches for an optimization. Ideally, the number of searches should be as large 
#'  as the computational resources permit.  
#'  
#' The \code{jumps} argument is the number of random swaps needed to escape a local maxima. A single swap appears to work well and
#' is the most efficient choice for the updating algorithm but the setting can be increased to any integer value, if required.
#' 
#' There are three classes of solution for the designs defined by the above arguments:   
#' \itemize{
#' 
#' \item Complete block designs constructed algebraically by direct substitution of complete treatment sets. Complete blocks can contain multiple 
#' sets of treatments but the equal block size restriction means that the number of blocks must divide the number of treatment replicates (see example).  
#'   
#' \item Lattice block designs with k replicates of v**2 treatments in k main blocks where each main block contains 
#' v incomplete blocks of size v and k < (v+2) for prime or prime-power v, or k < 5 for v = 10 or k < 4 for any other v. Lattice designs are constructed
#' by algebraic methods using Latin squares. The \code{\link[crossdes]{MOLS}} package is used if v is a prime-power.
#' 
#' \item  General block designs with arbitrary depth of nesting. These designs are optimized by making D-efficiency improving swaps between 
#' blocks nested within the constraints of higher level blocks, if any, until a local maxima is attained.
#' For repeated \code{searches}, a local maxima is escaped by one or more random swaps (= \code{jumps}) between blocks, again nested within
#'  the constraints of higher level blocks, if any.
#' The algorithm works from the top stratum downwards and the design of each nested stratum is conditional on 
#' the design of all higher level strata. The efficiency of nested designs may be slightly reduced relative to the efficiency of comparable non-nested designs
#' but, for all practical purposes, this reduction is of no practical consequence.
#'    
#' }
#' 
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A design matrix showing the allocation of treatments to blocks with the successive stratum blocks arranged in successive columns in standard block order.  \cr
#'  \item  A design matrix as above but with the bottom stratum blocks arranged in a plan view. \cr
#'  \item  A set of incidence matrices, one for each stratum, showing the number of times each treatment occurs in each block for each stratum. \cr
#'  \item  A table showing the achieved D- and A-efficiency factors for each nested blocks stratum together with an A-efficiency upper bound, where available. \cr
#' } 
#' 
#' @param treatments a set of cardinal numbers where the sum of the cardinals is the required number of treatments and the individual cardinals 
#' represent individual sets of equally replicated treatments
#' 
#' @param replicates a vector of replication numbers, one for each of the individual treatment sets defined above,  
#' where the replication numbers and the treatment sets are assumed to be in matching order.   
#' 
#' @param blocklevels an optional vector of integers where the first integer is the number of main blocks and the remaining integers, if any,
#' are the numbers of nested sub-blocks in a hierarchy of nested sub-blocks. The default is the HCF of the replication numbers.
#' 
#' @param seed an optional integer for initializing the random number generator. The default is a random integer less than 10000.
#' 
#' @param searches an optional integer for the maximum number of searches during optimization. The default is the maximum of 1 or (100 - sum of model terms). 
#' 
#' @param jumps an optional integer for the number of pairwise random treatment swaps used to escape a local maxima in a stratum. The default is a single swap.
#' 
#' @return  
#' \item{Design}{Data frame showing the optimized block and treatment factors in plot order}
#' \item{Plan}{Data frame showing a plan view of the treatments in the bottom stratum of the design}
#' \item{Incidences}{Blocks-by-treatments incidence matrices in each stratum of the design}
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
#'          
#' @export
#' 
blocks = function(treatments, replicates, blocklevels=HCF(replicates), searches=max(1,100-sum(treatments)-prod(blocklevels)),seed=sample(10000,1),jumps=1) { 

  #******************************************************** sizes ************************************************************************************ 
  Sizes=function(sizes,blocklevels) { 
    for  (i in 1:length(blocklevels)) {    
      newsizes=NULL
      for (z in 1: length(sizes)) 
        newsizes=c(newsizes, rep(sizes[z] %/% blocklevels[i], blocklevels[i]) + c( rep(1, sizes[z] %% blocklevels[i]), rep(0,(blocklevels[i]-sizes[z] %% blocklevels[i])))) 
      sizes=newsizes
    }   
    sizes 
  } 
     
  #******************************************************** Primality test ************************************************************************************
  isPrime=function(v) {
    if (v <= 3)  return(TRUE)
    else if (v %% 2 == 0 | v %% 3 == 0) return(FALSE)
    else if (v<25) return(TRUE)
    else 
      for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
        if( v %% (i-1) == 0 | v %% (i+1) == 0) return(FALSE) 
    return(TRUE)
  }      
  
  #******************************************************** Block contrasts ************************************************************************************
  BlockContrasts=function(MF,BF) {
    BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
    BM[cbind(1:length(BF),as.numeric(BF))]=1 # factor indicator matrix
    BM=BM[,rep(c(rep(TRUE,((nlevels(BF)/nlevels(MF))-1)),FALSE),nlevels(MF)),drop=FALSE]
    for (i in 1: nlevels(MF)) 
      BM[MF==i,]=scale(BM[MF==i,], center = TRUE, scale = FALSE) # scaling within main blocks
    BM
  } 
  
  #******************************************************** Treatment contrasts**********************************************************************************
  TreatContrasts=function(MF,TF) {
    TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
    TM[cbind(1:length(TF),as.numeric(TF))]=1 # factor indicator matrix	
    TM=TM[,c(rep(TRUE,(nlevels(TF)-1)),FALSE),drop=FALSE]	
    for (i in 1:nlevels(MF)) 
      TM[MF==i,]=scale(TM[MF==i,] , center = TRUE, scale = FALSE)
    TM
  }
  
  #******************************************************** Updates variance matrix ************************************************************************************
  UpDate=function(M11,M22,M12,ti,tj,bi,bj,TF,BF) {  
    m11=M11[ti,ti]+M11[tj,tj]-M11[tj,ti]-M11[ti,tj]
    m22=M22[bi,bi]+M22[bj,bj]-M22[bi,bj]-M22[bj,bi]
    m12=M12[ti,bi]-M12[tj,bi]-M12[ti,bj]+M12[tj,bj]   
    f = sqrt(2+m11+m22-2*m12)
    m = f/sqrt(1-2*m12-m11*m22+m12*m12)/2 
    Z1 = (M12[,bi]-M12[,bj]-M11[,ti]+M11[,tj])/f     
    Z2 = (M22[bi,]-M22[bj,]-M12[ti,]+M12[tj,])/f 
    W1 = (M11[,ti]-M11[,tj]+M12[,bi]-M12[,bj] - Z1*(m22-m11)/f)*m
    W2 = (M12[ti,]-M12[tj,]+M22[bi,]-M22[bj,] - Z2*(m22-m11)/f)*m
    M11 = M11 - tcrossprod(Z1) + tcrossprod(W1)
    M22 = M22 - tcrossprod(Z2) + tcrossprod(W2)
    M12 = M12 - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(M11=M11,M22=M22,M12=M12)
  }   
  
  #********************************************************Determinants of jumps using samples of increasing size***********************************************
  D_Max=function(M11,M22,M12,TF,MF,BF) {      
    relD=1
    mainSizes=tabulate(MF)
    nSamp=pmin(rep(8,nlevels(MF)),mainSizes)
    mainBlocks=split(rep(1:length(TF)),MF)  
    repeat {
      improved=FALSE
      for (k in 1:nlevels(MF)) {
        S=sort(sample(mainBlocks[[k]],nSamp[k]))  
        TT=2*M11[TF[S],TF[S],drop=FALSE]-tcrossprod(M11[cbind(TF[S],TF[S])]+rep(1,nSamp[k])) + tcrossprod(M11[cbind(TF[S],TF[S])]) + 1
        BB=2*M22[BF[S],BF[S],drop=FALSE]-tcrossprod(M22[cbind(BF[S],BF[S])]+rep(1,nSamp[k])) + tcrossprod(M22[cbind(BF[S],BF[S])]) + 1
        TB=M12[TF[S],BF[S],drop=FALSE]+t(M12[TF[S],BF[S],drop=FALSE])-tcrossprod(M12[cbind(TF[S],BF[S])]+rep(1,nSamp[k]))+tcrossprod(M12[cbind(TF[S],BF[S])]) + 2
        dMat=TB**2-TT*BB
        sampn=which.max(dMat)   
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        if ( !isTRUE(all.equal(dMat[i,j],1)) & dMat[i,j]>1) {
          improved=TRUE
          relD=relD*dMat[i,j]
          up=UpDate(M11,M22,M12,TF[S[i]],TF[S[j]], BF[S[i]], BF[S[j]], TF,BF)
          M11=up$M11
          M22=up$M22
          M12=up$M12
          TF[c(S[i],S[j])]=TF[c(S[j],S[i])]
        }
      } 
      if (improved) next
      if (sum(nSamp) < min(length(TF),512))
          nSamp=pmin(mainSizes,2*nSamp)
       else 
         break
    }  
    list(M11=M11,M22=M22,M12=M12,TF=TF,relD=relD)
  }
  
  #**************************** Calculates A-optimality *******************************************************
  optEffics=function(TF,BF,ntrts,nblks)   { 
    if (ntrts<=nblks) 
      e=eigen( (diag(ntrts)-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(ntrts-1)]     
    else       
      e=c(rep(1,(ntrts-nblks)),eigen((diag(nblks)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(nblks-1)])    
    aeff =1/mean(1/e) 
    deff = exp(sum(log(e))/(ntrts-1))
    c(deff,aeff)
  }
  
  #**************************** General optimization using multiple searches *******************************************************
  Optimise=function(TF,BF,MF,M11,M22,M12,searches,jumps)  {
    relD=1
    globrelD=0
    globTF=TF
    treps=as.integer(tabulate(TF))
    breps=as.integer(tabulate(BF))
    bound=NA
    if (identical(max(treps),min(treps)) & identical(max(breps),min(breps))  )
        bound=upper_bounds(length(TF),nlevels(TF),nlevels(BF)) 
    for (r in 1 : searches) {
      dmax=D_Max(M11,M22,M12,TF,MF,BF)  
      relD=relD*dmax$relD
      TF=dmax$TF
      M11=dmax$M11
      M22=dmax$M22
      M12=dmax$M12  
      if (!isTRUE(all.equal(relD,globrelD)) &  relD>globrelD) {
        globTF=TF
        globrelD=relD
        if ( !is.na(bound) & r<searches )
          if (isTRUE( all.equal(bound,  optEffics(globTF,BF,nlevels(TF),nlevels(BF))[2]))) break
      }
      if (r==searches) break
      for (iswap in 1 : jumps) {
        dswap=0
        while(isTRUE(all.equal(dswap,0)) | dswap<0) {
          s=sample(rep(1:length(TF))[MF==sample(nlevels(MF),1)],2)
          if ( identical(TF[s[1]],TF[s[2]]) | identical(BF[s[1]],BF[s[2]])  ) next
          dswap = (1+M12[TF[s[1]],BF[s[2]]]+M12[TF[s[2]],BF[s[1]]]-M12[TF[s[1]],BF[s[1]]]-M12[TF[s[2]],BF[s[2]]])**2-
            (2*M11[TF[s[1]],TF[s[2]]]-M11[TF[s[1]],TF[s[1]]]-M11[TF[s[2]],TF[s[2]]])*(2*M22[BF[s[1]],BF[s[2]]]-M22[BF[s[1]],BF[s[1]]]-M22[BF[s[2]],BF[s[2]]])  
        }
        relD=relD*dswap
        up=UpDate(M11,M22,M12,TF[s[1]],TF[s[2]], BF[s[1]], BF[s[2]],TF,BF)
        M11=up$M11
        M22=up$M22
        M12=up$M12
        TF[c(s[1],s[2])]=TF[c(s[2],s[1])]  
      } 
    }
   globTF
  } 
  
  #******************************************************** Initializes design***************************************************************************************
  GenOpt=function(TF,BF,MF,searches,jumps) {   
    TB=TreatContrasts(MF,TF)
    NB=BlockContrasts(MF,BF) 
    DD=crossprod(cbind(TB,NB))
    count=0
      while ( !identical(qr(DD)$rank,ncol(DD)) & count<200) 
      {
        rand=sample(1:length(TF))
        TF=TF[rand][order(MF[rand])] 
        TB=TB[rand,][order(MF[rand]),]
        DD=crossprod(cbind(TB , NB))
        count=count+1
      }
    if (count>199) stop("Unable to find a suitable starting design - perhaps the design is near singular?")
    V=chol2inv(chol(DD))
    M11=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))	
    M22=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
    M12=matrix(0,nrow=nlevels(TF),ncol=nlevels(BF))
    M11[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[1:(nlevels(TF)-1),1:(nlevels(TF)-1),drop=FALSE]
    M12[1:(nlevels(TF)-1),1:(ncol(V)-nlevels(TF)+1)]=V[1:(nlevels(TF)-1),nlevels(TF):ncol(V),drop=FALSE]
    M22[1:(ncol(V)-nlevels(TF)+1),1:(ncol(V)-nlevels(TF)+1)]=V[nlevels(TF):ncol(V),nlevels(TF):ncol(V),drop=FALSE]
    perm=order(order((1:nlevels(BF))%%(nlevels(BF)/nlevels(MF))==0)   )
    M12=M12[,perm]
    M22=M22[perm,perm]	
    TF=Optimise(TF,BF,MF,M11,M22,M12,searches,jumps)
  }
  
  #******************************************************** HCF of replicates************************************************************************************
  HCF=function(replevs)  {
    replevs=sort(replevs)
    v=(c(replevs[1],NULL))
    if (length(replevs)>1) 
      for (i in 2: length(replevs)) {
        v[2]=replevs[i] 
        while (v[2]%%v[1] != 0) v = c(v[2]%%v[1], v[1]) 
      }
    v[1]
  }
 
  #********************************************************  algorithm   ************************************************************************************
    optTF=function(Design,treatlevs,replevs,searches,jumps)  {
    nunits=nrow(Design)
    strata=ncol(Design)-2    
    ntrts=sum(treatlevs)
    hcf=HCF(replevs)
    TF=NULL
    for (i in 1: hcf) {
      s=sample(1:ntrts)
      TF=c(TF, rep(s , rep(replevs/hcf,treatlevs)[s]) )
    }
    TF=as.factor(TF)
    regreps=(identical(max(replevs),min(replevs) ) )
    index=0
    for (i in 1 : strata) { 
      bsizes=tabulate(Design[,i+1])
      if (all( bsizes %% (nunits/hcf)  == 0)) next 
      index=index+1
      v=sqrt(ntrts)
      nblocks=nlevels(Design[,i+1])
      reglat=(regreps & index==1 &  identical( nunits/nblocks , v ) &  identical( v, floor(v)) & identical(max(bsizes),min(bsizes))) 
      if (  reglat &  replevs[1]<4   ) {
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
      } else if ( reglat &  replevs[1]<(v+2)  & isPrime(v) ) { 
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
      } else if (reglat  &  replevs[1]<(v+2)  &  ntrts%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
        index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==ntrts)
        mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])			
        TF=c(rep(1:ntrts), rep(1:ntrts)[order(rep(0:(v-1),v))])
        for (i in 1: (replevs[1]-2))
          TF=c(TF, rep(1:ntrts)[order(    as.numeric(mols[,,i]) ) ]) 
        TF=as.factor(TF)
        levels(TF)=sample(1:ntrts)
      } else if (  reglat & v==10  & replevs[1]<5  ) {
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
        TF=GenOpt(as.factor(TF),Design[,(i+1)],Design[,i],searches,jumps)
      }
    }
   TF 
  }
  
  #******************************************************** replaces single rep treatments *************************************************************************** 
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
      if (replicates[i]>1) 
        trtlabs=c(trtlabs,  (index+1):(index+treatments[i]) )
      else extlabs=c(extlabs,  (index+1):(index+treatments[i]) )
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
  
  #******************************************************** A-efficiencies ************************************************************************************
  A_Efficiencies=function(Design)  {
    strata=ncol(Design)-2
    treps=tabulate(Design$Treatments)
    effics=matrix(1,nrow=strata,ncol=2)
    bounds=rep(NA,strata) 
    blocks=rep(0,strata)  
    for (i in 1:strata) { 
      blocks[i]=nlevels(Design[,i])
      breps=tabulate(Design[,i])
      if ( all(treps==treps[1]) & all(breps==breps[1]))
        bounds[i]=upper_bounds(nrow(Design),nlevels(Design$Treatments),blocks[i])    
      if (nlevels(Design$Treatments)>1 & nlevels(Design[,i])>1)
        effics[i,]=optEffics(Design$Treatments,Design[,i],nlevels(Design$Treatments),blocks[i])  
    }
    efficiencies=as.data.frame(cbind(blocks, effics, bounds))    
    colnames(efficiencies)=c("Blocks","D-Efficiencies","A-Efficiencies", "A-Upper Bounds")
    rnames=c("Main")
    if (strata>1)
      for (i in 1 : (strata-1)) rnames=c(rnames,paste("Sub",i))
    rownames(efficiencies)=rnames 
    efficiencies
  }
   
  #******************************************************** Randomizes blocks within strata************************************************************************************ 
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
 
  #******************************************************** Plan output************************************************************************************
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
  
 #******************************************************** Validates inputs************************************************************************************
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
  if (!all(treatments>=1)) 
     return("Treatments must be integers greater than zero")
  if (!all(replicates>=1)) 
    return("Replicates must be integers greater than zero")  
  if (all(replicates==1)) 
    return("Not all treatments can have only a single replication" )  
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
   if (  sum(treatments*replicates) < (prod(blocklevels) + sum(treatments)) ) 
     return("Design cannot be fitted :  too many blocks and treatments for the available plots")  
   return(TRUE)
 }
 
 #********************************************************design ************************************************************************************ 
 testout=testInputs(treatments,replicates,blocklevels,searches,seed,jumps) 
  if (!isTRUE(testout)) stop(testout)
  if (is.null(seed)) seed=sample(1:100000,1)
  set.seed(seed) 
 if (is.null(jumps)) jumps=1
  # omit any single replicate treatments here 
   treatlevs=treatments[replicates>1]
   replevs = replicates[replicates>1]
  if (is.null(blocklevels)) 
    blocklevels=HCF(replevs)
  nunits=sum(treatlevs*replevs) 
  if (is.null(searches)) 
    searches=1+2000%/%(sum(treatments)+prod(blocklevels))
 if (!all(blocklevels==1))
    blocklevels=blocklevels[blocklevels>1]
 else
   blocklevels=1
  strata=length(blocklevels)
 for (i in 1:strata)
  blocksizes=Sizes(nunits,blocklevels)
  facMat= matrix(nrow=prod(blocklevels),ncol=strata)
  for (r in 1 : strata) 
    facMat[,r]=gl(prod(blocklevels[1:r]),prod(blocklevels)/prod(blocklevels[1:r])  )  
  Design=facMat[rep(1:length(blocksizes),blocksizes),]
  Design=as.data.frame(cbind(rep(1,nunits), Design, rep(1:nunits)))
  Design[]=lapply(Design, factor) 
  Design=cbind(Design,optTF(Design,treatlevs,replevs,searches,jumps))  
  # add back single replicate treatments here 
  if (!all(replicates>1) )
   Design= fullDesign(Design,facMat,treatments,replicates,blocksizes,blocklevels) 
  Design=Design[,c(2:ncol(Design))] 
  # randomization
  Design=randBlocks(Design,facMat)
  stratumnames=c("Main_blocks")
  if (strata>1)
    for (i in 1:(strata-1))
      stratumnames=c(stratumnames,paste("Sub",i,"_blocks", sep=""))  
  colnames(Design)=c(stratumnames,"Sub-plots","Treatments")   
  rownames(Design) = NULL 
  Incidences=vector(mode = "list", length =strata )
  for (i in 1:strata)
    Incidences[[i]]=table( Design[,i] ,Design[,strata+2])  
  names(Incidences)=stratumnames
  list(Design=Design,Plan=Plan(Design),Incidences=Incidences,Efficiencies=A_Efficiencies(Design),Seed=seed,Searches=searches,Jumps=jumps) 
} 