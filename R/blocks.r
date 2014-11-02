#' @title Block designs 
#' 
#' @description
#' 
#' Constructs nested block designs for unstructured treatments with arbitrary replication, not necessarily all equal, and arbitrary depth of nesting.
#' 
#' @details
#' 
#' The \code{blocks} function constructs optimized nested block designs for unstructured treatment sets where treatments can have any arbitrary replication, not necessarily all equal, 
#' and blocks can have any feasible depth of nesting.
#' 
#' The treatments are defined by the \code{treatments} and \code{replicates} parameter lists where the \code{treatments} list partitions the
#' total number of treatments into sub-sets of equally replicated treatments and the \code{replicates} list provides the replication for each set. The two lists are of equal length and 
#' the treatment numbers and replication numbers are paired according to their index position in the two lists. Any feasible number of replication levels can be defined and the resulting treatments are 
#' labelled consecutively according to the ordering of the the treatment sets in the two lists.
#'  
#' The \code{blocklevels} list defines the blocks structure of the design where the first number in the list is the number of main blocks 
#' and the succesive numbers, if any, are the numbers of succesively nested blocks. The numbers are the nested blocks in each preceding block therefore the cumulative product of the numbers
#' is the total number of blocks in the design. The default value for the \code{blocklevels} list is a single number equal to the highest common factor (hcf) of
#' the replication numbers, which gives an orthogonal complete blocks design with the maximum possible number of othogonal blocks.   
#'  
#' Block sizes in any given stratum are equal if the cumulative number of blocks in that stratum exactly divides the total number of plots, otherwise they differ by not more than a single unit. 
#' 
#' Designs with t = v**2 equally replicated treatments, k = v+1 or fewer replicates and blocks of size v where v is a prime or prime-power if k>3 or is arbitrary otherwise,
#' are regular lattice designs and can be constructed algebraically. \code{blocks} constructs all regular lattice designs algebraically and also constructs the special 
#' non-prime lattice with t = 100, v = 10 and k< = 4 algebraically. Prime-power lattice designs require the \code{MOLS} function of the \code{crossdes} package.
#' All other designs are constructed algorithmically by a swapping algorithm that maximizes the determinant of the information matrix (D-optimality). 
#'  
#' Designs are fully randomized with each set of nested blocks randomized within the preceding set of blocks and with treatments fully randomized within the bottom set of blocks.
#'  
#' @param treatments a list of numbers for the number of treatments in each treatment replication set of the design. Each treatment number must have a matching
#' replication number in the \code{replicates} list.
#' 
#' @param replicates a list of numbers for the number of replicates in each treatment replication set in the design. Each replication number must 
#' have a matching treatment number in the \code{treatments} list.
#' 
#' @param blocklevels an optional list of nested blocks where the first number is the number of main blocks and where the  remaining numbers, if any,
#' are the numbers of blocks nested in each preceding block. The default is the hcf of the replication numbers which gives a maximum set of complete orthogonal main blocks.
#' 
#' @param searches an optional number for the number of local optima searched during optimization. The default is the minimum of 64 or the integer quotient of 4096 
#' divided by the number of plots.
#' 
#' @param seed an integer seed for initializing the random number generator where a design must be reproducible. The default is a random seed.
#' 
#' @return  
#' \item{Design}{Data frame showing the block and treatment factors for each sub-plot}
#' \item{Plan}{Data frame showing the allocation of treatments to sub-plots for each block in the design}
#' \item{Incidences}{List of blocks-by-treatments incidence matrices, one for each stratum in the design}
#' \item{Efficiencies}{Data frame showing the A-efficiency factor for each stratum in the design together with an upper bound, where available}
#' \item{seed}{Numerical seed for random number generator}
#'
#' @references
#' 
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. http://CRAN.R-project.org/package=crossdes
#'   
#' @examples
#' 
#' # 3 treatments with 2 reps, 2 treatments with 4 reps, 4 treatments with 3 reps 
#' # the hcf of the replication numbers is 1 and the default design is a completely randomized design 
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 50 treatments with 4 reps in 4 complete randomized blocks 
#' blocks(treatments=50,replicates=4)
#' 
#' # as above but with 4 main blocks and 5 nested blocks within each main block 
#' blocks(treatments=50,replicates=4,blocklevels=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments, one to each block
#' blocks(treatments=c(50,20),replicates=c(4,1),blocklevels=c(4,5))
#' 
#' # 64 treatments with 2 reps and 2 main blocks with five 2-level nested factors   
#' blocks(treatments=64,replicates=2,blocklevels=c(2,2,2,2,2,2),searches=4)
#' 
#' # concurrence matrices of 36 treatments with 3 reps and 3 main blocks with 6 nested blocks
#' crossprod(blocks(treatments=36,replicates=3,blocklevels=c(3,6))$Incidences[[2]])
#' 
#' # concurrence matrix for 13 treatments with 4 reps and 13 treatments with one rep in 13 blocks 
#' crossprod(blocks(c(13,13),c(4,1),13,searches=100)$Incidences[[1]])
#' 
#' @export
#' 

blocks = function(treatments, replicates, blocklevels=hcf, searches=min(64, floor(4096/nunits)), seed=NULL) { 
  
  Sizes=function(mainSizes,slevs) {
    newsizes=vector(length=slevs*length(mainSizes))
    for (z in 1: length(mainSizes)) {
      bsize = mainSizes[z] %/% slevs
      resid = mainSizes[z] %% slevs
      for (i in 1 : slevs) {
        newsizes[(z-1)*slevs+i] = bsize + (resid>0)
        resid=resid-1
      }
    }
    newsizes
  }
  
  isPrime=function(v) {
    if (v <= 3)  return(TRUE)
    else if (v %% 2 == 0 | v %% 3 == 0) return(FALSE)
    else if (v<25) return(TRUE)
    else 
      for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
        if( v %% (i-1) == 0 | v %% (i+1) == 0) return(FALSE) 
    return(TRUE)
  }      
  
  BlockContrasts=function(MF,BF) {
    BM=matrix(0,nrow=length(BF),ncol=nlevels(BF))
    BM[cbind(1:length(BF),as.numeric(BF))]=1 # factor indicator matrix
    BM=BM[,rep(c(rep(TRUE,((nlevels(BF)/nlevels(MF))-1)),FALSE),nlevels(MF)),drop=FALSE]
    for (i in 1: nlevels(MF)) 
      BM[(as.numeric(MF)==i),]=scale(BM[(as.numeric(MF)==i),], center = TRUE, scale = FALSE) # scaling within main blocks
    BM
  } 
  
  TreatContrasts=function(MF,TF) {
    TM=matrix(0,nrow=length(TF),ncol=nlevels(TF))
    TM[cbind(1:length(TF),as.numeric(TF))]=1 # factor indicator matrix	
    TM=TM[,c(rep(TRUE,(nlevels(TF)-1)),FALSE),drop=FALSE]	
    for (i in 1:nlevels(MF)) 
      TM[MF==i,]=scale(TM[MF==i,] , center = TRUE, scale = FALSE)
    TM
  }
  
  DetMat=function(M11,M22,M12,Samp,TF,RF) {
    detmat=vector("list", length(Samp))
    for (i in 1:length(Samp)) {
      trts=as.numeric(TF[Samp[[i]]])
      rblks=as.numeric(RF[Samp[[i]]])
      TT=M11[trts,trts,drop=FALSE]
      BB=M22[rblks,rblks,drop=FALSE]
      TB=M12[trts,rblks,drop=FALSE]
      TT=TT-crossprod(t(diag(TT)),t(rep(1,ncol(TT))))
      TT=TT+t(TT)
      BB=BB-crossprod(t(diag(BB)),t(rep(1,ncol(TT))))
      BB=BB+t(BB)
      TB=TB-crossprod(t(diag(TB)),t(rep(1,ncol(TT))))
      TB=1+TB+t(TB)
      detmat[[i]]=TB**2-TT*BB
    }
    detmat
  }
  
  UpDate=function(M11,M22,M12,si,sj,TF,RF) {
    M11T = M11[TF[si],]-M11[TF[sj],]
    M12T = M12[TF[si],]-M12[TF[sj],]
    M12B = M12[,RF[si]]-M12[,RF[sj]]
    M22B = M22[,RF[si]]-M22[,RF[sj]]
    # updating vectors
    m11=M11T[TF[si]]-M11T[TF[sj]]
    m22=M22B[RF[si]]-M22B[RF[sj]]
    m12=M12T[RF[si]]-M12T[RF[sj]]
    f = sqrt(2+m11+m22-2*m12)
    m = f/sqrt(1-2*m12-m11*m22+m12*m12)/2
    Z1 = (M12B-M11T)/f 
    Z2 = (M22B-M12T)/f 
    W1 = (M11T+M12B - Z1*(m22-m11)/f)*m
    W2 = (M12T+M22B - Z2*(m22-m11)/f)*m
    M11 = M11 - crossprod(t(Z1)) + crossprod(t(W1))
    M22 = M22 - crossprod(t(Z2)) + crossprod(t(W2))
    M12 = M12 - crossprod(t(Z1),t(Z2)) + crossprod(t(W1),t(W2))
    up=list(M11=M11,M22=M22,M12=M12)
  } # end of function
  
  
  
  
  Optimise=function(TF,BF,MF,M11,M22,M12,searches)   {
    # first stage finds an optima by optimizing samples of increasing size in powers of 2
    globrelD=1
    locrelD=1
    globTF=TF
    nunits=length(TF)
    mainSets=split(rep(1:nunits),MF)
    Samp=vector("list", nlevels(MF)) 
    for (r in 1 : searches) {
      nSamp=ceiling(  min(nunits,36)*tabulate(MF)/nunits) #sample size is smallest of 36 or nunits or has at least one sample per restriction
      repeat {
        relD=0
        for (i in 1:nlevels(MF)) Samp[[i]]=sort(sample(mainSets[[i]],nSamp[i]))						
        detmat=DetMat(M11,M22,M12,Samp,TF,BF)
        grelD=1
        for (i in 1:nlevels(MF)) {
          N=which.max(detmat[[i]])
          ti=1+(N-1)%%nrow(detmat[[i]])
          tj=1+(N-1)%/%nrow(detmat[[i]])
          relD=detmat[[i]][ti,tj]
          if (relD>grelD) {
            gsi=Samp[[i]][ti]
            gsj=Samp[[i]][tj]
            grelD=relD
          }	
        }			
        if (grelD>1.00000001) {
          up=UpDate(M11,M22,M12,gsi,gsj,TF,BF)
          M11=up$M11
          M22=up$M22
          M12=up$M12
          TF[c(gsi,gsj)]=TF[c(gsj,gsi)]
          locrelD=grelD*locrelD				
        } else if ( sum(nSamp) < min(nunits,512)) {
          nSamp=2*nSamp # doubles sample size
          if (sum(nSamp)>nunits) nSamp=tabulate(MF) # ensures sample size not greater than the population
        } else break
      } # repeat	
      if (locrelD>globrelD) {
        globrelD=locrelD
        globTF=TF
      }
      if (searches>r) {
        # escape local optima
        prop_change=1
        for (iswap in 1 : 6) {
          icount=0
          dswap=1
          while (icount<100 & (dswap<0.01 | dswap>.999)) {
            icount=icount+1
            s=sample(rep(1:nunits)[MF==sample(nlevels(MF),1)],2)
            # calculates the proportional change in the determinant of the design information due to swapping treatments on plots s1 and s2
            if ( (TF[s[1]]!=TF[s[2]]) & (BF[s[1]]!=BF[s[2]]) ) 	
              dswap=(1+M12[TF[s[1]],BF[s[2]]]+M12[TF[s[2]],BF[s[1]]]-M12[TF[s[1]],BF[s[1]]]-M12[TF[s[2]],BF[s[2]]])**2-
              (2*M11[TF[s[1]],TF[s[2]]]-M11[TF[s[1]],TF[s[1]]]-M11[TF[s[2]],TF[s[2]]])*
              (2*M22[BF[s[1]],BF[s[2]]]-M22[BF[s[1]],BF[s[1]]]-M22[BF[s[2]],BF[s[2]]])
          }
          #updates matrices
          if (icount<100) {
            prop_change=prop_change*dswap
            up=UpDate(M11,M22,M12,s[1],s[2],TF,BF)  
            M11=up$M11
            M22=up$M22
            M12=up$M12
            TF[c(s[1],s[2])]=TF[c(s[2],s[1])]	
          }
        } 
        locrelD=locrelD*prop_change
      } 
    } # next for
    globTF
  } # end of function
  
  FullRank=function(MF,TF,TB,NB,DD) {
    count=0
    while ( (qr(DD)$rank<ncol(DD))&(count<1000) ) {
      rand=sample(1:length(TF))
      TF=TF[rand][order(MF[rand])]
      TB=TB[rand,][order(MF[rand]),]
      DD=crossprod(cbind(TB,NB))
      count=count+1
    } # end of while
    fullrank=list(TF=TF,DD=DD,count=count)	
  }
  
  GenOpt=function(TF,NF,MF,searches)  {
    singular=FALSE
    T=TreatContrasts(MF,TF)
    N=BlockContrasts(MF,NF) 
    DD=crossprod(cbind(T,N))
    count=0
    # attempts to find a full rank randomisation	
    if (qr(DD)$rank<ncol(DD)) {
      fullrank=FullRank(MF,TF,T,N,DD)
      TF=fullrank$TF
      count=fullrank$count
      DD=fullrank$DD
    }
    if (count<100) {
      V=chol2inv(chol(DD))
      M11=matrix(0,nrow=nlevels(TF),ncol=nlevels(TF))	
      M22=matrix(0,nrow=nlevels(NF),ncol=nlevels(NF))
      M12=matrix(0,nrow=nlevels(TF),ncol=nlevels(NF))
      M11[1:(nlevels(TF)-1),1:(nlevels(TF)-1)]=V[1:(nlevels(TF)-1),1:(nlevels(TF)-1),drop=FALSE]
      M12[1:(nlevels(TF)-1),1:(ncol(V)-nlevels(TF)+1)]=V[1:(nlevels(TF)-1),nlevels(TF):ncol(V),drop=FALSE]
      M22[1:(ncol(V)-nlevels(TF)+1),1:(ncol(V)-nlevels(TF)+1)]=V[nlevels(TF):ncol(V),nlevels(TF):ncol(V),drop=FALSE]
      sortR=c(1:nlevels(NF))%%(nlevels(NF)/nlevels(MF))==0
      sortC=NULL
      sortN=NULL
      perm=order(order(c(sortR,sortC,sortN)))
      M12=M12[,perm]
      M22=M22[perm,perm]	
      TF=Optimise(TF,NF,MF,M11,M22,M12,searches)
    }	
    as.factor(TF)
  }
  
  if (missing(treatments) | missing(replicates) )  return(" Treatments or replicates not defined ")   
  if (is.null(treatments) | is.null(replicates))  return(" Treatments or replicates list is empty ") 	
  if (anyNA(treatments) | anyNA(replicates) ) return(" NA values not allowed")
  if (!all(is.finite(treatments)) | !all(is.finite(replicates)) | !all(!is.nan(treatments)) | !all(!is.nan(replicates))) return(" Treatments and replicates must contain only finite integers ")
  if ( length(treatments)!=length(replicates) ) return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
  if (is.null(seed)) seed=sample(1:100000,1)
  set.seed(seed)
  options(warn=-1)	
  # omit any single replicate treatments here 
  if (all(replicates==1)) {
    treatlevs=treatments
    replevs = replicates
  } else {
    treatlevs=treatments[replicates>1]
    replevs = replicates[replicates>1]
  }
  nunits=sum(treatlevs*replevs) 
  hcf=replevs[1]
  if (length(replevs)>1)
    for (i in 2: length(replevs)) {
      v=sort(c(replevs[i],hcf))
      while (v[2]%%v[1] != 0) v = c(v[2]%%v[1], v[1])
      hcf=v[1]        
    }
  if (anyNA(blocklevels) | anyNA(searches) ) return(" NA values not allowed") 
  if (!all(is.finite(blocklevels)) | !all(is.finite(searches)) | !all(!is.nan(blocklevels)) | !all(!is.nan(searches))) return(" Inputs must contain only finite integers ")
  if (searches<1)  return(" Repeats must be at least one ") 	
  if (  sum(treatments*replicates) < (prod(blocklevels) + sum(treatments)-1) ) return("Design cannot be fitted :  too many blocks and treatments for the available plots")	
  ntrts=sum(treatlevs)
  strata=length(blocklevels)	
  cumblocklevs=c(1,blocklevels)	
  for (i in 1 : strata) cumblocklevs[i+1]=cumblocklevs[i+1]*cumblocklevs[i]
  blocksizes=nunits
  
  for (i in 1:strata)
    blocksizes=Sizes(blocksizes,blocklevels[i])
  
  # treps is the vector of treatment replications for the minimum orthogonal block size
  treps=rep(replevs,treatlevs)/hcf  
  # facMat is the matrix of nested factor levels NOT sub-levels
  TF=as.factor(rep(rep(1:ntrts,treps),hcf))
  facMat= matrix(nrow=cumblocklevs[strata+1],ncol=strata)
  for (r in 1 : strata) 
    facMat[,r]=gl(cumblocklevs[r+1],cumblocklevs[strata+1]/cumblocklevs[r+1])
  desMat= matrix(1,nrow=nunits,ncol=(strata+1))
  for (r in 1 : strata) 
    desMat[,(r+1)]=as.factor(rep(facMat[,r],blocksizes)) 
  orthbsize=nunits/hcf 
  ortho=0
  for (i in 1 : strata) 
    if (all( tabulate(desMat[,i+1]) %% orthbsize == 0)) ortho=i else  break
  
  # genopt (current treatment factor, current block factor, restriction block factor in the next higher stratum   
  regreps=(max(replevs)==min(replevs))
  v=sqrt(ntrts)
  pp_trts=c(16,64,256,1024,4096,16384,81,729,6561,625,2401)	
  lattice = ( replevs[1]<4 |  (replevs[1]<(v+2) & isPrime(v) ) )
  
  if (ortho<strata) {
    for (i in (ortho+1) : strata) {
      lattice_i = (regreps & identical(v,ntrts%/%v) &  cumblocklevs[i+1]==v*replevs[1])
      if ( i==(ortho+1) & lattice & lattice_i ) {			
        TF=rep(1:(v*v))
        TF=c(TF, rep(1:(v*v))[order(rep(0:(v-1),v))])
        for (z in 1: (replevs[1]-2)) {
          set=NULL
          for (j in 0: (v-1)) 
            for (k in 0: (v-1)) 
              set=c(set,(j+k*z)%%v)
          TF=as.factor(c(TF, rep(1:(v*v))[order(set)]))		
        }				
      } else if (i==(ortho+1) & lattice_i & replevs[1]<(v+2) & ntrts%in%pp_trts ) {		
        prime= c(2,2,2,2,2,2,   3,3,3,  5,7)[which(pp_trts==ntrts)]
        ppower=c(2,3,4,5,6,7,   2,3,4,  2,2)[which(pp_trts==ntrts)]
        mols=crossdes::MOLS(prime,ppower)			
        TF=rep(1:(v*v))
        TF=c(TF, rep(1:(v*v))[order(rep(0:(v-1),v))])
        for (i in 1: (replevs[1]-2))
          TF=c(TF, rep(1:(v*v))[order(    as.numeric(mols[,,i]) ) ])
        TF=as.factor(TF)
      } else if (i==(ortho+1) & lattice_i & ntrts==100 & replevs[1]<=4) {			
        TF=rep(1:(v*v))
        TF=c(TF, rep(1:(v*v))[order(rep(0:(v-1),v))])
        tens1=c(
          1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
          5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8)
        tens2=c(
          1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
          6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6)
        if (replevs[1]>2)  TF=c(TF, rep(1:(v*v))[order(tens1)])	
        if (replevs[1]==4) TF=c(TF, rep(1:(v*v))[order(tens2)])
        TF=as.factor(TF)
      } else {		
        TF=GenOpt(TF,as.factor(desMat[,(i+1)]),as.factor(desMat[,i]),searches)
      }
    }
  }
  
  # add back single replicate treatments here 
  if ( !all(replicates>1) & !all(replicates==1) ) {
    nunits=sum(treatments*replicates)
    TF=c(TF,sample(c( (sum(treatlevs)+1):sum(treatments)) )	)
    fullblocksizes=nunits
    for (i in 1:strata)
      fullblocksizes=Sizes(fullblocksizes,blocklevels[i])	  
    BF=c( rep( 1:length(blocksizes),blocksizes),  rep( 1:length(blocksizes),(fullblocksizes-blocksizes) ) )
    # full TF in blocks
    TF=TF[order(BF)]
    # new desMat for blocksizes
    blocksizes=fullblocksizes	
    # desMat is a nested design matrix where column 1 is the mean and the remaining columns are the blocklevels in each stratum 
    desMat= matrix(1,nrow=nunits,ncol=(strata+1))
    for (r in 1 : strata) 
      desMat[,(r+1)]=as.factor(rep(facMat[,r],blocksizes))
    hcf=1
    ortho=0
  }	
  
  # randomize blocks within each nested stratum
  randMat= matrix(1,nrow=nunits,ncol=strata)
  rblocks=as.data.frame(facMat)
  rblocks[]=lapply(rblocks, factor)
  for (r in 1 : strata)
    levels(rblocks[,r])=sample(nlevels(rblocks[,r]))
  for (r in 1 : strata) 
    randMat[,r]=rep(as.numeric(levels(rblocks[,r]))[rblocks[,r]] ,blocksizes)
  dd=as.data.frame(cbind(randMat,sample(rep(1:nunits)),TF))
  dd$TF=as.factor(dd$TF)
  TF=dd[ do.call(order, dd), ]$TF
  
  # Design data frame
  plots=NULL
  for (i in 1: length(blocksizes))
    for (j in 1:blocksizes[i]) 
      plots=c(plots,j)
  Design=as.data.frame(cbind(desMat[,c(2:ncol(desMat))],plots,TF))
  Design[]=lapply(Design, factor)	
  designnames="Main"
  if (strata>1)
    for (i in 1:(strata-1))
      designnames=c(designnames,paste("Sub_",i,sep=""))
  colnames(Design)=c(designnames,"Sub_plots","Treatments")
    
  #Design plan layout
  index=1
  d=matrix(nrow=length(blocksizes),ncol=max(blocksizes))
  for (i in 1: length(blocksizes))
    for (j in 1:blocksizes[i]) {
      d[i,j]=TF[index]
      index=index+1
    }
  d[is.na(d)]  = " "
  Plan=as.data.frame(cbind(facMat,rep(" ",nrow(d)),d))
  designnames=c(designnames,"Sub_plots")
  for (i in 1:max(blocksizes))
    designnames=c(designnames,i)
  colnames(Plan)=designnames
  
  # Incidence matrix for each stratum
  Incidences=vector(mode = "list", length =strata )
  for (i in 1:strata)
    Incidences[[i]]=table(Design[,c(i,ncol(Design))])  
  
  # Efficiencies data frame
  bounds=rep(0,strata)
  r=1/sqrt(tabulate(TF))
  aeff=c(rep(0,strata))
  for (i in 1:strata) {	
    k=1/sqrt(tabulate(Design[,i]))
    U=crossprod(t(crossprod(diag(r,nrow = length(r)),table(TF,Design[,i]))),diag(k,nrow = length(k)))
    A=diag(length(r))-crossprod(t(U))
    aeff[i]=1/mean(1/eigen(A, symmetric=TRUE, only.values = TRUE)$values[1:length(r)-1])
    if (max(replicates)==min(replicates))
      bounds[i]=upper_bounds(nunits,ntrts,cumblocklevs[i+1]) 
    else
      bounds[i]=NA
  }
  effics=as.data.frame(cbind( cumblocklevs[2:(strata+1)], aeff, bounds))
  rnames=c("Main")
  if (strata>1)
    for (i in 1 : (strata-1)) rnames=c(rnames,paste("Sub",i))
  colnames(effics)=c("Blocks","A-Efficiencies", "Upper Bounds")
  rownames(effics)=rnames
  
  list(Design=Design,Plan=Plan,Incidences=Incidences,Efficiencies=effics,seed=seed)
} 