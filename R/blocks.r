#' @title Block designs 
#' 
#' @description
#' 
#' \code{blocks} function constructs nested blocks for unstructured treatments with arbitrary levels of replication and arbitrary depth of nesting
#' 
#' @details
#' 
#' \code{blocks} constructs nested block designs for unstructured treatment sets where treatments can have any arbitrary replication, not necessarily all equal, 
#' and blocks can have any feasible depth of nesting.
#' 
#' Treatment and replication numbers are defined by the \code{treatments} and \code{replicates} parameter lists. These lists must be of equal length
#' and each matching pair of numbers in the two lists represents a treatment set where the \code{treatments} list gives the number of treatments in the set and 
#' the \code{replicates} list gives the replication of the set.
#'  
#' Any number of treatment sets is allowed and the treatments are numbered consecutively according to the ordering of the treatment sets in the parameter lists
#' (see the examples). 
#'  
#' Blocks are defined by the \code{blocklevels} list which is a hierarchical list of nested blocks. The first number is the number of main blocks 
#' and the succesive numbers, if any, are the numbers of blocks nested in each preceding block. The cumulative product of the levels for any stratum
#' is the total number of blocks in that stratum. The default value for the \code{blocklevels} list is a single number equal to the highest common factor (hcf) of
#' the replication numbers, which gives an orthogonal blocks design with the maximum possible number of othogonal blocks.   
#'  
#' Block sizes in any given stratum will be equal if the cumulative number of blocks exactly divides the number of plots otherwise
#' they will be as near equal as possible and will never differ by more than a single unit. 
#' 
#' Lattice designs where the number of treatments is the square of the block size v and the number of replicates is k+2 or less and k mutually 
#' orthogonal latin squares (MOLS) of size v*v exist are constructed algebraically. Lattice designs exist for any 
#' v if k = 1 and for any prime or prime power v if k < v and for v = 10 and k = 2. Prime-power MOLS are constructed by using the MOLS 
#' function of the \code{crossdes} package (Sailer 2013).  
#' 
#' All other designs are constructed algorithmically by a swapping algorithm that maximizes the determinant of the information matrix (D-optimality). 
#'  
#' Designs are optimized hierarchically with the blocks of each new set optimized within the blocks of the preceding set.
#'  
#' Designs are fully randomized with treatments randomized within blocks and each set of nested blocks randomized within the
#' preceding set of blocks.
#'  
#' @param treatments A list of the number of treatments for each treatment set in the design. Each treatment number must have a matching
#' replication number in the \code{replicates} list
#' 
#' @param replicates A list of the number of replicates for each treatment set in the design. Each replication number must 
#' have a matching treatment number in the \code{treatments} list
#' 
#' @param blocklevels A hierarchical list of nested blocks where the first number is the number of main blocks and the  remaining numbers, if any,
#' are the numbers of blocks nested in each preceding block. The default is the hcf of the replication numbers.
#' 
#' @param searches The number of local optima searched during 
#' a design optimization. The default is the minimum of 64 or the integer quotient of 4096 divided by the number of plots.
#' 
#' @param seed An integer seed for initializing the random number generator where a design must be reproducible. The default is a random seed.
#' 
#' @return  
#' \item{Design}{Data frame showing the block and treatment factors for each plot}
#' \item{Plan}{Data frame showing the allocation of treatments to plots for each block in the design}
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
#' # the replication hcf is 1 and the default design is a completely randomized design 
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 50 treatments with 4 reps in 4 complete randomized blocks 
#' blocks(treatments=50,replicates=4)
#' 
#' # as above but with 4 main blocks and 5 nested blocks in each main block 
#' blocks(treatments=50,replicates=4,blocklevels=c(4,5))
#' 
#' # as above but with 20 additional single replicate treatments, one to each block
#' blocks(treatments=c(50,20),replicates=c(4,1),blocklevels=c(4,5))
#' 
#' # 64 treatments with 2 reps and 2 main blocks with five 2-level nested factors   
#' blocks(treatments=64,replicates=2,blocklevels=c(2,2,2,2,2,2),searches=12)
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
  if (missing(treatments) | missing(replicates) )  return(" Treatments or replicates not defined ")   
  if (is.null(treatments) | is.null(replicates))  return(" Treatments or replicates list is empty ") 	
  if (anyNA(treatments) | anyNA(replicates) ) return(" NA values not allowed")
  if (!all(is.finite(treatments)) | !all(is.finite(replicates)) | !all(!is.nan(treatments)) | !all(!is.nan(replicates))) return(" Treatments and replicates must contain only finite integers ")
  if ( length(treatments)!=length(replicates) ) return(paste("The number of treatments sets = " , length(treatments) , " does not equal the number of replication sets = " , length(replicates)))
  if (is.null(seed)) seed=sample(1:100000,1)
  set.seed(seed)
  # omit any single replicate treatments here unless all single replicate	
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
  blocklevels=c(1,blocklevels)		
  strata=length(blocklevels)
  cumblocklevs=blocklevels	
  for (i in 2 : strata) cumblocklevs[i]=cumblocklevs[i-1]*cumblocklevs[i]
  blocksizes=nunits		
  for (i in 2 : strata) {
    subSizes=vector(length=(length(blocksizes)*blocklevels[i]))
    for (z in 1: length(blocksizes)) {
      bsize = blocksizes[z] %/% blocklevels[i]
      bresid = blocksizes[z] %% blocklevels[i]
      for (w in 1 : blocklevels[i]) {
        subSizes[(z-1)*blocklevels[i]+w] = bsize + (bresid>0)
        bresid=bresid-1
      }
    }
    blocksizes=subSizes
  }
  # treps is the treatment replications for the minimum orthogonal block size excluding unreplicated treatments
  treps=rep(replevs,treatlevs)/hcf 
  # Trts is an initial minimum non-randomised complete blocks design
  Trts=as.factor(rep(rep(1:ntrts,treps),hcf))
  facMat= matrix(nrow=cumblocklevs[strata],ncol=strata)
  for (r in 1 : strata) 
    facMat[,r]=gl(cumblocklevs[r],cumblocklevs[strata]/cumblocklevs[r])	
  desMat= matrix(nrow=nunits,ncol=strata)
  for (r in 1 : strata) 
    desMat[,r]=rep(facMat[,r],blocksizes)		
  orthbsize=nunits/hcf 	
  ortho=1	
  for (i in 2 : strata) 
    if (all( tabulate(desMat[,i]) %% orthbsize == 0)) ortho=i else break
  v=sqrt(ntrts)		
  lattice= (max(replevs)==min(replevs) & identical(v,ntrts%/%v))	  
  
  # primality test function
  isPrime=function(v) {
    if (v <= 3)  return(TRUE)
    else if (v %% 2 == 0 | v %% 3 == 0) return(FALSE)
    else if (v<25) return(TRUE)
    else 
      for(i in  6*rep(1:floor((sqrt(v)+1)/6)) )
        if( v %% (i-1) == 0 | v %% (i+1) == 0) return(FALSE) 
    return(TRUE)
  }      
    
  #updates matrices
  UpDate=function(M11,M22,M12,si,sj,Trts,BF) {
    M11T = M11[Trts[si],]-M11[Trts[sj],]
    M12T = M12[Trts[si],]-M12[Trts[sj],]
    M12B = M12[,BF[si]]-M12[,BF[sj]]
    M22B = M22[,BF[si]]-M22[,BF[sj]]
    m11=M11T[Trts[si]]-M11T[Trts[sj]]
    m22=M22B[BF[si]]-M22B[BF[sj]]
    m12=1-(M12T[BF[si]]-M12T[BF[sj]])	
    f = sqrt(m11+m22+2*m12)
    m = f/sqrt(m12*m12-m11*m22)/2
    Z1 = (M12B-M11T)/f 
    Z2 = (M22B-M12T)/f 
    W1 = (M11T+M12B - Z1*(m22-m11)/f)*m
    W2 = (M12T+M22B - Z2*(m22-m11)/f)*m	
    M11 = M11 - crossprod(t(Z1)) + crossprod(t(W1))
    M22 = M22 - crossprod(t(Z2)) + crossprod(t(W2))
    M12 = M12 - crossprod(t(Z1),t(Z2)) + crossprod(t(W1),t(W2))	
    list(M11=M11,M22=M22,M12=M12)
  }	
    
  pp_trts=c(16,64,256,1024,4096,16384,81,729,6561,625,2401)	
  if (ortho<strata) {
    regreps=nunits/ntrts # replication for equireplicate design	
    for (i in (ortho+1) :strata) {	
      
      if ( i==(ortho+1) & lattice  & cumblocklevs[i]==v*regreps &  ( regreps<=3   |  (regreps<=(v+1) & isPrime(v) ) ) ) {
        mols=c( rep(sample(0:(v-1)),v),rep(sample(v:(2*v-1)),each=v))
        square=vector(length=(v*v))
        if (regreps>2) {
          rows=sample(0:(v-1))
          cols=sample(0:(v-1))
          for (i in 1: (regreps-2)) {
            index=1
            for (j in rows) 
              for (k in cols) {
                square[index]=(j+k*i)%%v	
                index=index+1
              }			
            mols=c(mols,(square + v*(i+1)))
          }	
        }
        Trts=rep(sample(1:(v*v)),regreps)[order(mols)]
        rand=sample(1:(v*v*regreps))
        Trts=as.factor(Trts[rand][order(rep(1:(v*regreps),each=v)[rand])]) # randomizes plots within sub-blocks	
        
      } else if (i==(ortho+1) & lattice & cumblocklevs[i]==v*regreps & regreps<=(v+1) & ntrts%in%pp_trts ) {								
        prime=c(2,2,2,2,2,2,   3,3,3,  5,7)[which(pp_trts==ntrts)]
        ppower=c(2,3,4,5,6,7,   2,3,4,  2,2)[which(pp_trts==ntrts)]
        mols=c( rep(sample(0:(v-1)),v),rep(sample(v:(2*v-1)),each=v))	
        fullmols=crossdes::MOLS(prime,ppower)[sample(rep(1:v)),sample(rep(1:v)),]	
        for (i in 1: (regreps-2)) mols=c(mols,(as.numeric(fullmols[,,i]) + v*(i+1) -1))
        Trts=rep((1:(v*v)),regreps)[order(mols)]
        rand=sample(1:(v*v*regreps))
        Trts=as.factor(Trts[rand][order(rep(1:(v*regreps),each=v)[rand])]) # randomizes plots within sub-blocks	
        
      } else if (i==(ortho+1) & lattice & ntrts==100 & regreps<=4 & cumblocklevs[i]==regreps*v) {	
        mols=c( rep(sample(0:9),10),rep(sample(10:19),each=10))
        tens=c(
          1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
          5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8, 
          1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
          6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6)
        rand=sample(1:200)
        tens=tens[rand][order(rep(1:20,each=10)[rand],rep(1:10,20)[rand])]
        if (regreps>2) mols=c(mols,(tens[1:100] + 20))
        if (regreps==4) mols=c(mols,(tens[101:200] + 30))
        Trts=rep(sample(1:100),regreps)[order(mols)]
        rand=sample(1:(100*regreps))
        Trts=as.factor(Trts[rand][order(rep(1:(10*regreps),each=10)[rand])]) # randomize plots in sub-blocks
        # algorithmic
      } else {	
        MF=as.factor(desMat[,(i-1)])	
        BF=as.factor(desMat[,i])	
        M=matrix(0,nrow=length(MF),ncol=nlevels(MF))
        M[cbind(rep(1:length(MF)),MF)]=1
        IM=diag(length(MF))-crossprod(crossprod(diag(1/sqrt(as.numeric(table(MF))),nrow=nlevels(MF)),t(M)))
        T=matrix(0,nrow=length(MF),ncol=nlevels(Trts))
        T[cbind(rep(1:length(MF)),Trts)]=1
        T=crossprod(IM,T)[,rep(1:(nlevels(Trts)-1))]
        B=matrix(0,nrow=length(MF),ncol=nlevels(BF))
        B[cbind(rep(1:length(MF)),BF)]=1
        B=crossprod(IM,B)[,rep( c(rep(TRUE,((nlevels(BF)/nlevels(MF))-1)),FALSE) ,nlevels(MF))]
        DD=crossprod(cbind(T,B))
        ndim=ncol(DD)
        count=0
        # attempts to find a full rank randomisation
        while ( (qr(DD)$rank<ndim)&(count<1000) ) {
          rand=sample(1:length(Trts))
          Trts=Trts[rand][order(MF[rand])]
          T=T[rand,][order(MF[rand]),]
          DD=crossprod(cbind(T,B))
          count=count+1
        } 
        if (count<1000) {
          V=chol2inv(chol(DD))
          M11=matrix(0,nrow=nlevels(Trts),ncol=nlevels(Trts))	
          M22=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
          M12=matrix(0,nrow=nlevels(Trts),ncol=nlevels(BF))
          M11[1:(nlevels(Trts)-1),1:(nlevels(Trts)-1)]=V[1:(nlevels(Trts)-1),1:(nlevels(Trts)-1),drop=FALSE]
          M12[1:(nlevels(Trts)-1),1:(ncol(V)-nlevels(Trts)+1)]=V[1:(nlevels(Trts)-1),nlevels(Trts):ncol(V),drop=FALSE]
          M22[1:(ncol(V)-nlevels(Trts)+1),1:(ncol(V)-nlevels(Trts)+1)]=V[nlevels(Trts):ncol(V),nlevels(Trts):ncol(V),drop=FALSE]
          perm=order(order(c(1:nlevels(BF))%%(nlevels(BF)/nlevels(MF))==0))
          M12=M12[,perm]
          M22=M22[perm,perm]
          globrelD=1
          locrelD=1
          globTF=Trts         
          mainSets=split(rep(1:nunits),MF)
          for (r in 1 : searches) {					
            nSamp=ceiling(  min(nunits,36)*tabulate(MF)/nunits) # assuming initial sample size is smallest of 36 or nunits or with at least one sample per restriction block
            repeat {
              relD=0						
              for (i in 1:nlevels(MF)) {
                Samp=sort(sample(mainSets[[i]],nSamp[i]))
                TT=M11[Trts[Samp],Trts[Samp],drop=FALSE]
                BB=M22[BF[Samp],BF[Samp],drop=FALSE]
                TB=M12[Trts[Samp],BF[Samp],drop=FALSE]
                TT=TT-crossprod(t(diag(TT)),t(rep(1,ncol(TT))))
                TT=TT+t(TT)
                BB=BB-crossprod(t(diag(BB)),t(rep(1,ncol(TT))))
                BB=BB+t(BB)
                TB=TB-crossprod(t(diag(TB)),t(rep(1,ncol(TT))))
                TB=1+TB+t(TB)
                detmat=TB**2-TT*BB
                maxD=max(detmat)						
                if (maxD>relD & maxD>1.00001) {
                  N=which.max(detmat)
                  ti=1+(N-1)%%nrow(detmat)
                  tj=1+(N-1)%/%nrow(detmat)
                  si=Samp[ti]
                  sj=Samp[tj]	
                  relD=maxD
                }
              }
              if (relD>0) {
                up=UpDate(M11,M22,M12,si,sj,Trts,BF)  
                M11=up$M11
                M22=up$M22
                M12=up$M12
                Trts[c(si,sj)]=Trts[c(sj,si)]								
                locrelD=relD*locrelD	
              } else if ( sum(nSamp) < min(nunits,512)) {
                nSamp=2*nSamp # doubles sample size
                if (sum(nSamp)>nunits) nSamp=tabulate(MF) # ensures sample size not greater than the population
              } else break
            } # repeat	
            
            if (locrelD>globrelD) {
              globrelD=locrelD
              globTF=Trts
            }						
            # next search
            if (searches>r) {
              # escape local optima
              prop_change=1
              for (iswap in 1 : 6) {
                icount=0
                dswap=0
                while (icount<100 & dswap<0.01) {
                  icount=icount+1
                  s1=sample(1:length(Trts),1)
                  if (all(!( (MF==MF[s1]) & (Trts!=Trts[s1]) & (BF!=BF[s1])))) next		
                  s2 = sample( rep(1:length(Trts))[ ( (MF==MF[s1]) & (Trts!=Trts[s1]) & (BF!=BF[s1]) ) ],1)
                  # calculates the proportional change in the determinant of the design information due to swapping treatments on plots s1 and s2
                  dswap=(1+M12[Trts[s1],BF[s2]]+M12[Trts[s2],BF[s1]]-M12[Trts[s1],BF[s1]]-M12[Trts[s2],BF[s2]])**2-
                    (2*M11[Trts[s1],Trts[s2]]-M11[Trts[s1],Trts[s1]]-M11[Trts[s2],Trts[s2]])*(2*M22[BF[s1],BF[s2]]-M22[BF[s1],BF[s1]]-M22[BF[s2],BF[s2]])
                }
                #updates matrices
                prop_change=prop_change*dswap
                up=UpDate(M11,M22,M12,s1,s2,Trts,BF)	
                M11=up$M11
                M22=up$M22
                M12=up$M12
                Trts[c(s1,s2)]=Trts[c(s2,s1)]	
              } 
              locrelD=locrelD*prop_change
            } 
          } # for r
          Trts=globTF
        } #if count
      } # if class			
    }	# for ortho
  } # if ortho
  
  # add back single replicate treatments here 
  if ( !all(replicates>1) & !all(replicates==1) ) {
    nunits=sum(treatments*replicates)	
    fullblocksizes=nunits
    for (i in 2 : strata) {
      subSizes=vector(length=(length(fullblocksizes)*blocklevels[i]))
      for (z in 1: length(fullblocksizes)) {
        bsize = fullblocksizes[z] %/% blocklevels[i]
        bresid = fullblocksizes[z] %% blocklevels[i]
        for (w in 1 : blocklevels[i]) {
          subSizes[(z-1)*blocklevels[i]+w] = bsize + (bresid>0)
          bresid=bresid-1
        }
      }
      fullblocksizes=subSizes
    }
    # append and reorder full Trts in full blocks
    Trts=c(Trts,sample(c((ntrts+1):(ntrts+sum(treatments[replicates==1])))))
    Trts=Trts[order(c( rep( 1:length(blocksizes),blocksizes),  rep( 1:length(blocksizes),(fullblocksizes-blocksizes) ) ) )]
    # single replicates are appended at the end of the new treatments but may need re-ordering according to the treatments list order
    newlabels=c( rep(1:sum(treatments))[ rep((replicates>1),treatments) ],rep(1:sum(treatments))[rep((replicates==1),treatments)])
    for (i in 1 :length(Trts)) 
      Trts[i]=newlabels[Trts[i]]	
    
    desMat= matrix(1,nrow=nunits,ncol=strata)
    for (r in 1 : strata) 
      desMat[,r]=rep(facMat[,r],fullblocksizes)
    blocksizes=fullblocksizes	
    ortho=0
    ntrts=sum(treatments)	   
  }	
  rand=sample(1:nunits)
  Trts=Trts[rand][order(desMat[,strata][rand])]	
  Design=as.data.frame(cbind(desMat[,c(2:ncol(desMat))],Trts))
  Design[]=lapply(Design, factor)	
  designnames="Main"
  if (strata>2)
    for (i in 1:(strata-2))
      designnames=c(designnames,paste("Sub",i))
  plannames=c("Blocks",designnames)
  designnames=c(designnames,"Treatments")
  colnames(Design)=designnames
  
  Incidences=vector(mode = "list", length =(strata-1) )
  for (i in 1:(strata-1))
    Incidences[[i]]=table(Design[,c(i,ncol(Design))])	
  maxb=max(blocksizes)
  plots=matrix(nrow=length(blocksizes),ncol=maxb)
  count=1
  for (i in 1:length(blocksizes)) {
    for (j in 1:blocksizes[i]) {
      plots[i,j]=Trts[count]
      count=count+1
    }
    if (blocksizes[i]<maxb) plots[i,maxb]="" 
  }	
  blank=rep("",length(blocksizes))
  plannames=c(plannames,"Treatments")
  for (i in 1:max(blocksizes)) plannames=c(plannames,paste("Plot",i))	
  Plan=as.data.frame(cbind(blank,facMat[,c(2:ncol(facMat))],blank,plots))
  colnames(Plan)=plannames
  
  bounds=rep(0,(strata-1))
  r=1/sqrt(tabulate(Trts))
  
  aeff=c(rep(0,(strata-1)))
  for (i in 1:(strata-1)) {	
    k=1/sqrt(tabulate(Design[,i]))
    U=crossprod(t(crossprod(diag(r,nrow = length(r)),table(Trts,Design[,i]))),diag(k,nrow = length(k)))
    A=diag(length(r))-crossprod(t(U))
    aeff[i]=1/mean(1/eigen(A, symmetric=TRUE, only.values = TRUE)$values[1:length(r)-1])
    if (max(replicates)==min(replicates))
      bounds[i]=upper_bounds(nunits,ntrts,cumblocklevs[i+1]) 
    else
      bounds[i]=NA
  }
  effics=as.data.frame(cbind( cumblocklevs[2:strata], aeff, bounds))
  rnames=c("Main")
  if (strata>2)
    for (i in 1 : (strata-2)) rnames=c(rnames,paste("Sub",i))
  colnames(effics)=c("Blocks","A-Efficiencies", "Upper Bounds")
  rownames(effics)=rnames
  list(Design=Design,Plan=Plan,Incidences=Incidences,Efficiencies=effics,seed=seed)
} 
