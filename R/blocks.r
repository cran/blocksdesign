#' @title Block designs for unstructured treatment sets
#'
#' @description
#'
#' Constructs randomized nested block designs for unstructured treatment sets.
#'
#' @details
#'
#' Constructs randomized nested block designs for any arbitrary
#' number of unstructured treatments and any arbitrary feasible depth of nesting.
#' 
#' \code{treatments} is a set of numbers that partitions the total number of treatments into
#' equi-replicate treatment sets.
#' 
#' \code{replicates} is a set of treatment replication numbers that defines the replication number of each equi-replicate
#' \code{treatments} set.
#' 
#' \code{blocks} is a sequence of nested block levels in descending order of block sizes from largest to smallest
#' where each block level is the number of blocks nested within each preceding block
#' (assuming a single top-level super-block). The block sizes will automatically be as nearly equal as
#'  possible and will never differ in size by more than a single plot for any particular level of nesting. 
#'
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number 
#' of unreplicated treatments using the \code{treatments} and \code{replicates} formulae. 
#' However, it will usually be preferable to find an efficient block design
#' for just the replicated treatment sets and then to add the unreplicated treatments individually by hand. 
#'
#' Resolvable incomplete block designs for r replicates of p*p treatments arranged in blocks of size p where
#' r < p+2 for prime or prime power p or r < 4 for general p are square lattices and 
#' these designs are constructed automatically from Latin squares or MOLS.
#'
#' Resolvable incomplete block designs for r replicates of (p-1)*p treatments arranged in blocks of size p-1
#' where r < p+1 for prime or prime power p are rectangular lattices and these designs are constructed automatically
#' by reducing an appropriate algebraic square lattice, see Cochran and Cox, Experimental Designs, 2nd Edition, 
#' Page 417 (Shrikhande method).
#'
#' Outputs:
#'
#' \itemize{
#' \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order.\cr
#' \item  A table showing the replication number of each treatment in the design. \cr
#' \item  A table showing the block levels and the achieved D-efficiency and A-efficiency factor for each nested level together
#'   with A-efficiency upper bounds, where available. \cr
#' \item  A plan showing the allocation of treatments to blocks in the bottom level of the design.\cr
#' }
#'
#' @param treatments  the required number of treatments partitioned into equally replicated treatment sets.
#'
#' @param replicates  the replication number for each partitioned treatment set.
#'
#' @param blocks the number of nested blocks in each level of nesting from the top level down.
#' 
#' @param seed an integer initializing the random number generator.
#'
#' @param searches the maximum number of local optima searched for a design optimization. 
#'
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima.
#'
#' @return
#' \item{Treatments}{A table showing the replication number of each treatment in the design.}
#' \item{Design}{Data frame giving the optimized block and treatment design in plot order.}
#' \item{Plan}{Data frame showing a plan view of the treatment design in the bottom level of the design.}
#' \item{Blocks_model}{The D-efficiencies and the A-efficiencies of the blocks in each nested level of the 
#'  design together with A-efficiency upper-bounds, where available.}
#' \item{seed}{Numerical seed used for random number generator.}
#' \item{searches}{Maximum number of searches used for each level.}
#' \item{jumps}{Number of random treatment swaps used to escape a local maxima.}
#'
#' @references
#'
#' Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
#'
#' @examples
#' 
#' ## The number of searches in the following examples have been limited for fast execution.  
#' ## In practice, the number of searches may need to be increased for optimum results.
#' ## Designs should be rebuilt several times to check that a near-optimum design has been found.  
#' 
#' # 12 treatments with 4 replicates and 1 control treatment with 8 replicates 
#' # the blocks automatically default to 4 complete randomized blocks each of size 14
#' blocks(treatments=list(12,1),replicates=list(4,8))
#' 
#' # 12 treatments x 4 replicates in 4 complete blocks with 4 sub-blocks of size 3
#' # rectangular lattice see Plan 10.10 Cochran and Cox 1957.
#' \donttest{blocks(treatments=12,replicates=4,blocks=list(4,4))}
#'
#' # 3 treatments x 2 replicates + 2 treatments x 4 replicates in two complete randomized blocks
#' blocks(treatments=list(3,2),replicates=list(2,4),blocks=2,searches=10)
#'
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block
#' blocks(treatments=50,replicates=4,blocks=list(4,5))
#'
#' # as above but with 20 additional single replicate treatments, one single treatment per sub-block
#' \donttest{blocks(treatments=list(50,20),replicates=list(4,1),blocks=list(4,5))}
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,blocks=4)
#'
#' # 128 treatments x 2 replicates with two main blocks and 3 levels of nesting
#'  \donttest{blocks(128,2,list(2,2,2,2))}
#' 
#' #' # 64 treatments x 4 replicates with 4 main blocks nested blocks of size 8 (lattice square)
#' blocks(64,4,list(4,8)) 
#' 
#' # 100 treatments x 4 replicates with 4 main blocks nested blocks of size 10 (lattice square)
#' blocks(100,4,list(4,10)) 
#' 
#' @export
#' @importFrom stats coef anova lm model.matrix as.formula setNames 
#'
 blocks = function(treatments,replicates,blocks=NULL,searches=NULL,seed=NULL,jumps=1) {
   options(contrasts=c('contr.SAS','contr.poly'))
   options(warn=0)
   tol = .Machine$double.eps ^ 0.5
   if (is.list(treatments)) treatments=unlist(treatments)
   if (is.list(blocks)) blocks=unlist(blocks)
   if (is.list(replicates)) replicates=unlist(replicates)
 # ***********************************************************************************************
 # Tests for and constructs square lattice designs in top 2 levels of a lattice design.
 #  Further nested levels for levels 3... etc can be nested within the levels of the lattice blocks
 #  Allocates v*v treatments to blocks assuming a lattice design with r complete replicate blocks
 #  and nested blocks of size v. Returns a simple lattice for r=2 or a triple lattice for 
 #  r=3 for any size of v. Returns a lattice for any r < v + 2 if v is prime or any prime power 
 # where primes=c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97) and
 #  powers=c(12,7,5,4,3,3,3,3,rep(2,17)). Returns a triple lattice for v=10 and r=4. 
 # ***********************************************************************************************
   squarelattice=function(v,u) {
     PP=isPrimePower(v)
     p=PP$base
     q=PP$power
     primes=c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97) 
     powers=c(12,7,5,4,3,3,3,3,rep(2,17))
     ppLat = ( p%in%primes)
     if (ppLat) qmax=powers[which(primes == p)] else qmax=0
     
     if (ppLat & q>1 & q<=qmax & u>1)   {
       mols=MOLS(p,q,u)    #  lattice designs for prime powers where q > 1 and r > 3  
     } else if ( (isPrime(v) & u<v) | u < 2) {
       # cyclic designs for any prime v with r < v+2 or any v with r<4 
       mols=sapply(0:u,function(z){ sapply(0:(v-1), function(j){ (rep(0:(v-1))*z +j)%%v}) })
       mols=data.frame(row=rep(1:v,v),mols+1) 
     } else if (v==10 & u<3 & u>0) {
       #  lattice design for v = 10 and r =3 or r= 4 
       square1=       c(1,    8,    9,    2,    0,    5,    3,    6,    4,     7,
                        9,    1,    0,    4,    2,    7,    8,    5,    3,     6,
                        0,    3,    1,    6,    8,    9,    7,    4,    2,     5,
                        3,    7,    4,    1,    5,    2,    6,    9,    8,     0,
                        8,    9,    5,    0,    1,    6,    4,    2,    7,     3,
                        2,    6,    3,    7,    4,    1,    5,    8,    0,     9,
                        5,    2,    6,    3,    7,    4,    1,    0,    9,     8,
                        4,    0,    7,    5,    3,    8,    9,    1,    6,     2,
                        7,    5,    8,    9,    6,    0,    2,    3,    1,     4,
                        6,    4,    2,    8,    9,    3,    0,    7,    5,     1)
       square2 =      c(1,    3,    5,    4,    2,    6,    7,    8,    9,     0,
                        3,    4,    8,    6,    7,    9,    1,    5,    0,     2,
                        5,    6,    7,    0,    9,    1,    4,    3,    2,     8,
                        9,    8,    1,    2,    3,    0,    5,    6,    7,     4,
                        2,    0,    4,    1,    6,    7,    8,    9,    5,     3,
                        8,    1,    2,    3,    0,    5,    9,    4,    6,     7,
                        0,    5,    9,    8,    1,    2,    3,    7,    4,     6,
                        4,    9,    6,    7,    5,    8,    2,    0,    3,     1,
                        7,    2,    0,    9,    4,    3,    6,    1,    8,     5,
                        6,    7,    3,    5,    8,    4,    0,    2,    1,     9)
       
       mols=data.frame(row=rep(1:v,v),col=rep(1:v,each=v),square1+1)
       if (u==2) mols=cbind(mols,square2+1)
     } else return(NULL)
     
     TF=factor(sapply(1:(u+2),function(i){seq_len(v*v)[order(as.numeric(mols[[i]]))]}))
     TF=data.frame(Reps=rep(sample(1:(u+2)),each=(v*v)),Blocks=rep(sample(1:(v*(u+2))),each=v),plots=sample(1:(v*v*(u+2))),TF)
     TF=TF[ do.call(order, TF), ]
     names(TF)[4] = "treatments"
     rownames(TF)=NULL
     return(TF[,4,drop=FALSE])
   }
   # **************************************************************************************************
  # Tests for and constructs rectangularlattice designsin top 2 levels of a rectangular lattice.
  # Further nested levels for levels 3... etc can be nested within the levels of the rectangular lattice blocks
  # ***************************************************************************************************
   rectlattice=function(v,r) {
    LT=squarelattice(v,r-1)
    if (is.null(LT)) return(NULL)
    
    LT=split(LT[,1], factor(rep(1:((r+1)*v), each=v))) 
    drop=factor((v*(v-1)+1) : (v*v))
    dropblock=which(sapply (1:length(LT), function(i) all(drop%in%LT[[i]]))   )
    droprep=(dropblock-1)%/%v + 1
    omitblocks=((droprep-1)*v + 1):(droprep*v)
    LT[omitblocks]=NULL
    LT=unlist(LT)
    TF=data.frame(droplevels(LT[!LT%in%drop]))
    names(TF)[1] = "treatments"
    rownames(TF)=NULL
    return(TF)
   }
  # ***************************************************************************************************
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the 
  # best swap for D-efficiency improvement. Sampling is used initially when many feasible swaps are 
  # available but later a full search is used to ensure steepest ascent optimization.
  # ***************************************************************************************************
  DMax=function(VTT,VBB,VTB,TF,MF,BF,TM,BM) { 
    fBF=interaction(MF,BF, lex.order = TRUE)
    locrelD=1
    mainSizes=tabulate(MF)
    nSamp=pmin(rep(8,nlevels(MF)), mainSizes)
    repeat {
      kmax=1
      for (k in 1: nlevels(MF)) {
        s=sort(sample(seq_len(nrow(TF))[MF==levels(MF)[k]], nSamp[k])) 
        TB=VTB[TF[s,],fBF[s],drop=FALSE]
        TT=VTT[TF[s,],TF[s,],drop=FALSE]
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
        up=UpDate(VTT,VBB,VTB,TF[pi,],TF[pj,],fBF[pi],fBF[pj])
        VTT=up$VTT
        VBB=up$VBB
        VTB=up$VTB
        TF[c(pi,pj),]=TF[c(pj,pi),]
        TM[c(pi,pj),]=TM[c(pj,pi),]
      } else if (sum(nSamp) == nrow(TF)) break
      else nSamp=pmin(mainSizes,2*nSamp)
    }
    list(VTT=VTT,VBB=VBB,VTB=VTB,TF=TF,locrelD=locrelD,TM=TM)
  } 
  
  # *************************************************************************************************
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb   
  # ***************************************************************************************************
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
  # *******************************************************************************************************************
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # nested blocks only
  # ********************************************************************************************************************
  Optimise=function(TF,MF,BF,VTT,VBB,VTB,TM,BM) {
    fBF=interaction(MF,BF, lex.order = TRUE)
    globrelD=0
    relD=1
    globTF=TF
    blocksizes=tabulate(fBF)
    if (regReps & max(blocksizes)==min(blocksizes) )
      blocksEffBound=A_bound(n=nrow(TF),v=nlevels(TF),b=nlevels(fBF)) 
    else blocksEffBound=1
    for (r in 1:searches) {
      dmax=DMax(VTT,VBB,VTB,TF,MF,BF,TM,BM)
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
            reff=blockEstEffics(globTF,fBF)[2]
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
            available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
            z= seq_len(nunits)[MF==MF[s1] & fBF!=fBF[s1] & available ] 
            if (length(z)==0) next
            if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z)
            v11=VTT[c(TF[s[1],],TF[s[2],]),c(TF[s[1],],TF[s[2],])]
            v22=VBB[c(fBF[s[1]],fBF[s[2]]),c(fBF[s[1]],fBF[s[2]])]
            v12=VTB[c(TF[s[1],],TF[s[2],]),c(fBF[s[1]],fBF[s[2]])]
            Dswap=(1-2*sum(diag(v12))+sum(v12))**2-(2*sum(diag(v22))-sum(v22))*(2*sum(diag(v11))-sum(v11))
            if (Dswap>.1| counter>1000) break
          }
          if (counter>1000) return(globTF) # finish with no non-singular swaps
          relD=relD*Dswap 
          up=UpDate(VTT,VBB,VTB,TF[s[1],],TF[s[2],], fBF[s[1]], fBF[s[2]])
          VTT=up$VTT
          VBB=up$VBB
          VTB=up$VTB
          TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]  
          TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
        } #end of jumps
      }
    }
    return(globTF)
  }
  # ********************************************************************************************************************
  # Random swaps
  # ********************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank) {
    candidates=NULL
    while (is.null(candidates)) {
      if (rank<(nrow(TF)-1)) 
        s1=sample(pivot[(1+rank):nrow(TF)],1) else 
          s1=pivot[nrow(TF)]
        available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
        candidates = (1:nrow(TF))[ MF==MF[s1] & BF!=BF[s1] & available ] 
    }
    if ( length(candidates)>1) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  } 
  # *******************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with 
  # positive selection are used to to increase design rank
  # ********************************************************************************************************************
  NonSingular=function(TF,BF,TM,BM,restrict) {
    fullrank=ncol(cbind(TM,BM))
    Q=qr(t(cbind(TM,BM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:100) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      swap=Swaps(TF,restrict,BF,pivot,rank)
      TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
      Q=qr(t(cbind(TM,BM)))
      if (Q$rank>rank) {
        TF[c(swap[1],swap[2]),]=TF[c(swap[2],swap[1]),]
        rank=Q$rank
        pivot=Q$pivot
      } else TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
    }
    stop(" Unable to find an initial non-singular choice of treatment design for this choice of block design")
  }
  # *****************************************************************************************************************
  # Optimizes the nested Blocks assuming a possible set of Main block constraints and assuming a randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # The model matrix for the full set of treatment or block effects is represented here by SAS contrast which assume an overall 
  # mean equal to the last treatment or block level which is omittted to give a full rank model. 
  # The variance matrices are calculated by taking contrasts WITHIN each term to give the contrast of each term from the last term.
  # The first term of each model is a zero vector which is moved to the last term to represent the deviations of the last term from itself. 
  # The covariances of the non null terms and the variances and covariances of the last term are invariant and are equated to zero. 
  # The model.matrix function MUST define SAS type contrasts and the options setting for this function are SAS contrasts.
  # For nested blocks, the full set of nested block contrasts including the final column of zeroes
  # must be nested within each main block which is achieved using the reorder function shown below 
  # *****************************************************************************************************************
  blocksOpt=function(TF,MF,BF) {
      fBF=interaction(MF,BF, lex.order = TRUE)
      TM=model.matrix(as.formula("~treatments"),TF)[,-1,drop=FALSE] # drops mean contrast
      TM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)}))
      BM=model.matrix(as.formula(~fBF))[,-1,drop=FALSE]
      BM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(BM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)}))
      BM=BM[,-c((1:nlevels(MF))*nlevels(BF)) ,drop=FALSE]
      if ((ncol(cbind(TM,BM))+1) > nrow(TM)) stop( paste("Too many parameters: plots =",nrow(TM)," parameters = ",ncol(cbind(TM,BM)))) 
      nonsing=NonSingular(TF,BF,TM,BM,MF)
      TF=nonsing$TF
      TM=nonsing$TM
      V=chol2inv(chol(crossprod(cbind(TM,BM))))
      VTT = matrix(0, nrow=nlevels(TF[,1]), ncol=nlevels(TF[,1]) )
      VBB = matrix(0, nrow=nlevels(fBF),ncol=nlevels(fBF))
      VTB = matrix(0, nrow=nlevels(TF[,1]), ncol=nlevels(fBF) )
      VTT[1:ncol(TM),1:ncol(TM)] = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
      VBB[1:ncol(BM),1:ncol(BM)] = V[(ncol(TM)+1):ncol(V),(ncol(TM)+1):ncol(V),drop=FALSE]
      VTB[1:ncol(TM),1:ncol(BM)] = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
      reorder=as.vector(rbind( matrix(1:(nlevels(fBF)-nlevels(MF)),nrow=(nlevels(BF)-1),ncol=nlevels(MF)), nlevels(fBF)-nlevels(MF) + seq(1,nlevels(MF))))
      VBB=VBB[reorder,reorder] 
      VTB=VTB[,reorder] 
      TF=Optimise(TF,MF,BF,VTT,VBB,VTB,TM,BM)
    return(TF)
  }
  # ***********************************************************************************************
  # Block design efficiency factors
  # *********************************************************************************************** 
  blockEstEffics=function(TF,BF) {
    if (nlevels(BF)==1) return(c(1,1))
    if (is.data.frame(TF)) TF=TF[,1]
    if (is.data.frame(BF)) BF=BF[,1]
    TM  = scale(model.matrix(~ TF))[,-1]
    BMO = scale(model.matrix(~ BF))[,-1]
    TMO = qr.Q(qr(TM)) # orthogonal basis for TM 
    BMO = qr.Q(qr(BMO)) # orthogonal basis for BM
    if (nlevels(TF)<=nlevels(BF)) 
      E=eigen(diag(ncol(TMO))-tcrossprod(crossprod(TMO,BMO)),symmetric=TRUE,only.values = TRUE)
    else E=eigen(diag(ncol(BMO))-tcrossprod(crossprod(BMO,TMO)),symmetric=TRUE,only.values = TRUE)
    Deff=prod(E$values)**(1/ncol(TMO))
    if (nlevels(TF)<=nlevels(BF)) Aeff=ncol(TMO)/sum(1/E$values)
    else Aeff=ncol(TMO)/(ncol(TMO)-ncol(BMO) + sum(1/E$values) )
    return(list(Deffic= round(Deff,5), Aeffic=round(Aeff,5)))
  }
  # *******************************************************************************************************************************
  # Finds efficiency factors for block designs
  # *******************************************************************************************************************************
  BlockEfficiencies=function(blkDesign,TF,orthogSize,blocks,hcf) {
    
    Design=data.frame(lapply(2:ncol(blkDesign),function(i) {interaction(blkDesign[,1:i], lex.order = TRUE)}),plots=factor(1:nunits),TF)
    sizes=lapply(1:(ncol(Design)-2),function(i) {tabulate(Design[,i])}) 
    effics=matrix(NA,nrow=(ncol(Design)-2),ncol=2)
    
    for (i in 1:(ncol(Design)-2))  
      if (hcf%%prod(blocks[1:i])==0 )  
        effics[i,]=1  else 
        effics[i,]=unlist(blockEstEffics(Design[,ncol(Design)],Design[,i]))
    
    sizes=lapply(1:(ncol(Design)-2),function(i) {tabulate(Design[,i])}) 
    bounds=sapply(1:(ncol(Design)-2),function(i){( all(floor(sizes[[i]]/orthogSize)==sizes[[i]]/orthogSize))})
    bounds=sapply(1:length(bounds),  function(i){ if(bounds[i]==FALSE) bounds[i]=NA else bounds[i]=1} )
    if (regReps)
      for (i in 1:(ncol(Design)-2))  {
        if (nunits%%nlevels(Design[,i])==0 )
          bounds[i]=A_bound(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,i]) )
      }
    blocklevs=unlist(lapply(1:(ncol(Design)-2), function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(labels(blocks),blocklevs,effics,round(bounds,5)))
    colnames(efficiencies)=c("Level","Blocks","D-Efficiency","A-Efficiency", "A-Bound")
    return(efficiencies)
  }
  # *******************************************************************************************************************************
  # Finds nested block sizes where all block sizes are as equal as possible within each level of nesting 
  # *******************************************************************************************************************************
  BlockSizes=function(blocks,nunits) {
  blocksizes=nunits
  for (i in 1:length(blocks))
    blocksizes=unlist(lapply(1:length(blocksizes), function(j) {
      nestblocksizes=rep(blocksizes[j]%/%blocks[i],blocks[i])
      blocksizes=nestblocksizes+rep(c(1,0),c(blocksizes[j]-sum(nestblocksizes), blocks[i]-blocksizes[j]+sum(nestblocksizes)))
    }))
  blocksizes
  }
  # *******************************************************************************************************************************
  #  Expands set of 'blocks' factor levels in left to right 'natural' order with first set changing slowest and last changing fastest
  # *******************************************************************************************************************************
  blocksGrid=function(blocks){expand.grid(lapply(length(blocks):1,function(i) {factor(seq(blocks[i]),                                                             
                              labels=lapply(1:blocks[i], function(j){paste0("Blocks_",j)}))}))[length(blocks):1]}
  
  # *************************************************************************************************************************
  # Main design program which tests input variables, omits any single replicate treatments, optimizes design, replaces single
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # *************************************************************************************************************************
  if (missing(treatments)|is.null(treatments)) stop(" Treatments missing or not defined ")
  if (missing(replicates)|is.null(replicates)) stop(" Replicates missing or not defined ")
  if (!is.null(seed)) set.seed(seed)
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>25) 
    stop(" maximum number of jumps is 25 ")
  if (any(replicates%%1!=0)|any(replicates<1)) stop(" replication numbers must be integers")
  if (anyNA(treatments)|any(is.nan(treatments))|any(!is.finite(treatments))|any(treatments<1)) 
    stop(" treatments parameter invalid")
  if (length(replicates)!=length(treatments)) 
    stop("treatments and replicates parameters must both be the same length")
  hcf=HCF(replicates) 
  if (is.null(blocks)) blocks=hcf
  if (anyNA(blocks)|any(is.nan(blocks))|any(!is.finite(blocks))|any(blocks%%1!=0)|any(blocks<1)|is.null(blocks)) 
    stop(" blocks invalid")
  blocks=blocks [!blocks %in% 1]
  if (length(blocks)==0) blocks=1
  ntrts=sum(treatments)
  TF=data.frame(unlist(lapply(1:hcf,function(i){sample(rep(factor(1:ntrts), rep(replicates/hcf,treatments)))})))
  colnames(TF)="treatments"
  nunits=nrow(TF)
  if (is.null(searches)) searches=1+5000%/%nunits
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  if (is.null(names(blocks))) 
    names(blocks)=unlist(lapply(1:length(blocks), function(j) {paste0("Level_",j-1)}))
  if (prod(blocks)*2>nunits) stop("Too many blocks for the available plots  - each block must contain at least two plots")
  blkdesign=blocksGrid(blocks)
  colnames(blkdesign)=labels(blocks)
  blkdesign=cbind("Null"=factor(rep("Blocks_1",nrow(blkdesign))),  blkdesign)
  blocksizes=BlockSizes(blocks,nunits)
  blkDesign=blkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  orthogSize=sum(treatments*replicates)/hcf
  
  if (length(blocks)>1) regBlocks = (floor(nunits/prod(blocks[1:2]))==nunits/prod(blocks[1:2])) else 
    regBlocks = floor(nunits/blocks)==nunits/blocks
  regReps = isTRUE(all.equal(max(replicates), min(replicates)))
  
  if (length(blocks)>1) 
    k=nunits/prod(blocks[1:2]) else k=nunits/blocks # average block size
  if (length(blocks)>1) 
    b=prod(blocks[1:2]) else b=blocks # number of blocks
  
  sk=sqrt(ntrts)  # size of blocks of a square lattice 
  rk=(sqrt(1+4*ntrts)-1)/2 #  size of blocks of a rectangular lattice 
  
  Lattice=isTRUE(regBlocks & regReps & identical(sk,floor(sk)) & replicates[1]<(sk+2) & k==sk & 
             identical(floor(nunits/blocks[1]),nunits/blocks[1]))
  rectLattice=isTRUE(regBlocks & regReps & identical(rk,floor(rk)) & replicates[1]<(rk+2) & k==rk & 
                 identical(floor(nunits/blocks[1]),nunits/blocks[1]))
  if (Lattice) {
    lsTF=squarelattice(sk,(replicates[1]-2)) 
    if (length(blocks)>2 & !is.null(lsTF))
      for (i in 3:length(blocks)) 
        lsTF=blocksOpt(lsTF,interaction(blkDesign[1:i], lex.order = TRUE),blkDesign[,i+1]) 
  } else if (rectLattice) {
    lsTF=rectlattice(rk+1,replicates[1])
    if (length(blocks)>2 & !is.null(lsTF))
      for (i in 3:length(blocks)) 
        lsTF=blocksOpt(lsTF,interaction(blkDesign[1:i], lex.order = TRUE),blkDesign[,i+1]) 
  } else lsTF=NULL
  
  if (!is.null(lsTF)) TF=lsTF else if (is.null(lsTF)) {
      for (i in 1:length(blocks)) {
        if (hcf%%prod(blocks[1:i])!=0 ) { 
           TF=blocksOpt(TF,interaction(blkDesign[1:i], lex.order = TRUE),blkDesign[,i+1]) 
        }
      }
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - try a simpler block or treatment design")
  Efficiencies=BlockEfficiencies(blkDesign,TF,orthogSize,blocks,hcf)
  Design=data.frame(blkDesign,plots=factor(1:nunits),treatments=TF)[,-1,drop=FALSE]
  V = split( Design[,ncol(Design)], interaction(Design[,1:(ncol(Design)-2)],lex.order = TRUE)) 
  V = lapply(V, function(x){ length(x) =max(blocksizes); x })
  Plan = data.frame(blkdesign[,-1 ,drop=FALSE],rep("",length(V)),matrix(unlist(V),nrow=length(V),byrow=TRUE))
  
  colnames(Plan)=c(colnames(blkdesign[,-1 ,drop=FALSE]),"Blocks.Plots:", c(1:max(blocksizes)))
  row.names(Plan)=NULL
  row.names(Design)=NULL
  row.names(Efficiencies)=NULL
  Treatments=count(Design[,ncol(Design),drop=FALSE])

  list(Treatments=Treatments,Blocks_model=Efficiencies,Design=Design,Plan=Plan,seed=seed,searches=searches,jumps=jumps)
}