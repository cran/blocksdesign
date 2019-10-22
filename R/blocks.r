#' @title Block designs for unstructured treatment sets
#'
#' @description
#'
#' Constructs randomized nested block designs for unstructured treatment sets and any feasible depth of nesting.
#'
#' @details
#'
#' Constructs randomized nested block designs with arbitrary depth of nesting for arbitrary unstructured treatment sets.
#' 
#' The \code{treatments} parameter is a set of numbers that gives a partition of the total number of treatments
#'  while the \code{replicates} parameter is a matching set of numbers that defines the 
#' replication of each treatment set in the partition.
#' 
#' The \code{blocks} parameter, if any, defines the number of blocks for each level of nesting from the highest
#' to the lowest. The first number, if any, is the number of nested row blocks in the first-level of nesting, 
#' the second number, if any, is the number of nested row blocks in
#' the second-level of nesting and so on down to any required feasible depth of nesting.
#'
#' Block sizes are as nearly equal as possible and will never differ by more than a single plot in any 
#' particular block classification. 
#'
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number 
#' of single unreplicated treatments to give augmented blocks that never differ in size by more than a single plot.
#' However, it may sometimes be preferable to find an efficient 
#' block design for the replicated treatments and then add the unreplicated treatments to the design heuristically. 
#'
#' Square lattice designs are resolvable incomplete block designs for r replicates of p*p treatments 
#' arranged in blocks of size p where r < p+2 for prime or prime power p or r < 4 for general p. 
#' Square lattice designs are constructed algebraically from Latin squares or MOLS.
#'
#' Rectangular lattice designs are resolvable incomplete block designs for r replicates of (p-1)*p treatments 
#' arranged in blocks of size p-1 where r < p+1 for prime or prime power p. Rectangular lattice designs are
#' constructed algebraically by reducing an algebraic square lattice, see 
#' Cochran and Cox, Experimental Designs. 2nd Edition. Page 417 (Shrikhande method).
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
#' @param treatments  a partition of the total required number of treatments into equally replicated treatment sets, possibly a complete
#' partition into individual treatments.
#'
#' @param replicates  a set of treatment replication numbers with one replication number for each partitioned treatment set, possibly a 
#' complete set of treatment replication numbers.
#'
#' @param blocks the number of blocks nested in each preceding block for each level of nesting from the top-level block downwards.
#' 
#' @param seed an integer initializing the random number generator. The default is a random seed.
#'
#' @param searches the maximum number of local optima searched for a design optimization. The default number decreases
#' as the design size increases.
#'
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#'
#' @return
#' \item{Treatments}{A table showing the replication number of each treatment in the design.}
#' \item{Design}{Data frame giving the optimized block and treatment design in plot order.}
#' \item{Plan}{Data frame showing a plan view of the treatment design in the bottom level of the design.}
#' \item{blocks_model}{The D-efficiencies and the A-efficiencies of the blocks in each nested level of the 
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
#' # 12 treatments x 4 replicates in 4 complete blocks with 4 sub-blocks of size 3
#' # rectangular lattice see Plan 10.10 Cochran and Cox 1957.
#' \donttest{blocks(treatments=12,replicates=4,blocks=c(4,4))}
#'
#' # 3 treatments x 2 replicates + 2 treatments x 4 replicates in two complete randomized blocks
#' blocks(treatments=c(3,2),replicates=c(2,4),searches=10)
#'
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block
#' blocks(treatments=50,replicates=4,blocks=c(4,5))
#'
#' # as above but with 20 additional single replicate treatments, one single treatment per sub-block
#' \donttest{blocks(treatments=c(50,20),replicates=c(4,1),blocks=c(4,5))}
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,blocks=4)
#'
#' # 128 treatments x 2 replicates with two main blocks and 3 levels of nesting
#'  \donttest{blocks(128,2,c(2,2,2,2))}
#' 
#' #' # 64 treatments x 4 replicates with 4 main blocks nested blocks of size 8 (lattice square)
#' blocks(64,4,c(4,8)) 
#' 
#' # 100 treatments x 4 replicates with 4 main blocks nested blocks of size 10 (lattice square)
#' blocks(100,4,c(4,10)) 
#' 
#' @export
#' @importFrom stats coef anova lm model.matrix as.formula setNames 
#'
 blocks = function(treatments,replicates,blocks=NULL,searches=NULL,seed=NULL,jumps=1) {
   
   
  # ***********************************************************************************************
  # Block design efficiency factors
  # *********************************************************************************************** 
   blockEstEffics=function(TF,BF) {
     if (nlevels(BF)==1) return(c(1,1))
     if (is.data.frame(TF))TF=TF[,1]
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

     return(list(Deffic= round(Deff,7), Aeffic=round(Aeff,7)))
   }
   
 # ***********************************************************************************************
 # Tests for and constructs squarelattice designs
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
     df=data.frame(Reps=rep(sample(1:(u+2)),each=(v*v)),Blocks=rep(sample(1:(v*(u+2))),each=v),Plots=sample(1:(v*v*(u+2))),TF)
     df=df[ do.call(order, df), ]
     return(df[,4])
   }
   
   # **************************************************************************************************
  # Tests for and constructs rectangularlattice designs
  # ***************************************************************************************************
  rectlattice=function(v,r) {
    
    LT=squarelattice(v,r-1)
    if (is.null(LT)) return(NULL)
    LT=split(LT, factor(rep(1:((r+1)*v), each=v))) 
    drop=factor((v*(v-1)+1) : (v*v))
    dropblock=which(sapply (1:length(LT), function(i) all(drop%in%LT[[i]]))   )
    droprep=(dropblock-1)%/%v + 1
    omitblocks=((droprep-1)*v + 1):(droprep*v)
    LT[omitblocks]=NULL
    LT=unlist(LT)
    LT=(droplevels(LT[!LT%in%drop]))
    return(LT)
  }
  # ***************************************************************************************************
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the 
  # best swap for D-efficiency improvement. Sampling is used initially when many feasible swaps are 
  # available but later a full search is used to ensure steepest ascent optimization.
  # ***************************************************************************************************
     DMax=function(VTT,VBB,VTB,TF,MF,fBF,BF,TM) {  
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
    fBF=interaction(MF:BF)
    globrelD=0
    relD=1
    globTF=TF
    blocksizes=tabulate(fBF)
    if (regReps & max(blocksizes)==min(blocksizes) )
      blocksEffBound=A_bound(nrow(TF),nlevels(TF),nlevels(fBF)) 
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
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # *****************************************************************************************************************
  blocksOpt=function(TF,MF,BF) {
      TM=model.matrix(as.formula("~Treatments"),TF)[,-1,drop=FALSE] # drops mean contrast
      TM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(TM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)}))
      if (nlevels(MF)==1) BM=model.matrix(as.formula(~BF))[,-1,drop=FALSE] else
        BM=model.matrix(as.formula(~MF+MF:BF))[,-c(1:nlevels(MF)),drop=FALSE]
      #centres sub-blocks within main blocks
      BM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(BM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)}))
      # reorder within main blocks
      BM=BM[ ,as.numeric(matrix(1:(nlevels(BF)*nlevels(MF)),nrow=nlevels(BF),ncol=nlevels(MF),byrow=TRUE)[-nlevels(BF),]),drop=FALSE]
      if ((ncol(cbind(TM,BM))+1) > nrow(TM)) stop( paste("Too many parameters: plots =",nrow(TM)," parameters = ",ncol(cbind(TM,BM)))) 
      nonsing=NonSingular(TF,BF,TM,BM,MF)
      TF=nonsing$TF
      TM=nonsing$TM
      V=chol2inv(chol(crossprod(cbind(TM,BM))))
      VTT = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
      VBB = V[ (ncol(TM)+1):ncol(V) , (ncol(TM)+1):ncol(V),drop=FALSE]
      VTB = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
      VTT=rbind ( cbind( VTT, rep(0,ncol(TM) )) , rep(0,(1+ncol(TM)) ))
      VBB=rbind(  cbind( VBB, matrix(0, nrow=ncol(BM),ncol=nlevels(MF)) ),matrix(0, nrow=nlevels(MF),ncol=(ncol(BM)+nlevels(MF))))
      VTB=rbind(  cbind( VTB, matrix(0, nrow=ncol(TM),ncol=nlevels(MF)) ),rep(0,(ncol(BM)+nlevels(MF)  )))
      reorder=c(rbind( matrix(seq_len(ncol(BM)), nrow=(nlevels(BF)-1), ncol=nlevels(MF)), seq(ncol(BM)+1, nlevels(BF)*nlevels(MF))))
      VBB=VBB[reorder,reorder] 
      VTB=VTB[,reorder] 
      TF=Optimise(TF,MF,BF,VTT,VBB,VTB,TM,BM)
    return(TF)
  }
  # *******************************************************************************************************************************
  # Finds efficiency factors for block designs
  # *******************************************************************************************************************************
  BlockEfficiencies=function(Design) {
    effics=matrix(NA,nrow=(ncol(Design)-2),ncol=2)
    for (i in seq_len(ncol(Design)-2)) 
      effics[i,]=unlist(blockEstEffics(Design[,ncol(Design)],Design[,i]))
    bounds=rep(NA,(ncol(Design)-2))
    if (regReps)
      for (i in seq_len(ncol(Design)-2))
        if (nunits%%nlevels(Design[,i])==0 )
          bounds[i]=A_bound(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,i]) )
    blocklevs=unlist(lapply(1:(ncol(Design)-2), function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(labels(blocks),blocklevs,effics,bounds))
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
  # *************************************************************************************************************************
  # Main design program which tests input variables, omits any single replicate treatments, optimizes design, replaces single
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # *************************************************************************************************************************
  options(contrasts=c('contr.SAS','contr.poly'))
  tol=.Machine$double.eps^0.5
  if (missing(treatments)|is.null(treatments)) stop(" Treatments missing or not defined ")
  if (missing(replicates)|is.null(replicates)) stop(" Replicates missing or not defined ")
  if (!is.null(seed)) set.seed(seed)
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) 
    stop(" maximum number of jumps is 10 ")
  if (any(replicates%%1!=0)|any(replicates<1)) stop(" replication numbers must be integers")
  if (anyNA(treatments)|any(is.nan(treatments))|any(!is.finite(treatments))|any(treatments<1)) 
    stop(" treatments parameter invalid")
  if (length(replicates)!=length(treatments)) 
    stop("treatments and replicates parameters must both be the same length")
  hcf=HCF(replicates) 
  if (is.null(blocks)) blocks=hcf
  if (anyNA(blocks)|any(is.nan(blocks))|any(!is.finite(blocks))|any(blocks%%1!=0)|any(blocks<1)|is.null(blocks)) 
    stop(" blocks invalid")
  ntrts=sum(treatments)
  Treatments=factor(1:ntrts)
  Treatments=unlist(lapply(1:hcf,function(i){sample(rep(Treatments, rep(replicates/hcf,treatments)))}))
  nunits=length(Treatments)
  if (is.null(searches)) 
    if (nunits<1000) searches=10000%/%nunits else if (nunits<5000) searches=5000%/%nunits else searches=1
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  
  if (is.null(names(blocks))) 
    names(blocks)=unlist(lapply(1:length(blocks), function(j) {paste0("Level_",j-1)}))
  
  if (prod(blocks)*2>nunits) stop("Too many blocks for the available plots  - each block must contain at least two plots")
  blocksizes=BlockSizes(blocks,nunits)
  
  
  grid=expand.grid(lapply(length(blocks):1,function(i) {factor(seq(blocks[i]),                                                             
                                                               labels=lapply(1:blocks[i], function(j){paste0("Blocks_",j)}))}))
  grid=grid[,length(blocks):1,drop=FALSE]
  colnames(grid)=labels(blocks)
  blkdesign=cbind("block0"=rep("Blocks_1",prod(blocks)),grid )
  
  
  blkDesign=blkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  TF=NULL
  regBlocks=isTRUE(all.equal(max(blocksizes), min(blocksizes)))
  regReps  =isTRUE(all.equal(max(replicates), min(replicates)))
  k=nunits/prod(blocks)  # average block size
  
  
  if (regReps & regBlocks & replicates[1]==blocks[1] & length(blocks)==2) { 
  s=sqrt(ntrts)  # dimension of a square lattice 
  Lattice=(identical(s,floor(s)) & replicates[1]<(s+2) & k==s)
  if (Lattice) TF=squarelattice(s,(replicates[1]-2)) 
  s=(sqrt(1+4*ntrts)+1)/2 # dimension of a rectangular lattice 
 rectLat=(identical(s,floor(s)) & replicates[1]<(s+1) & k==(s-1))
  if (rectLat) TF=rectlattice(s,replicates[1])
  }
  
  if (is.null(TF)) {
   if ( hcf%%prod(blocks)==0 ) 
     TF=data.frame(Treatments)
   else {
     attempts=0
     while (is.null(TF) & attempts<10) {
       attempts=attempts+1
       TF=data.frame(Treatments)
         for ( i in 1:length(blocks)) 
           TF=blocksOpt(TF,interaction(blkDesign[1:i]),blkDesign[,i+1]) 
     }
   }
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - try a simpler block or treatment design")
  
    NestDesign=data.frame(sapply(2:ncol(blkDesign),function(i) {interaction(blkDesign[,1:i])}),Plots=factor(1:nunits),TF)
    Efficiencies=BlockEfficiencies(NestDesign) 
    Design=data.frame(blkDesign,Plots=factor(1:nunits),Treatments=TF)[,-1,drop=FALSE]
    V = split(Design[,ncol(Design)],NestDesign[,(ncol(NestDesign)-2)])
    V = lapply(V, function(x){ length(x) =max(blocksizes); x })
    Plan = data.frame(blkdesign[,-1 ,drop=FALSE],rep("",length(V)),matrix(unlist(V),nrow=length(V),byrow=TRUE))
    colnames(Plan)=c(colnames(blkdesign[,-1 ,drop=FALSE]),"Blocks.Plots:", c(1:max(blocksizes)))
    row.names(Plan)=NULL
    row.names(Design)=NULL
    row.names(Efficiencies)=NULL
    Treatments=count(Design[,ncol(Design),drop=FALSE])
  list(Treatments=Treatments,blocks_model=Efficiencies,Design=Design,Plan=Plan,seed=seed,searches=searches,jumps=jumps)
}