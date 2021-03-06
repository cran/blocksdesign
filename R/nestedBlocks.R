#' @title nestedBlocks
#' @description
#' Internal function for optimizing nested block designs for a given set of nested block levels 
#' \code{BF} and a required set of treatments \code{TF} 
#' @param TF is a treatment factor with a treatment factor level for each plot level
#' @param BF is a data frame with a column for each set of nested blocks and a row for each plot 
#' @param searches is the number of optimizations searched
#' @keywords internal
#' 
nestedBlocks=function(TF,BF,searches,seed,jumps) {
  options(contrasts=c('contr.treatment','contr.poly'),warn=0)
  tol = .Machine$double.eps ^ 0.5
  
  # ***********************************************************************************************
  #  efficiency factors for TF and individual block factors BF
  # *********************************************************************************************** 
  blockEstEffics=function(TF,BF) {
      TM = qr.Q(qr(  scale(model.matrix(~ TF))[,-1]  )) # orthogonal basis for TM 
      BM = qr.Q(qr(  scale(model.matrix(~ BF))[,-1]  )) # orthogonal basis for BM
      if (nlevels(TF)<=nlevels(BF)) 
        E=eigen(diag(ncol(TM))-tcrossprod(crossprod(TM,BM)),symmetric=TRUE,only.values = TRUE) else
          E=eigen(diag(ncol(BM))-tcrossprod(crossprod(BM,TM)),symmetric=TRUE,only.values = TRUE)
      Deff=exp(sum(log(E$values))/ncol(TM))
      if (nlevels(TF)<=nlevels(BF)) Aeff=ncol(TM)/sum(1/E$values) else
         Aeff=ncol(TM)/(ncol(TM)-ncol(BM) + sum(1/E$values) )
    return(list(Deffic= round(Deff,7), Aeffic=round(Aeff,7)))
  }
  # *******************************************************************************************************************************
  # Finds efficiency factors for unstructured treatment factor TF and nested block designs BF
  # *******************************************************************************************************************************
  BlockEfficiencies=function(Design){
    TF=Design[,ncol(Design)]
    regreps=table(TF)
    regReps=isTRUE(max(regreps)==min(regreps))
    sizes = lapply(1:(ncol(Design)-2),function(i) {table(Design[,i])}) 
    regBlocks=sapply(1:length(sizes),function(i) {max(sizes[[i]])==min(sizes[[i]])})
    bounds=sapply(1:(ncol(Design)-2),function(i) {
      if (regBlocks[i] & regReps) A_bound(length(TF),nlevels(TF),nlevels(Design[,i])) else 1}) 
    
    blocklevs=unlist(lapply(1:(ncol(Design)-2), function(j) {nlevels(Design[,j])}))
    Effics  = t(sapply(1:(ncol(Design)-2),function(i) {
     if (nlevels(Design[,i])>1)
      blockEstEffics(TF,Design[,i]) else
        list(Deffic= 1, Aeffic=1)
      }) )
    efficiencies=data.frame(1:(ncol(Design)-2),blocklevs,as.numeric(Effics[,1]),as.numeric(Effics[,2]),round(bounds,7))
    colnames(efficiencies)=c("Level","Blocks","D-Efficiency","A-Efficiency", "A-Bound")
    return(efficiencies)
  }
  # ***************************************************************************************************
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the 
  # best swap for D-efficiency improvement. Sampling is used initially when many feasible swaps are 
  # available but later a full search is used to ensure steepest ascent optimization.
  # ***************************************************************************************************
  DMax=function(VTT,VBB,VTB,TF,MF,BF,TM,BM) { 
    locrelD=1
    mainSizes=tabulate(MF)
    nSamp=pmin(rep(8,nlevels(MF)), mainSizes)
    repeat {
      kmax=1
      for (k in 1: nlevels(MF)) {
        s=sort(sample(seq_len(length(TF))[MF==levels(MF)[k]], nSamp[k])) 
        TB=VTB[TF[s],BF[s],drop=FALSE]
        TT=VTT[TF[s],TF[s],drop=FALSE]
        BB=VBB[BF[s],BF[s],drop=FALSE]
        
        TT=sweep(sweep(2*TT,1,diag(TT)),2,diag(TT))
        BB=sweep(sweep(2*BB,1,diag(BB)),2,diag(BB))
        TBBT=sweep(sweep(TB+t(TB),1,diag(TB)),2,diag(TB))
        dMat=(1+TBBT)**2 - TT*BB
        
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nrow(dMat)
        j=1+(sampn-1)%/%nrow(dMat)
        if (dMat[i,j]>kmax) {kmax=dMat[i,j]; pi=s[i]; pj=s[j]} 
      }
      if (kmax>(1+tol)) {
        locrelD=locrelD*kmax
        up=UpDate(VTT,VBB,VTB,TF[pi],TF[pj],BF[pi],BF[pj])
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
  # *************************************************************************************************
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb 
  # UpDate applies for an unstructured treatment design TF where each row of TF is a single treatment factor
  # and a nested blocks design BF where each row of BF is a fully nested blocks factor
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
  Optimise=function(TF,MF,BF,VTT,VBB,VTB,TM,BM,regReps) {
    globrelD=0
    relD=1
    globTF=TF
    blocksizes=tabulate(BF)
    if (regReps & max(blocksizes)==min(blocksizes) )
      blocksEffBound=A_bound(n=length(TF),v=nlevels(TF),b=nlevels(BF)) 
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
            if (nlevels(BF)>1)
              reff=blockEstEffics(globTF,BF)[2] else
                reff=1
            if (isTRUE(all.equal(blocksEffBound,reff))) return(globTF)
          }  
        }
      } 
      if (r<searches) {
        for (iswap in 1:jumps) {
          counter=0
          repeat {
            counter=counter+1
            s1=sample(seq_len(length(TF)),1)
            z= seq_len(length(TF))[MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1]] 
            if (length(z)==0) next
            if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z)
            v11=VTT[c(TF[s[1]],TF[s[2]]),c(TF[s[1]],TF[s[2]])]
            v22=VBB[c(BF[s[1]],BF[s[2]]),c(BF[s[1]],BF[s[2]])]
            v12=VTB[c(TF[s[1]],TF[s[2]]),c(BF[s[1]],BF[s[2]])]
            Dswap=(1-2*sum(diag(v12))+sum(v12))**2-(2*sum(diag(v22))-sum(v22))*(2*sum(diag(v11))-sum(v11))
            if (Dswap>.1| counter>1000) break
          }
          if (counter>1000) return(globTF) # finish with no non-singular swaps
          relD=relD*Dswap 
          up=UpDate(VTT,VBB,VTB,TF[s[1]],TF[s[2]],BF[s[1]],BF[s[2]])
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
  # ********************************************************************************************************************
  # Random swaps
  # ********************************************************************************************************************    
  Swaps=function(TF,MF,BF,pivot,rank) {
    candidates=NULL
    while (is.null(candidates)) {
      if (rank<(length(TF)-1)) 
        s1=sample(pivot[(1+rank):length(TF)],1) else 
          s1=pivot[length(TF)]
        candidates = (1:length(TF))[ MF==MF[s1] & BF!=BF[s1] & TF!=TF[s1] ] 
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
    for ( i in 1:1000) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      swap=Swaps(TF,restrict,BF,pivot,rank)
      TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
      Q=qr(t(cbind(TM,BM)))
      if (Q$rank>rank) {
        TF[c(swap[1],swap[2])]=TF[c(swap[2],swap[1])]
        rank=Q$rank
        pivot=Q$pivot
      } else TM[c(swap[1],swap[2]),]=TM[c(swap[2],swap[1]),]
    }
    stop("Unable to find a non-singular solution for this design - try a simpler design")
  }
  # *****************************************************************************************************************
  # Optimizes the nested Blocks assuming a possible set of Main block constraints and assuming a randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # The model matrix for the full set of treatment or block effects is represented here by treatment contrasts which assume an overall 
  # mean equal to the first treatment or block level which is omitted to give a full rank model. 
  # The variance matrices are calculated by taking contrasts WITHIN each term to give the contrast of each term from the first term.
  # The include vector assigns the the rows and columns of the reduced inverse variance matrix to the appropriate positions for the full
  # inverse variance matrix with zero rows and columns for the constant row treatmnent and block effects 
  # The model.matrix assumes treatment type contrasts and the options setting for this function is treatment contrasts.
  # *****************************************************************************************************************
  blocksOpt=function(TF,MF,BF,regReps) {
   
    exclude=sapply(1:nlevels(MF),function(i) {(which(table(MF,BF)[i,]>0))[1]})# the level of the first nested block in each main block 
    include=(1:nlevels(BF))[!(1:nlevels(BF))%in%exclude] # nested blocks omitting first block in each main block)
    # center treatment effects and nested block effects within main blocks so that treatments are optimized within nested blocks 
    TM = model.matrix(~TF)[,-1,drop=FALSE] # Omits the constant treatment mean 
    TM = do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TM[MF==levels(MF)[i],],scale = FALSE)}))
    BM = model.matrix(~BF)[,include,drop=FALSE] # Omits the first nested blocks contrast in each main block
    BM = do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(BM[MF==levels(MF)[i],],scale = FALSE)}))
    if ((ncol(cbind(TM,BM))+1) > nrow(TM)) stop( paste("Too many parameters: plots =",nrow(TM)," parameters = ",ncol(cbind(TM,BM))))
    nonsing=NonSingular(TF,BF,TM,BM,MF)
    TF=nonsing$TF
    TM=nonsing$TM
    V = chol2inv(chol(crossprod(cbind(TM,BM))))
    VTT = matrix(0, nrow=nlevels(TF), ncol=nlevels(TF))
    VBB = matrix(0, nrow=nlevels(BF), ncol=nlevels(BF))
    VTB = matrix(0, nrow=nlevels(TF), ncol=nlevels(BF))
    VTT[2:nlevels(TF),2:nlevels(TF)] = V[1:(nlevels(TF)-1),1:(nlevels(TF)-1),drop=FALSE]
    VBB[include,include]             = V[nlevels(TF):ncol(V),nlevels(TF):ncol(V),drop=FALSE]
    VTB[2:nlevels(TF),include]       = V[1:(nlevels(TF)-1), nlevels(TF):ncol(V), drop=FALSE]
    TF=Optimise(TF,MF,BF,VTT,VBB,VTB,TM,BM,regReps)
    return(TF)
  }
  
  # *****************************************************************************************************************
  # nestedBlocks - valid for nested block designs ONLY
  # *****************************************************************************************************************
  TF=TF[,1]
  trtreps=as.matrix(table(TF))[,1]
  regReps=isTRUE(diff(range(trtreps))<tol)
  # the Base block is a single super-block and is the null 'restriction' factor for the first set of BF blocks 
  BF=cbind("Base" = factor(rep(1,nrow(BF))),BF)
  # orthogonality test for proportional treatment representation in each block
  orthBlocks = sapply(1:ncol(BF),function(i) {isTRUE(diff(range(scale(table(TF,BF[,i])/trtreps, center=TRUE,scale=FALSE)))<tol)})
  orthblocks=sum(orthBlocks)
  # fully randomized in bottom level of orthogonal blocks
  TF= unlist(lapply(1:nlevels(BF[,orthblocks]), function(j) {sample(TF[BF[,orthblocks]==levels(BF[,orthblocks])[j]])}))

  if(orthblocks<ncol(BF)) { 
    for (i in (orthblocks+1):ncol(BF)) {
      if (isTRUE(diff(range(table(BF[,i])))<tol) & regReps) {
        b=nlevels(BF[,i]) # number of incomplete blocks
        k=length(TF)/b # block size of incomplete blocks
        sk=sqrt(nlevels(TF))  # required block size for a square lattice 
        rk=(sqrt(1+4*nlevels(TF))-1)/2 #  required block size for a rectangular lattice
        ## necessary conditions for square or rectangular lattice designs 
        if (sk%%1==0 & k==sk & trtreps[1]<(sk+2) ) {
          TL=squarelattice(sk,trtreps[1]) 
        } else if (rk%%1==0 & k==rk & trtreps[1]<(rk+2)) { 
          TL=rectlattice(rk+1,trtreps[1]) 
        } else TL=NULL
      } else TL=NULL
      if (is.null(TL)) TF=blocksOpt(TF,BF[,(i-1)],BF[,i],regReps) else
        TF=TL 
    }
  }
  
  Design=data.frame(BF,plots=factor(1:length(TF)),treatments=TF)[,-1,drop=FALSE]
  Efficiencies=BlockEfficiencies(Design)
  blocksizes=table(BF[,ncol(BF)])
  BF=unique(BF[,-1 ,drop=FALSE])
  V = split( Design[,ncol(Design)], Design[,(ncol(Design)-2)],lex.order = TRUE) 
  V = lapply(V, function(x){ length(x) =max(blocksizes); x })
  Plan = data.frame(BF,rep("",length(V)),matrix(unlist(V),nrow=length(V),byrow=TRUE))
  colnames(Plan)=c(colnames(BF),"Blocks.Plots:", c(1:max(blocksizes)))
  row.names(Plan)=NULL
  row.names(Design)=NULL
  row.names(Efficiencies)=NULL
  list(Effic=Efficiencies,Design=Design,Plan=Plan)
}

