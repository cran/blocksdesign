#' @title Efficiency factors
#' 
#' @description
#' Finds efficiency bounds for block designs.
#' 
#' @details
#' efficiency factors of regular block designs 
#' 
#' @param TF the treatments factor data frame
#' 
#' @param BF the block factors data frame
#' 
#' @param treatments_model the treatments model formula for the treatments data frame. Default is a full factorial model. 
#' 
#' @export
  blockEfficiencies=function(TF,BF,treatments_model=NULL) {
  Effics = matrix(nrow=0, ncol=6)
  TF=data.frame(TF)
  BF=data.frame(BF)
  nunits=nrow(TF)
  maxrank=nunits-ncol(TF)
  
  addblocks =lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse="+")})
  multblocks=lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse=".")})
  multBF=BF
  if (ncol(multBF)>1)
    for (i in 2:ncol(multBF)) multBF[,i]= droplevels(interaction( multBF[,i-1],BF[,i] ))
  colnames(multBF)=multblocks
  multBF=cbind(mean=factor(rep(1,nrow(multBF))),multBF)
  
  if (is.null(treatments_model)) treatments_model=paste("~",paste(colnames(TF),collapse="*")) 
  TM=scale(model.matrix(as.formula(treatments_model),TF), center = TRUE, scale = FALSE)[,-1,drop=FALSE]
  Dbase=determinant(crossprod(TM),logarithm = TRUE)$modulus/ncol(TM)

  for (i in 1:ncol(BF)) {
    aBM=scale(model.matrix(as.formula(paste("~",addblocks[[i]])),BF), center = TRUE, scale = FALSE)[,-1,drop=FALSE]
    Qa=qr(aBM)
    dfa=Qa$rank
    if (Qa$rank+ncol(TM)<nunits) {
      aBM = aBM[,Qa$pivot[1:Qa$rank],drop=FALSE]
      aBM = qr.Q(qr(aBM)) # orthogonal basis 
      H=diag(rep(1,nrow(aBM)))-tcrossprod(aBM) 
      addDeffic=determinant(crossprod(t(crossprod(TM,H)),TM),logarithm = TRUE)$modulus/ncol(TM)
    } else addDeffic=NA
    dfm=nlevels(multBF[,i+1])-1
    if (nlevels(multBF[,i+1]) < nunits) {
      mBM=scale(model.matrix(as.formula(paste("~",multblocks[[i]])),multBF), center = TRUE, scale = FALSE)[,-1,drop=FALSE]
      Q=qr(mBM)
      dfm=Q$rank
     if (Q$rank+ncol(TM)<nunits) {
      mBM=mBM[,Q$pivot[1:Q$rank],drop=FALSE]
      mBM = qr.Q(qr(mBM)) # orthogonal basis
      H=diag(rep(1,nrow(mBM)))-tcrossprod(mBM) 
      multDeffic=determinant(crossprod(t(crossprod(TM,H)),TM),logarithm = TRUE)$modulus/ncol(TM)
    } else multDeffic=NA
    } else multDeffic=NA
    Effics=rbind(Effics,
                c(addblocks[[i]], dfa,  round(exp(addDeffic-Dbase),4),
                   multblocks[[i]], dfm, round(exp(multDeffic-Dbase),4)))
  }
  Effics=data.frame(Effics)
  colnames(Effics)= c("Additive_model","df","D-Effic","Multiplicative model","df","D-Effic")

  return(Effics)
}
