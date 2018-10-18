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
  Effics = matrix(nrow=0, ncol=7)
  TF=data.frame(TF)
  BF=data.frame(BF)
  nunits=nrow(TF)
  maxrank=nunits-ncol(TF)
  factors =lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse=",")}) 
  if (is.null(treatments_model)) treatments_model=paste("~",paste(colnames(TF),collapse="*")) 
  TM=scale(model.matrix(as.formula(treatments_model),TF), center = TRUE, scale = FALSE)[,-1,drop=FALSE]
  Dbase=determinant(crossprod(TM),logarithm = TRUE)$modulus/ncol(TM)
  eig=eigen(crossprod(TM), only.values = TRUE)[[1]]
  DTeffic=prod(eig)**(1/length(eig))
  ATeffic=sum(1/eig)
  for (i in 1:ncol(BF)) {
    addBM=scale(model.matrix(as.formula(paste("~", paste0(colnames(BF)[1:i],collapse="+"))),BF), center = TRUE, scale = FALSE)
    Qadd=qr(addBM)
    if (Qadd$rank+ncol(TM)<nunits) {
      addBM = addBM[,Qadd$pivot[1:Qadd$rank],drop=FALSE]
      addBM = qr.Q(qr(addBM)) # orthogonal basis 
      H=diag(rep(1,nrow(addBM)))-tcrossprod(addBM) 
      Info=crossprod(t(crossprod(TM,H)),TM)
      eig=eigen(Info, only.values = TRUE)[[1]]
      addDeffic=prod(eig)**(1/length(eig))/DTeffic
      addAeffic=ATeffic/sum(1/eig)
    } else {
      addDeffic=NA
      addAeffic=NA
    }
    
    twoFact=scale(model.matrix(as.formula(paste("~",paste0("(",paste0(colnames(BF)[1:i],collapse="+"),")^2"))),BF), center = TRUE, scale = FALSE)
    Qtwo=qr(twoFact)
    if (Qtwo$rank+ncol(TM)<nunits) {
      twoFact = twoFact[,Qtwo$pivot[1:Qtwo$rank],drop=FALSE]
      twoFact = qr.Q(qr(twoFact)) # orthogonal basis 
      H=diag(rep(1,nrow(twoFact)))-tcrossprod(twoFact) 
      Info=crossprod(t(crossprod(TM,H)),TM)
      eig=eigen(Info, only.values = TRUE)[[1]]
      twoFactDeffic=prod(eig)**(1/length(eig))/DTeffic
      twoFactAeffic=ATeffic/sum(1/eig)
    } else {
      twoFactDeffic=NA
      twoFactAeffic=NA
    }
    Effics=rbind(Effics,
                 c(factors[i],Qadd$rank,round(addDeffic,4),round(addAeffic,4),Qtwo$rank, round(twoFactDeffic,4),round(twoFactAeffic,4)))
  }
  Effics=data.frame(Effics)
  colnames(Effics)= c("Factors", "Main effects df","D-Effic.", "A-Effic.","2-factor effects df","D-Effic.","A-Effic.")
  return(Effics)
}
