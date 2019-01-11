#' @title D-Efficiency factors
#' 
#' @description
#' Finds D-efficiency for general treatment and block designs.
#' 
#' @details
#' efficiency factors of regular block designs 
#' 
#' @param TF the treatments factor data frame
#' 
#' @param BF the block factors data frame
#' 
#' @param treatments_model a model formula for the required treatments design where the default 
#' formula assumes a fully crossed factorial treatment model.
#'  
#' @export
blockEfficiencies=function(TF,BF,treatments_model){
  TF=data.frame(TF)
  BF=data.frame(BF)
  nunits=nrow(TF)
  maxrank=nunits-ncol(TF)
  prodfactors =unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse="*")}) )
  Add_levs=integer(length=ncol(BF))
  Add_D_Effic=numeric(length=ncol(BF))
  Add_A_Effic=numeric(length=ncol(BF))
  Int_levs=integer(length=ncol(BF))
  Int_D_Effic=numeric(length=ncol(BF))
  Int_A_Effic=numeric(length=ncol(BF))
  
  if (ncol(TF)==1 & is.factor(TF[,1])) {
    for (i in 1:ncol(BF)) {
    Add_levs[i]=nlevels(BF[,i])
    if (nlevels(BF[,i])+nlevels(TF)<nunits) {
      effics=EstEffics(TF[,1],BF[,i])
      Add_D_Effic[i]=effics$Deffic
      Add_A_Effic[i]=effics$Aeffic
      } 
    inti=droplevels(interaction(BF[,c(1:i)]))
    Int_levs[i]=nlevels(inti)
    if (nlevels(inti)>nlevels(BF[,i]) &  nlevels(inti)+nlevels(TF)<nunits   ) {
        effics=EstEffics(TF[,1],inti)
        Int_D_Effic[i]=effics$Deffic
        Int_A_Effic[i]=effics$Aeffic
      } else if ( nlevels(inti)+nlevels(TF)<nunits   ) {
        Int_D_Effic[i]=Add_D_Effic[i]
        Int_A_Effic[i]=Add_A_Effic[i]
      } 
    }
  } else {
    TM=scale(model.matrix(as.formula(paste("~",treatments_model)),TF), center = TRUE, scale = FALSE)
    qTM=qr(TM)
    TM = TM[,qTM$pivot[1:qTM$rank],drop=FALSE]
    q = eigen(crossprod(TM))
    Z = crossprod(t(q$vectors),diag(sqrt(1/q$values)))
    for (i in 1:ncol(BF)) {
      addBM = scale(model.matrix(as.formula(paste("~", paste0(colnames(BF)[1:i],collapse="+"))),BF), 
                center = TRUE, scale = FALSE)
      Qadd=qr(addBM)
      if (Qadd$rank+ncol(TM)<nunits) {
        addBM = addBM[,Qadd$pivot[1:Qadd$rank],drop=FALSE]
        addBM = qr.Q(qr(addBM)) # orthogonal basis 
        TB=crossprod(TM,addBM)
        U=crossprod(Z,TB)
        l=eigen(diag(1,ncol(TM))-tcrossprod(U))$values
        Add_D_Effic[i]=exp(mean(log(l)))
        Add_A_Effic[i]=1/mean(1/l)
      } else {
        Add_D_Effic[i]=NA
        Add_A_Effic[i]=NA
      }
    twoFact=scale(model.matrix(as.formula(paste("~",paste0("(",paste0(colnames(BF)[1:i],collapse="+"),")^2"))),BF), 
                  center = TRUE, scale = FALSE)
    Qtwo=qr(twoFact)
    if (Qtwo$rank+ncol(TM)<nunits) {
        twoFact = twoFact[,Qtwo$pivot[1:Qtwo$rank],drop=FALSE]
        twoFact = qr.Q(qr(twoFact)) # orthogonal basis 
        TB=crossprod(TM,twoFact)
        U=crossprod(Z,TB)
        l=eigen(diag(1,ncol(TM))-tcrossprod(U))$values
        Int_D_Effic[i]=exp(mean(log(l)))
        Int_A_Effic[i]=1/mean(1/l)
      } else {
        Int_D_Effic[i]=NA
        Int_A_Effic[i]=NA
      }
    Add_levs[i]=Qadd$rank
    Int_levs[i]=Qtwo$rank
    }
  }
  Effics=data.frame(Blocks=colnames(BF),Add_levs,Add_D_Effic,Add_A_Effic,Interactions=prodfactors,Int_levs,Int_D_Effic,Int_A_Effic)
  return(Effics)
}