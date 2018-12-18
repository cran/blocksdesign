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
  tol = .Machine$double.eps ^ 0.5
  EstEffics=function(TF,BF) {
    k=nlevels(BF)
    if (k==1) return(c(1,1))
      M1= diag(nlevels(TF))-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))
      M2= diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))
    if (nlevels(TF)<=k)
      e=eigen( M1,symmetric=TRUE, only.values = TRUE)$values[1:(nlevels(TF)-1)] 
    else
      e=c(rep(1,(nlevels(TF)-k)),eigen(M2,symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])
    if (all(e>tol)) {
    D=exp(mean(log(e)))
    A=1/mean(1/e)
    } else {
      D=NA
      A=NA
    }
    return(list(Deffic= round(D,7), Aeffic=round(A,7)))
  }

  TF=data.frame(TF)
  BF=data.frame(BF)
  nunits=nrow(TF)
  maxrank=nunits-ncol(TF)
  addfactors =unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse="+")}) )
  prodfactors =unlist(lapply(1:ncol(BF),function(j){ paste0(colnames(BF)[1:j],collapse="*")}) )
  if (ncol(TF)==1 & is.factor(TF[,1])) {
    Effics=matrix(nrow=ncol(BF),ncol=8)
    for (i in 1:ncol(BF)) {
    Effics[i,1]=addfactors[i]
    Effics[i,2]=nlevels(BF[,i])
    if (nlevels(BF[,i])+nlevels(TF)<nunits) {
      effics=EstEffics(TF[,1],BF[,i])
      Effics[i,3]=effics$Deffic
      Effics[i,4]=effics$Aeffic
    } else {
      Effics[i,3]=NA
      Effics[i,4]=NA
    }
    inti=droplevels(interaction(BF[,c(1:i)]))
    Effics[i,5]=prodfactors[i]
    Effics[i,6]=nlevels(inti)
    if (nlevels(inti)>nlevels(BF[,i]) &  nlevels(inti)+nlevels(TF)<nunits   ) {
      effics=EstEffics(TF[,1],inti)

      Effics[i,7]=effics$Deffic
      Effics[i,8]=effics$Aeffic
    } else if ( nlevels(inti)+nlevels(TF)<nunits   ) {
      Effics[i,7]=Effics[i,3]
      Effics[i,8]=Effics[i,4]
    } else {
      Effics[i,7]=NA
      Effics[i,8]=NA
    }
    }
  } else {
    Effics=matrix(nrow=0,ncol=8)
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
      addDeffic=exp(mean(log(l)))
      addAeffic=1/mean(1/l)
    } else {
      addDeffic=NA
      addAeffic=NA
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
      twoFactDeffic=exp(mean(log(l)))
      twoFactAeffic=1/mean(1/l)
    } else {
      twoFactDeffic=NA
      twoFactAeffic=NA
    }
  
    Effics=rbind(Effics,
                 c(addfactors[i],Qadd$rank,round(addDeffic,7),round(addAeffic,7),
                   prodfactors[i],Qtwo$rank,round(twoFactDeffic,7),round(twoFactAeffic,7)))
    
  }
  }

 Effics=data.frame(Effics)
 colnames(Effics) = c("Additive model","levels","D-Efficiency","A-Efficiency",
                       "Interaction model","levels","D-Efficiency","A-Efficiency")
  return(Effics)
}