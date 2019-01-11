#' @title D-Efficiency factors
#' 
#' @description
#' Finds D and A-efficiency for an unstructured treatmnent set TF with blocks factor BF
#' 
#' @details
#' efficiency factors of regular block designs 
#' 
#' @param TF a treatments factor data frame
#' 
#' @param BF a block factors data frame
#'  
#' @export
  EstEffics=function(TF,BF) {
    tol = .Machine$double.eps ^ 0.5
    if (is.data.frame(TF))TF=TF[,1]
    if (is.data.frame(BF)) BF=BF[,1]
    k=nlevels(BF)
    if (k==1) return(c(1,1))
    if (nlevels(TF)<=k) {
      M1= diag(nlevels(TF))-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))
      e=eigen( M1,symmetric=TRUE, only.values = TRUE)$values[1:(nlevels(TF)-1)] 
    } else {
      M2= diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))
      e=c(rep(1,(nlevels(TF)-k)),eigen(M2,symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])
    }
    if (all(e>tol)) {
    D=exp(mean(log(e)))
    A=1/mean(1/e)
    } else {
      D=NA
      A=NA
    }
    return(list(Deffic= round(D,7), Aeffic=round(A,7)))
  }

 