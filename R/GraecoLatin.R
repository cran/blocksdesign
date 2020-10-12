#' @title Graeco-Latin squares
#' 
#' @description
#' Constructs mutually orthogonal Graeco-Latin squares for the following N:
#' 
#' i) any odd valued N
#' 
#' ii) any prime-power N = p**q where p and q can be  chosen from 
#'    \tabular{rrrrrrrr}{
#'    \bold{prime p} \tab                                       \bold{maximum q}\cr
#'    2 \tab                                                                13\cr
#'    3 \tab                                                                 8\cr
#'    5 \tab                                                                 6\cr                                                      
#'    7 \tab                                                                 5\cr
#'    11 \tab                                                                4\cr
#'    13 17 19 23 \tab                                                       3\cr
#'    29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 \tab                   2\cr
#'    Any prime >97 \tab                                                     1
#'    }
#' 
#' iii) any even valued N <= 30 except for 6 or 2
#' 
#' @details Plans are given for pairs of MOLS classified by rows and columns.
#' The output is a single data frame of size \eqn{p**q x (r+2)} for the required set of MOLS
#' with a column for the rows classification, a column for the columns classification and a 
#' column for each treatment set from the required set of MOLS.
#' 
#' Also see the function \code{MOLS} which will generate complete sets of MOLS for prime-power design sizes.
#'  
#' @seealso \code{\link{MOLS}}
#'   
#' @references
#' Street, A. P. & Street, D. J. (1987). Combinatorics of Experimental Design, Chapters 6 and 7.
#' Clarendon Press, Oxford.
#'  
#' @param N any suitable integer N
#' 
#' @return Data frame of factor levels for rows, columns and treatment sets
#'  
#' @examples
#' X=GraecoLatin(8) 
#' table(X[,3],X[,4])
#' X=GraecoLatin(9) 
#' table(X[,3],X[,4])
#' X=GraecoLatin(32)
#' table(X[,3],X[,4])
#' X=GraecoLatin(97)
#' table(X[,3],X[,4])
#'  
#' @export
GraecoLatin=function(N) {
  # odd valued N
  oddN=function(N) {
  v=0:(N-1)
  row=rep(1:N,each=N)
  col=rep(1:N,N)
  s1=sapply(1:N,function(i) {(v+i)%%N + 1})
  s2=sapply(1:N,function(i) {(v+2*i)%%N + 1})  
  mols=data.frame(row,col,as.numeric(s1),as.numeric(s2))
  colnames(mols)=c("Row","Col",paste("T", 1:2, sep = ""))
  return(mols)
  }
 # N = 10 or14 
  eulerN=function(g,d) {
  dim=length(g)
  G=matrix(g,ncol=dim,nrow=dim,byrow=TRUE)
  for (i in 2:dim) G[i,] = (G[(i-1),c(dim,1:(dim-1))]+1)%%dim
  colmarg=sapply(1:dim, function(i) G[i,(d+i-1)%%dim+1] )
  rowmarg=colmarg[((0:(dim-1))+(dim-d))%%dim+1 ]
  for (i in 1:dim) G[i,(d+i-1)%%dim+1]=dim
  G=cbind(G,colmarg)
  G=rbind(G,c(rowmarg,dim))
  G=G+1
  mols=data.frame(rep(1:(dim+1),each=(dim+1)),rep(1:(dim+1),(dim+1)),as.numeric(G),as.numeric(t(G)))
  mols
  }
  # N = 18, 22, 26 or 30
  EulerN=function(L1,l2,l3,l4) {
    N=length(l4)
    L2=matrix(l2,nrow=4,ncol=N)
    L3=matrix(l3,nrow=4,ncol=N)
    L4=matrix(l4,nrow=N,ncol=N,byrow=TRUE)
    for(i in 2 :N)
      L2[,i]= (L2[,(i-1)]+1)%%N
    for(i in 2 :N)
      L3[,i]= (L3[,(i-1)]+1)%%N
    for(i in 2:N){
      if (L4[(i-1),N]>(N-1)) L4[i,1]=L4[(i-1),N] else L4[i,1]=(L4[(i-1),N]+1)%%N
      for(j in 1:(N-1))
        if (L4[(i-1),j]>(N-1)) L4[i,(j+1)]=L4[(i-1),j] else L4[i,(j+1)]=(L4[(i-1),j]+1)%%N
    }
    S=rbind(cbind(L1,L2),cbind(t(L3),L4))
    mols=data.frame(rows=rep(1:(N+4),each=(N+4)),cols=rep(1:(N+4),(N+4)),S1=as.numeric(S),S2=as.numeric(t(S)))
    return(mols)
  }
  ## Composite N divisible by 2**q for q>1 
  compN=function(N) {
    q1=1
    while(N%%2**q1==0) q1=q1+1 
    q1=q1-1 
    if (q1<2) stop("N is not divisible by 2**q for q>1 ")
    N1=2**q1
    N2=N/N1
    S1=MOLS(2,q1,2) 
    if (isPrime(N2)) S2=MOLS(N2,1,2) else if (N2%%2!=0) S2=oddN(N2) else
      stop("cannot find a suitable pair of MOLS for n2 and q2")
    
    S1g=matrix(S1[,3],nrow=N1,ncol=N1)
    S1l=matrix(S1[,4],nrow=N1,ncol=N1)
    S2g=matrix(S2[,3],nrow=N2,ncol=N2)
    S2l=matrix(S2[,4],nrow=N2,ncol=N2)
    Sg= as.numeric(kronecker(matrix(1,nrow=N1,ncol=N1),S2g) + kronecker((S1g-1)*N2,matrix(1,nrow=N2,ncol=N2)))
    Sl= as.numeric(kronecker(matrix(1,nrow=N1,ncol=N1),S2l) + kronecker((S1l-1)*N2,matrix(1,nrow=N2,ncol=N2)))
    mols=data.frame( rep(1:N,each=N),rep(1:N,N),Sg,Sl)
    colnames(mols)=c("Row","Col",paste("T", 1:2, sep = ""))
    return(mols)
  }
  ##  Any prime or prime-power N
  X=isPrimePower(N)
  if (!is.null(X)) { 
    mols=MOLS(X$base,X$power,2) 
    return(mols)
  }
  ## Anyy odd valued non-prime power N
  if (N%%2!=0) {
    mols=oddN(N) 
    return(mols)
  }
  ##  Euler' sizes for N = 4n + 2 and n<=7
  if (N==10) {
    mols=eulerN(c(0,2,4,7,1,8,5,3,6),6)
    return(mols)
    }
  if (N==14) {
    mols=eulerN(c(0,9,5,8,10,2,4,6,12,3,11,7,1),12)
    return(mols)
    }
  if (N==18) {
    L1=matrix(c(14,17,15,16,16,15,17,14,17,14,16,15,15,16,14,17),nrow=4,ncol=4)
    l2=c(13,3,12,2)
    l3=c(8,5,4,6)
    l4=c(0,7,13,12,11,10,2,1,9,17,16,15,14,3)
    mols=EulerN(L1,l2,l3,l4)
    return(mols)
  }
  if (N==22) {
    L1=matrix(c(18,21,19,20,20,19,21,18,21,18,20,19,19,20,18,21),nrow=4,ncol=4)
    l2=c(7,5,4,16)
    l3=c(2,3,10,12)
    l4=c(0,9,17,16,15,14,8,6,11,1,4,7,13,21,20,19,18,5)
    mols=EulerN(L1,l2,l3,l4)
    return(mols)
  }
  if (N==26) {
    L1=matrix(c(22,25,23,24,24,23,25,22,25,22,24,23,23,24,22,25),nrow=4,ncol=4)
    l2=c(10,20,7,3)
    l3=c(3,6,4,2)
    l4=c(0,9,21,20,19,18,15,12,10,13,16,1,11,14,8,7,5,25,24,23,22,17)
    mols=EulerN(L1,l2,l3,l4)
    return(mols)
  }
  if (N==30) {
    L1=matrix(c(26,29,27,28,28,27,29,26,29,26,28,27,27,28,26,29),nrow=4,ncol=4)
    l2=c(13,8,15,18)
    l3=c(2,3,6,16)
    l4=c(0,7,25,24,23,22,4,11,15,18,21,12,14,1,8,20,19,13,17,5,10,29,28,27,26,9)
    mols=EulerN(L1,l2,l3,l4)
    return(mols)
  }
  ## Composite N divisible by 2**q for q>2 
  if (N==12) {mols=compN(12)
  return(mols)}
  if (N==20) {mols=compN(20)
  return(mols)}
  if (N==24) {mols=compN(24)
  return(mols)}
  if (N==28) {mols=compN(28)
    return(mols)}
}
