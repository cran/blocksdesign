#' @title Block designs 
#' 
#' @description
#' 
#' Constructs randomized nested and crossed block designs for factorial or unstructured treatment sets where blocks can have any feasible depth of nesting.
#' 
#' @details
#' 
#' \code{blocksdesign} constructs nested block designs with arbitrary depth of nesting. The top or zero-level stratum is assumed to be a single super-block 
#' and the blocks of each succesive strata are nested hierarchically within the blocks of each preceding stratum. Strata are optimized sequentially from the top
#' down with the blocks of each stratum optimized conditionally within the blocks of the immediately preceding stratum.
#' 
#' The treatment design can be either a single unstructured treatment set with an arbitrary number of treatments and an arbitrary number of replications per
#'  treatment or a factorial treatment set with an arbitrary combination of quantitative or qualitative level treatment factors and an arbitrary factorial model.
#' 
#' For an unstructured treatment set, the \code{treatments} parameter is a set of cardinal numbers that partition the total number of treatments into
#' sub-sets of equally replicated treatments. In this mode, the \code{replicates} parameter is a set of replication numbers with one replication 
#' number for each equally replicated treatment set. 
#' 
#' For a factorial treatment set, the \code{treatments} parameter is a data frame containing a suitable factorial treatment design where the rows 
#' are individual treatment combinations and the columns are individual treatment factors, either qualitative or quantitative.
#' In this mode, the \code{replicates} parameter is a single number gving the number of complete repeats of the data frame in the full treatment design. 
#' 
#' The model for a factorial treatment design is defined by a model formula which is based on the notation used in the \code{\link[stats]{model.matrix}} package.
#' The default model is a complete set of factorial effects but the \code{models} formula can be used to define any feasible reduced set of factorial effects. 
#' The \code{models} formula has no effect on single unstructured treatment sets.   
#' 
#' The \code{rows} parameter, if any, defines the nested row blocks in each nested stratum taken in order from the highest to the lowest.
#' The first number, if any, is the number of rows in the blocks of the first-level stratum, the second, if any, is the number of rows in the blocks of
#'  the second-level stratum and so on for all the required strata in the design. 
#'  
#' The \code{columns} parameter, if any, defines the nested column blocks in each nested stratum taken in order from the highest to the lowest. 
#' The first number, if any, is the number of columns in the blocks of the first-level stratum, the second, if any, is the number of columns in the blocks of
#'  the second-level stratum and so on for all the required strata in the design. 
#'  
#' The default block design for any particular stratum is a crossed row-and-column design but if a single column block is selected for that stratum, 
#' the design will reduce to a nested blocks design with a single set of nested blocks in that stratum. 
#'   
#' The \code{rows} and \code{columns} parameters, if defined, must be of equal length and if a simple nested blocks design is required for 
#' any particular choice of stratum, the number of columns for that stratum should be set to unity. Any required combination of simple or 
#' crossed blocks can be obtained by appropriate choice of the levels of the \code{rows} and \code{columns} parameters.
#' 
#' If both the \code{rows} parameter and the \code{columns} parameter are null, the default block design will be a single set of orthogonal
#' main blocks equal in number to the highest common factor of the replication numbers. If the \code{rows} parameter is defined but the \code{columns} parameter
#' is null, the design will comprise a set of nested blocks in each stratum, as defined by the \code{rows} parameter.
#' 
#' Block sizes are always as nearly equal as possible and will never differ by more than a single plot in any particular block classification. Row blocks and 
#' column blocks must always contain at least two plots per block and this restriction will constrain the permitted numbers of rows and columns in the various 
#' strata of a design.
#' 
#' Unreplicated treatments are allowed and any simple nested block design can be augmented by any number of single unreplicated treatments to give augmented 
#' blocks that never differ in size by more than a single plot. General crossed block designs are more complex and currently the algorithm will only accommodate 
#' single unreplicated treatments in a crossed block design if the block sizes of the replicated part of the design are all equal in each stratum of the design.
#' 
#'  For any particular stratum, the algorithm will first optimize the row blocks in that stratum conditional on the blocks of any immediately preceding stratum, 
#'  and will then optimise the columns blocks, if any, conditional on the rows blocks. 
#'  
#'  Special designs:
#'  
#' Trojan designs are row-and-column designs for p replicates of v*p treatments arranged in p-rows and p-columns where v<p and where every row x column intersection 
#' block contain v plots. Trojan designs have orthogonal rows and columns and optimal rows x columns intersection  blocks and exist whenever p is prime or prime-power. 
#' \code{blocksdesign} constructs these designs algebraically from mutually orthogonal Latin squares (MOLS).  
#'  
#' Square lattice designs are efficient resolvable incomplete block designs for r replicates of p*p treatments arranged in blocks of size p where r < p+2
#' for prime or prime power p or r < 4 for general p. \code{blocksdesign} constructs these designs algebraically from Latin squares or MOLS.
#' 
#'  Lattice designs and Trojan designs based on prime-power MOLS require the \code{\link[crossdes]{MOLS}} package.
#'     
#'  All other designs are constructed algorithmically.
#'  
#' Warnings:
#'  
#'  Row-and-column designs may contain useful treatment information in the individual row-by-column intersection blocks but \code{blocksdesign} does not currently 
#'  optimize the efficiency of these blocks in row-and-column designs except for the special case of Trojan designs.
#'
#'  Row-and-column design with 2 complete treatment replicates, 2 complete rows and 2 complete columns will always confound one treatment contrast in the 
#'  rows-by-columns interaction and for these designs, it is impossible to nest a non-singular block design in the rows-by-columns intersections. 
#'  Instead, we recommend a randomized nested blocks design with four incomplete main blocks. 
#'  
#'  Outputs:
#'  
#'  The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
#'  \item  A table showing the replication number of each treatment in the design. \cr
#'  \item  A table showing the block levels and the achieved D-efficiency and A-efficiency (unstructured treatment 
#'  designs only) factors for each stratum together with A-efficiency upper bounds, where available . \cr
#'  \item  A plan showing the allocation of treatments to blocks or to rows and to columns in the bottom stratum of the design (unstructured treatment 
#'  designs only).\cr
#' } 
#' 
#' @param treatments  either a data frame with columns for the individual treatment factors or a set of cardinal numbers giving a partition 
#' of the total required number of treatments into sets of equally replicated treatments.
#' 
#' @param replicates  either a single replication number if the \code{treatments} parameter is a data frame or a set of replication numbers, one
#' per replication set, if the \code{treatments} parameter is a partition into equally replicated treatment sets
#' 
#' @param rows  a set of nested row block levels for the row blocks in each succesive stratum of the blocks design taken in order from the highest to the lowest. 
#' The default is the hcf of the replication numbers.
#' 
#' @param columns a set of nested column block levels for the column blocks in each succesive stratum of the blocks design taken in order from the highest to the lowest. 
#' The \code{rows} and the \code{columns} parameters, if both present, must be of equal length. The null default is a single column block for each nested stratum.  
#' 
#' @param model  a model equation for the treatment factors in the design where the equation is defined by the model.matrix notation
#' in the {\link[stats]{model.matrix}} package. If undefined, the model is a full factorial model. 
#' 
#' @param seed  an integer initializing the random number generator. The default is a random seed.
#' 
#' @param searches  the maximum number of local optima searched for a design optimization. The default is 1 plus the floor of 10000 divided by the number of plots.
#' 
#' @param jumps  the number of pairwise random treatment swaps used to escape a local maxima. The default is a single swap.
#' 
#' @return  
#' \item{Treatments}{The treatment factors defined by the \code{treatments} inputs in standard factorial order.}
#' \item{model.matrix}{The model.matrix used to define the \code{treatments} design.}
#' \item{Design}{Data frame giving the optimized block and treatment factors in plot order.}
#' \item{Plan}{Data frame for single factor designs showing a plan view of the treatment design in the bottom stratum of the design. A NULL plan is returned for multi-factor designs.}
#' \item{Efficiencies}{The achieved D-efficiencies and A-efficiencies (unstructured treatment designs only) for each stratum of the design together with an A-efficiency upper-bound, where available}
#' \item{seed}{Numerical seed for random number generator}
#' \item{searches}{Maximum number of searches in each stratum}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima}
#' 
#' @references
#' 
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. https://CRAN.R-project.org/package=crossdes
#' 
#' Edmondson R. N. (1998). Trojan square and incomplete Trojan square designs for crop research. Journal of Agricultural Science, Cambridge, 131, pp.135-142
#' 
#' Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
#' 
#' @examples
#' 
#' # First-order model for five 2-level factors with 2 main and 2 x 2 nested row-and-column blocks 
#' treatments =data.frame( F1=gl(2,16), F2=gl(2,8,32),  F3=gl(2,4,32), F4=gl(2,2,32) , F5=gl(2,1,32))
#' blocks(treatments=treatments,model="~ F1+F2+F3+F4+F5",rows=c(2,2),columns=c(1,2),searches=5)
#' 
#' # Full factorial model for two 2-level factors with three replicates in 6 randomized blocks 
#' treatments =data.frame( f1=gl(2,6,12), f2=gl(2,3,12))
#' blocks(treatments=treatments,rows=6,searches=5) # incomplete blocks with .6667 efficiency
#' 
#' # Quadratic regression for one 6-level numeric factor in 2 randomized blocks
#' blocks(treatments=data.frame(X=c(1:6)),model=" ~ (X + I(X^2))",rows=2,searches=5) 
#' 
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' TF=data.frame( F1=gl(2,16), F2=gl(2,8,32),  F3=gl(2,4,32), F4=gl(2,2,32) , F5=gl(2,1,32) )
#' blocks(treatments=TF,model=" ~ (F1+F2+F3+F4+F5)*(F1+F2+F3+F4+F5)",rows=4,searches=5)
#' 
#' # Second-order model for four qualitative 3-level factors in 9 randomized blocks
#' TF=data.frame( F1=gl(3,27), F2=gl(3,9,81),  F3=gl(3,3,81), F4=gl(3,1,81)  )
#' \dontrun{blocks(treatments=TF,model=" ~ (F1+F2+F3+F4)*(F1+F2+F3+F4)",rows=9)}
#' 
#' # Second-order model for two qualitative and two quantitative factors in 4 randomized blocks 
#' TF=data.frame(F1=gl(2,36), F2=gl(3,12,72), V1=rep(rep(1:3,each=4),6), V2=rep(1:4,18))
#' modform=" ~ F1*F2 + V1*V2 + I(V1^2) + I(V2^2) + F1:V1 + F1:V2 + F2:V1 + F2:V2"
#' blocks(treatments=TF,model=modform,rows=4,searches=10) 
#' 
#' # Unequal replication with hcf = 1 gives  default design with 1 fully randomised main block 
#' # 3 treatments x 2 replicates + 2 treatments x 4 replicates + 4 treatments x 3 replicates  
#' blocks(treatments=c(3,2,4),replicates=c(2,4,3))
#' 
#' # 4 treatments x 4 replicates with 2 main rows each containing two complete replicates  
#' blocks(treatments=4,replicates=4,rows=2)
#' 
#' # 50 treatments x 4 replicates with 4 main blocks and 5 nested sub-blocks in each main block 
#' blocks(treatments=50,replicates=4,rows=c(4,5))
#' 
#' # as above but with 20 single replicate treatments giving one extra treatment per sub-block
#' blocks(treatments=c(50,20),replicates=c(4,1),rows=c(4,5))
#' 
#' # 6 replicates of 6 treatments in 4 blocks of size 9 (non-binary block design)
#' blocks(treatments=6,replicates=6,rows=4)
#' 
#' # 4 replicates of 13 treatments arranged in a 13 x 4 Youden rectangle 
#' blocks(treatments=13,replicates=4,rows=13,columns=4)
#' 
#' # 64 treatments x 2 replicates with nested 8 x 8 row-and-column designs in two main blocks 
#' blocks(treatments=64,replicates=2,rows=c(2,8),columns=c(1,8)) 
#' 
#' # 64 treatments x 2 replicates with two main blocks and a 4 x 4 row-and-column in each main block
#' blocks(treatments=64,replicates=2,rows=c(2,4),columns=c(1,4),searches=12) 
#' 
# # 64 treatments x 4 replicates with succesively nested 2 x 2 row-and-column blocks
#' \dontrun{ blocks(treatments=64,replicates=4,rows=c(2,2,2,2),colums=c(2,2,2,2) }
#'      
#'                  
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames
#' 
blocks = function(treatments,replicates=1,rows=NULL,columns=NULL,model=NULL,searches=NULL,seed=sample(10000,1),jumps=1) { 
  options(contrasts=c('contr.SAS','contr.poly'))
  tol=.Machine$double.eps^0.5
  
  # ******************************************************************************************************************************************************** 
  # Finds the highest common factor (hcf) of a set of numbers omitting any zero values (Euclidean algorithm)
  # ********************************************************************************************************************************************************
  HCF=function(replicates)  {
    replicates=sort(replicates[replicates>0])
    for (i in  seq_len(length(replicates)))
      while (!isTRUE(all.equal(replicates[i]%%replicates[1],0))) replicates[c(1,i)] = c(replicates[i]%%replicates[1], replicates[1])
      replicates[1]
  }    
  # ******************************************************************************************************************************************************** 
  # Tests a given number for primality and returns TRUE or FALSE
  # ********************************************************************************************************************************************************
  isPrime=function(v) {
    if (v < 4) return(TRUE) 
    if ( isTRUE(all.equal(v %% 2,0)) ||  isTRUE(all.equal(v %% 3,0)) ) return(FALSE) 
    if (v<25) return(TRUE)
    for(i in  6*seq_len(length(floor((sqrt(v)+1)/6)))        )
      if ( isTRUE(all.equal(v %% (i-1) , 0)) ||   isTRUE(all.equal(v %% (i+1) , 0)) ) return(FALSE) 
    return(TRUE)
  } 
  # ********************************************************************************************************************************************************
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
  # ********************************************************************************************************************************************************
  factUpDate=function(MTT,MBB,MTB,t,b) {
    Mttt=crossprod(MTT,t)
    Mbbb=crossprod(MBB,b)
    Mtbt=crossprod(MTB,t)
    Mtbb=crossprod(t(MTB),b)
    tMt=as.numeric(crossprod(t,Mttt))
    bMb=as.numeric(crossprod(b,Mbbb))
    tMb=as.numeric(crossprod(b,Mtbt))
    f1=(Mttt+Mtbb)/sqrt(2)
    f2=(Mbbb+Mtbt)/sqrt(2)
    g1=(Mtbb-Mttt)/sqrt(2)
    g2=(Mbbb-Mtbt)/sqrt(2)
    a=(tMt+bMb+2*tMb)/2
    b=(tMt+bMb-2*tMb)/2
    c=(bMb-tMt)/2
    d=(1+a)*(1-b)+c*c
    MTT=MTT- (tcrossprod(f1)*(1-b) - tcrossprod(g1)*(1+a) + (tcrossprod(g1,f1)+tcrossprod(f1,g1))*c)/d
    MBB=MBB- (tcrossprod(f2)*(1-b) - tcrossprod(g2)*(1+a) + (tcrossprod(g2,f2)+tcrossprod(f2,g2))*c)/d
    MTB=MTB- (tcrossprod(f1,f2)*(1-b) - tcrossprod(g1,g2)*(1+a) + (tcrossprod(g1,f2)+tcrossprod(f1,g2))*c)/d
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }
  # ******************************************************************************************************************************************************** 
  # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb   
  # ********************************************************************************************************************************************************
  UpDate=function(MTT,MBB,MTB,ti,tj,bi,bj) { 
    mtt=MTT[ti,ti]+MTT[tj,tj]-2*MTT[tj,ti]
    mbb=MBB[bi,bi]+MBB[bj,bj]-2*MBB[bi,bj]
    mtb=1-MTB[ti,bi]+MTB[tj,bi]+MTB[ti,bj]-MTB[tj,bj]  
    TBbij=MTB[,bi]-MTB[,bj]
    TBtij=MTB[ti,]-MTB[tj,]
    TTtij=MTT[,ti]-MTT[,tj]
    BBbij=MBB[bi,]-MBB[bj,]
    Z1 = (TBbij-TTtij)/sqrt(2*mtb+mtt+mbb)   
    Z2 = (BBbij-TBtij)/sqrt(2*mtb+mtt+mbb)
    W1 = ( sqrt(2*mtb+mtt+mbb)*(TTtij+TBbij) - (mbb-mtt)*Z1) /(2*sqrt(mtb**2-mtt*mbb))
    W2 = ( sqrt(2*mtb+mtt+mbb)*(TBtij+BBbij) - (mbb-mtt)*Z2) /(2*sqrt(mtb**2-mtt*mbb))
    MTT = MTT - tcrossprod(Z1) + tcrossprod(W1)
    MBB = MBB - tcrossprod(Z2) + tcrossprod(W2)
    MTB = MTB - tcrossprod(Z1,Z2) + tcrossprod(W1,W2) 
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }  
  # ******************************************************************************************************************************************************** 
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
  # ********************************************************************************************************************************************************
  DMax=function(MTT,MBB,MTB,TF,Restrict,BF,TM){  
    locrelD=1
    mainSizes=tabulate(Restrict)
    nSamp=pmin( rep(8,nlevels(Restrict)), mainSizes)
    repeat {
      improved=FALSE
      for (k in 1: nlevels(Restrict)) {
        s=sort(sample( seq_len(nrow(TF)) [Restrict==levels(Restrict)[k]], nSamp[k])) 
        TB=MTB[TF[s,1],BF[s],drop=FALSE]-tcrossprod(MTB[cbind(TF[s,1],BF[s])],rep(1,nSamp[k]))
        dMat=(TB+t(TB)+1)**2-
          (2*MTT[TF[s,1],TF[s,1],drop=FALSE]- tcrossprod(MTT[cbind(TF[s,1],TF[s,1])]+rep(1,nSamp[k])) + tcrossprod(MTT[cbind(TF[s,1],TF[s,1])]) + 1)*
          (2*MBB[BF[s],BF[s],drop=FALSE]- tcrossprod(MBB[cbind(BF[s],BF[s])]+rep(1,nSamp[k])) + tcrossprod(MBB[cbind(BF[s],BF[s])]) + 1)
        if ( max(dMat, na.rm = TRUE)<1+tol ) next
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nSamp[k]
        j=1+(sampn-1)%/%nSamp[k]
        improved=TRUE
        locrelD=locrelD*dMat[i,j]
        up=UpDate(MTT,MBB,MTB,TF[s[i],1],TF[s[j],1],BF[s[i]],BF[s[j]])
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TF[c(s[i],s[j]),]=TF[c(s[j],s[i]),]
        TM[c(s[i],s[j]),]=TM[c(s[j],s[i]),]
      } 
      if (improved) next
      if (sum(nSamp) == nrow(TF))  break
      nSamp=pmin(mainSizes,2*nSamp)
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,TF=TF,locrelD=locrelD,TM=TM)
  } 
  # ******************************************************************************************************************************************************** 
  # Maximises the design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
  # ********************************************************************************************************************************************************
  factDMax=function(MTT,MBB,MTB,TF,Main,BF,TM,BM) {  
    locrelD=1
    mainSizes=tabulate(Main)
    nSamp=pmin( rep(8,nlevels(Main)), mainSizes)
    repeat {
      improved=FALSE
      for (k in 1: nlevels(Main)) {
        s=sort(sample( seq_len(nrow(TF)) [Main==levels(Main)[k]], nSamp[k])) 
        TMB=crossprod(t(crossprod(t(TM[s,]),MTB)),t(BM[s,]))
        TMT=crossprod(t(crossprod(t(TM[s,]),MTT)),t(TM[s,]))
        BMB=crossprod(t(crossprod(t(BM[s,]),MBB)),t(BM[s,]))
        TMB=sweep(TMB,1,diag(TMB))
        TMT=sweep(TMT,1,diag(TMT))
        BMB=sweep(BMB,1,diag(BMB))
        dMat=(1+TMB+t(TMB))**2 - (TMT + t(TMT))*(BMB + t(BMB))
        if ( max(dMat, na.rm = TRUE)<1+tol ) next
        sampn=which.max(dMat)
        i=1+(sampn-1)%%nrow(dMat)
        j=1+(sampn-1)%/%nrow(dMat)
        improved=TRUE
        locrelD=locrelD*dMat[i,j]
        t=TM[s[i],]-TM[s[j],]
        b=BM[s[j],]-BM[s[i],]
        up=factUpDate(MTT,MBB,MTB,t,b)
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TF[c(s[i],s[j]),]=TF[c(s[j],s[i]),]
        TM[c(s[i],s[j]),]=TM[c(s[j],s[i]),]
      } 
      if (improved) next
      if (sum(nSamp) == nrow(TF))  break
      nSamp=pmin(mainSizes,2*nSamp)
    }
    list(MTT=MTT,MBB=MBB,MTB=MTB,TF=TF,locrelD=locrelD,TM=TM)
  }  
  # ********************************************************************************************************************************************************
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,Main,BF,Restrict,MTT,MBB,MTB,TM,BM) {
    globrelD=0
    relD=1
    globTF=TF
    breps=tabulate(BF)  
    if ( regReps && max(breps)==min(breps)  && simpleTF)
      bound=upper_bounds( nrow(TF), nlevels(TF[,1]), nlevels(BF) ) else bound=NA
    
    if (!is.na(bound) && isTRUE(all.equal(bound,EstEffics(globTF[,1],BF)[2]))) return(globTF)
    if ( !simpleTF && isTRUE(all.equal(1,FactEstEffics(globTF,Main,BF))) ) return(globTF)
    for (r in 1:searches) {
      if (ncol(treatments)>1 | !is.factor(treatments[,1])) 
        dmax =factDMax(MTT,MBB,MTB,TF,Restrict,BF,TM,BM) 
      else
        dmax =DMax(MTT,MBB,MTB,TF,Restrict,BF,TM)
      if (dmax$locrelD>(1+tol)) {
        relD=relD*dmax$locrelD
        TF=dmax$TF
        TM=dmax$TM
        MTT=dmax$MTT
        MBB=dmax$MBB
        MTB=dmax$MTB 
        if (relD>globrelD) {
          globTF=TF
          globrelD=relD
          if (!is.na(bound) && isTRUE(all.equal(bound,EstEffics(globTF[,1],BF)[2]))) return(globTF)
          if (!simpleTF && isTRUE(all.equal(1,FactEstEffics(globTF,Main,BF))) ) return(globTF)
        }
      }
      if (r==searches) return(globTF)
      for (iswap in 1:jumps) {
        counter=0
        repeat {  
          counter=counter+1
          s1=sample(seq_len(nunits),1)
          available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
          z= seq_len(nunits)[Main==Main[s1] & Restrict==Restrict[s1] & BF!=BF[s1] & available]   
          if (length(z)==0) next
          if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z)
          if (ncol(treatments)>1 | !is.factor(treatments[,1])) {
            TMB=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],MTB)),BM[s[2],]-BM[s[1],] )
            TMT=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],MTT)),TM[s[2],]-TM[s[1],])
            BMB=crossprod(t(crossprod(BM[s[1],]-BM[s[2],],MBB)),BM[s[2],]-BM[s[1],])
            Dswap=(1+TMB)**2-TMT*BMB
          } else 
            Dswap = (1+MTB[TF[s[1],],BF[s[2]]]+MTB[TF[s[2],],BF[s[1]]]-MTB[TF[s[1],],BF[s[1]]]-MTB[TF[s[2],],BF[s[2]]])**2-
            (2*MTT[TF[s[1],],TF[s[2],]]-MTT[TF[s[1],],TF[s[1],]]-MTT[TF[s[2],],TF[s[2],]])*(2*MBB[BF[s[1]],BF[s[2]]]-MBB[BF[s[1]],BF[s[1]]]-MBB[BF[s[2]],BF[s[2]]]) 
          if (Dswap>.1| counter>1000) break
        }
        if (counter>1000) return(globTF) # finish with no non-singular swaps
        relD=relD*Dswap 
        if (ncol(treatments)>1 | !is.factor(treatments[,1])) 
          up=factUpDate(MTT,MBB,MTB, TM[s[1],]-TM[s[2],], BM[s[2],]-BM[s[1],] )
        else 
          up=UpDate(MTT,MBB,MTB,TF[s[1],],TF[s[2],], BF[s[1]], BF[s[2]])
        MTT=up$MTT
        MBB=up$MBB
        MTB=up$MTB
        TF[c(s[1],s[2]),]=TF[c(s[2],s[1]),]  
        TM[c(s[1],s[2]),]=TM[c(s[2],s[1]),]  
      } #end of jumps
    }
    globTF
  } 
  # ******************************************************************************************************************************************************** 
  # Random swaps
  # ********************************************************************************************************************************************************    
  Swaps=function(TF,MF,BF,restrict,pivot,rank,nunits) {
    candidates=NULL
    while (is.null(candidates)) {
      if (rank<(nunits-1)) s1=sample(pivot[(1+rank):nunits],1) else s1=pivot[nunits]
      available=!apply( sapply(1:ncol(TF),function(i) {TF[,i]==TF[s1,i]}),1,all)
      candidates = (1:nunits)[ MF==MF[s1] & restrict==restrict[s1] & BF!=BF[s1] & available==TRUE]
    }
    if ( length(candidates)>1 ) s2=sample(candidates,1) else s2=candidates[1] 
    return(c(s1,s2))
  }
  # ******************************************************************************************************************************************************** 
  # Initial randomized starting design. If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************  
  NonSingular=function(TF,MF,BF,TM,BM,restrict) {
    fullrank=ncol(BM)+ncol(TM)
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    times=0
    while (rank<fullrank && times<1000) {
      times=times+1
      s=Swaps(TF,MF,BF,restrict,pivot,rank,nunits)
      tindex=seq_len(length(BF))
      tindex[c(s[1],s[2])]=tindex[c(s[2],s[1])]
      newQ=qr(t(cbind(BM,TM[tindex,])))
      if (newQ$rank>rank) {   
        TF=TF[tindex,,drop=FALSE]
        TM=TM[tindex,,drop=FALSE]
        rank=newQ$rank
        pivot=newQ$pivot
      }
    }
    if (times>999) stop(" Unable to find an initial non-singular choice of treatment design for this choice of block design") 
    list(TF=TF,TM=TM)
  }  
  # *******************************************************************************************************************************************************
  # Information matrix for single treatment factor design 
  # ******************************************************************************************************************************************************** 
  infoMatSimple=function(TM,BM,MF,BF) { 
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    indicv=seq( (ncol(BM)+1), (ncol(BM)+ncol(TM)))  
    MTT=matrix(0,nrow=(ncol(TM)+1),ncol=(ncol(TM)+1))
    MTT[seq_len(ncol(TM)),seq_len(ncol(TM))]=V[indicv,indicv,drop=FALSE]
    MBB=matrix(0,nrow=nlevels(BF),ncol=nlevels(BF))
    MBB[seq_len(ncol(BM)),seq_len(ncol(BM))]=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
    MTB=matrix(0,nrow=(ncol(TM)+1),ncol=nlevels(BF))
    MTB[seq_len(ncol(TM)),seq_len(ncol(BM))]=V[indicv,seq_len(ncol(BM)),drop=FALSE]
    reorder=c(rbind( matrix(seq_len(ncol(BM)),nrow=(nlevels(BF)/nlevels(MF))-1,ncol=nlevels(MF)),seq(ncol(BM)+1,nlevels(BF))))
    MBB=MBB[reorder,reorder] 
    MTB=MTB[,reorder]
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }
  # *******************************************************************************************************************************************************
  # Information matrix for complex treatment factor design 
  # ******************************************************************************************************************************************************** 
  infoMatComplex=function(TM,BM) { 
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    indicv=seq( (ncol(BM)+1), (ncol(BM)+ncol(TM)))  
    MBB=V[seq_len(ncol(BM)),seq_len(ncol(BM)),drop=FALSE]
    MTB=V[indicv,seq_len(ncol(BM)),drop=FALSE]
    MTT=V[indicv,indicv,drop=FALSE] 
    list(MTT=MTT,MBB=MBB,MTB=MTB)
  }
  # *******************************************************************************************************************************************************
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design. 
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ******************************************************************************************************************************************************** 
  blocksOpt=function(TF,MF,BF,restrict) { 
    BM=Contrasts(MF,BF)[, -seq( nlevels(BF)/nlevels(MF), nlevels(BF) , by=nlevels(BF)/nlevels(MF) ) ,drop=FALSE]
    TM=model.matrix(as.formula(model),TF)[,-1,drop=FALSE] # drops mean contrast
    TM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(TM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    nonsing=NonSingular(TF,MF,BF,TM,BM,restrict)
    TF=nonsing$TF
    TM=nonsing$TM
    V=chol2inv(chol(crossprod(cbind(BM,TM))))
    indicv=seq( (ncol(BM)+1), (ncol(BM)+ncol(TM))  )
    if (simpleTF) inf=infoMatSimple(TM,BM,MF,BF) else inf=infoMatComplex(TM,BM)
    TF=Optimise(TF,MF,BF,restrict,inf$MTT,inf$MBB,inf$MTB,TM,BM)
    TF
  }  
  # ******************************************************************************************************************************************************** 
  # Contrasts for factor NF centered within the levels of factor MF to ensure that NF information is estimated within the levels of factor MF only  
  # ********************************************************************************************************************************************************
  Contrasts=function(MF,NF) {
    NM=matrix(0,nrow=length(NF),ncol=nlevels(NF))
    NM[cbind(seq_len(length(NF)),NF)]=1 # factor indicator matrix  
    NM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(NM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
  }
  # ******************************************************************************************************************************************************** 
  # Finds row and column sizes in each stratum of a design 
  # ********************************************************************************************************************************************************  
  Sizes=function(blocksizes,stratum) {
    nblocks=length(blocksizes)
    newblocksizes=NULL
    for (j in 1:nblocks) {
      rowsizes=rep(blocksizes[j]%/%rows[stratum],rows[stratum])
      resid=blocksizes[j]-sum(rowsizes)
      if (resid>0) 
        rowsizes[1:resid]=rowsizes[1:resid]+1
      rowcolsizes=vector(mode = "list", length =rows[stratum])
      for ( z in 1:rows[stratum])
        rowcolsizes[[z]]=rep(rowsizes[z]%/%columns[stratum] , columns[stratum])
      shift=0
      for (z in seq_len(rows[stratum])) {
        resid=rowsizes[z]-sum(rowcolsizes[[z]])
        if (resid>0) {
          rowcolsizes[[z]][(shift:(shift+resid-1))%%columns[stratum]+1]=rowcolsizes[[z]][(shift:(shift+resid-1))%%columns[stratum]+1]+1
          shift=shift+resid
        }
      }
      newblocksizes=c(newblocksizes,unlist(rowcolsizes))
    }
    newblocksizes
  }
  # ******************************************************************************************************************************************************** 
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  FactEstEffics=function(TF,MF,BF) {
    if (nlevels(MF)==nlevels(BF)) return(1)
    TX=model.matrix(as.formula(model),TF)[,-1,drop=FALSE] # drops mean contrast
    TM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TX[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    BM=Contrasts(MF,BF)[, -seq( nlevels(BF)/nlevels(MF), nlevels(BF) , by=nlevels(BF)/nlevels(MF) ) ,drop=FALSE]
    TB=crossprod(TM,BM)
    RI=backsolve(  chol(crossprod(TM)) ,diag(ncol(TM)))
    QI=backsolve(chol(crossprod(BM)),diag(ncol(BM)))
    U=crossprod(t(crossprod(RI,TB)),QI)
    return(round(det( diag(ncol(TM))-tcrossprod(U))**(1/ncol(TM)),6))
  }
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for row-and-column designs 
  # ********************************************************************************************************************************************************     
  FactRowColEffics=function(Design) {
    TF=Design[, c( (ncol(Design)-ncol(TF)+1):ncol(Design)),drop=FALSE]
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata))  
      for (j in 1:3) 
        effics=c(effics,FactEstEffics(TF,Design[,3*(i-1)+1],Design[,3*(i-1)+1+j]))
    names =unlist(lapply(1:strata, function(j) {c(paste("Rows",j),paste("Columns",j),paste("Rows x Columns",j))}))
    blocks=unlist(lapply(1:strata, function(j) {c( nlevels(Design[,3*(j-1)+2]),nlevels(Design[,3*(j-1)+3]),nlevels(Design[,3*(j-1)+4]))}))
    efficiencies=data.frame(cbind(names,blocks,effics))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies")
    efficiencies
  } 
  # ******************************************************************************************************************************************************** 
  # Calculates D and A-efficiency factors for treatment factors TF assuming block factor BF
  # ********************************************************************************************************************************************************
  FactBlocksEffics=function(Design) { 
    TF=Design[, c( (ncol(Design)-ncol(TF)+1):ncol(Design)),drop=FALSE]
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata)) 
      effics=c(effics,FactEstEffics(TF,Design[,i],Design[,i+1]))
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocks=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocks,effics))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies")
    efficiencies
  }
  # ******************************************************************************************************************************************************** 
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  EstEffics=function(TF,BF) {
    k=nlevels(BF)
    if (k==1) return(c(1,1))
    if (nlevels(TF)<=k) 
      e=eigen( (diag(nlevels(TF))-crossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(nlevels(TF)-1)] else    
        e=c(rep(1,(nlevels(TF)-k)), 
            eigen((diag(k)-tcrossprod(t(table(TF, BF)*(1/sqrt(tabulate(TF))) ) * (1/sqrt(tabulate(BF))))), symmetric=TRUE, only.values = TRUE)$values[1:(k-1)])  
      return(round(c(mean(e)*prod(e/mean(e))^(1/length(e)),1/mean(1/e)),6))
  }
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for block designs 
  # ********************************************************************************************************************************************************   
  BlockEfficiencies=function(Design) {
    effics=matrix(NA,nrow=strata,ncol=2)
    for (i in seq_len(strata))    
      effics[i,]=EstEffics(Design[,ncol(Design)],Design[,i])  
    bounds=rep(NA,strata)
    if (regReps) 
      for (i in seq_len(strata))  
        if (nunits%%nlevels(Design[,i])==0 ) 
          bounds[i]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,i]) )
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocks=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocks,effics,bounds))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
  }
  # ******************************************************************************************************************************************************** 
  # Finds efficiency factors for row-and-column designs 
  # ********************************************************************************************************************************************************     
  RowColEfficiencies=function(Design) {
    Design=Design[,-(ncol(Design)-1)]
    effics=matrix(NA,nrow=(3*strata),ncol=2)
    for (i in 1:strata) {
      effics[3*(i-1)+1,]=EstEffics(Design[,ncol(Design)],Design[,3*(i-1)+1])  
      effics[3*(i-1)+2,]=EstEffics(Design[,ncol(Design)],Design[,3*(i-1)+2])  
      effics[3*(i-1)+3,]=EstEffics(Design[,ncol(Design)],Design[,3*(i-1)+3])  
    }
    bounds=rep(NA,(3*strata))
    if (max(replicates)==min(replicates)) {
      for (i in seq_len(strata))  {
        if (nunits%%nlevels(Design[,3*(i-1)+1])==0)
          bounds[3*(i-1)+1]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+1]))
        if (nunits%%nlevels(Design[,3*(i-1)+2])==0)
          bounds[3*(i-1)+2]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+2]))
        if (nunits%%nlevels(Design[,3*(i-1)+3])==0) 
          bounds[3*(i-1)+3]=upper_bounds(nunits,nlevels(Design[,ncol(Design)]),nlevels(Design[,3*(i-1)+3]))
      }
    }
    names =unlist(lapply(1:strata, function(j) {c(paste("Rows",j),paste("Columns",j),paste("Rows x Columns",j))}))
    blocks=unlist(lapply(1:strata, function(j) {c( nlevels(Design[,3*(j-1)+1]),nlevels(Design[,3*(j-1)+2]),nlevels(Design[,3*(j-1)+3]))}))
    efficiencies=data.frame(cbind(names,blocks,effics,bounds))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
  }
  # ********************************************************************************************************************************************************   
  # Efficiency factors for unreplicated randomized designs 
  # ********************************************************************************************************************************************************     
  UnrepEfficiencies=function() {
    efficiencies=data.frame(cbind("Stratum_1",1,1,1,1))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    efficiencies
  }
  # *******************************************************************************************************************************************************
  # Returns v-1 cyclically generated v x v Latin squares plus the rows array (first) and the columns array (last). If v is prime, the squares are MOLS  
  # ********************************************************************************************************************************************************  
  cMOLS=function(v) {
    mols=lapply(0:(v-1),function(z){do.call(rbind, lapply(0:(v-1), function(j){ (rep(0:(v-1))*z +j)%%v} ))})
    mols[[v+1]]=t(mols[[1]])
    mols=mols[c(1,length(mols),2:(length(mols)-1))]
    mols=lapply(1:length(mols),function(z){ (z-1)*v+mols[[z]] })
    mols
  }
  # *******************************************************************************************************************************************************
  # Tests for and constructs  balanced lattice designs
  # ******************************************************************************************************************************************************** 
  lattice=function(v,r) { 
    TF=vector(length=r*r*v)
    if ( r<4 | (isPrime(v) && r<(v+2)) ) {
      TF=rep(seq_len(v*v),r)[order(unlist(cMOLS(v))[1:(r*v*v)])]
    } else if (r<(v+2)  && (v*v)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      TF[1:(2*v*v)]=c(seq_len(v*v),seq_len(v*v)[order(rep(0:(v-1),v))]   )
      if (r>2) {
        index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(v*v))
        mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])	
        for (i in 1:(r-2)  )
          TF[ c( (v*v*(i-1)+1) : (v*v*i)) + 2*v*v  ] = seq_len(v*v)[order(as.numeric(mols[,,i]))]
      }
    } else if (v==10  && r==4) {
      square1=c(1, 8, 9, 4, 0, 6, 7, 2, 3, 5, 8, 9, 1, 0, 3, 4, 5, 6, 7, 2, 9, 5, 0, 7, 1, 2, 8, 3, 4, 6, 2, 0, 4, 5, 6, 8, 9, 7, 1, 3, 0, 1, 2, 3, 8, 9, 6, 4, 5, 7, 
                5, 6, 7, 8, 9, 3, 0, 1, 2, 4, 3, 4, 8, 9, 7, 0, 2, 5, 6, 1, 6, 2, 5, 1, 4, 7, 3, 8, 9, 0, 4, 7, 3, 6, 2, 5, 1, 0, 8, 9, 7, 3, 6, 2, 5, 1, 4, 9, 0, 8)
      square2=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 3, 0, 4, 9, 6, 7, 2, 1, 8, 5, 5, 4, 8, 6, 7, 3, 0, 2, 1, 9, 4, 1, 6, 7, 0, 5, 9, 3, 2, 8, 2, 6, 7, 5, 9, 8, 4, 0, 3, 1, 
                6, 7, 9, 8, 1, 4, 3, 5, 0, 2, 7, 8, 1, 2, 4, 0, 6, 9, 5, 3, 8, 9, 5, 0, 3, 2, 1, 4, 6, 7, 9, 5, 0, 3, 2, 1, 8, 6, 7, 4, 0, 3, 2, 1, 8, 9, 5, 7, 4, 6)
      TF=c(seq_len(100),seq_len(100)[order(rep(0:9,10))],seq_len(100)[order(square1)],seq_len(100)[order(square2)])
    }
    TF=data.frame(factor(TF))
  }
  # *******************************************************************************************************************************************************
  # Tests for balanced trojan designs and constructs available designs
  # ******************************************************************************************************************************************************** 
  trojan=function(r,k) { 
    TF=vector(length=r*r*k)
    if (isPrime(r)) { 
      for (z in 1:k)
        for (y in 0:(r-1)) 
          for (x in 0:(r-1)) 
            TF[(x + y*r)*k + z]=(y+x*z)%%r + (z-1)*r +1
    } else if ((r*r)%in% c(16,64,256,1024,4096,16384,81,729,6561,625,2401)) {
      index=which(c(16,64,256,1024,4096,16384,81,729,6561,625,2401)==(r*r))
      mols=crossdes::MOLS(c(2,2,2,2,2,2,3,3,3,5,7)[index],c(2,3,4,5,6,7,2,3,4,2,2)[index])	
      for (i in 1:r)
        for (j in 1:r)
          for (x in 1:k)
            TF[x+(j-1)*k+(i-1)*k*r]=mols[i,j,x]+(x-1)*r
    }
    TF=data.frame(factor(TF))
  }
  # ******************************************************************************************************************************************************** 
  # Design data frames for rows columns and row.column blocks
  # ********************************************************************************************************************************************************     
  dataframesBlocks=function(cumrows,cumblocks,rows,columns,strata) {
    rowdesign=data.frame(lapply(1:strata,function(i) {gl(cumrows[i],cumrows[strata]/cumrows[i],labels=unlist(lapply(1:cumrows[i], function(j) {paste0("Block_",j)})))}))
    colnames(rowdesign)=unlist(lapply(1:ncol(rowdesign), function(j) {paste0("Stratum_",j)}))
    list(blkdesign=NULL,rowdesign=rowdesign,coldesign=NULL)
  }
  # ******************************************************************************************************************************************************** 
  # Design data frames for rows columns and row.column blocks
  # ********************************************************************************************************************************************************     
  dataframesRowCol=function(cumrows,cumblocks,rows,columns,strata) {
    blkdesign=as.data.frame(lapply(1:(strata+1),function(i) {gl(cumblocks[i],cumblocks[strata+1]/cumblocks[i],labels=unlist(lapply(1:cumblocks[i], function(j) {paste0("Block_",j)})))}))
    colnames(blkdesign)=unlist(lapply(0:(ncol(blkdesign)-1), function(j) {paste0("Stratum_",j,":Blocks")}))
    rowdesign=data.frame(lapply(1:strata,function(i) {gl(rows[i],cumblocks[strata+1]/cumblocks[i]/rows[i],cumblocks[strata+1],
                                                         labels=unlist(lapply(1:rows[i], function(j) {paste0("Row_",j)}))) }))
    coldesign=data.frame(lapply(1:strata,function(i) {gl(columns[i],cumblocks[strata+1]/cumblocks[i+1],cumblocks[strata+1],
                                                         labels=unlist(lapply(1:columns[i], function(j) {paste0("Col_",j)}))) }))
    rowdesign=data.frame(lapply(1:ncol(rowdesign), function(i){ interaction(blkdesign[,i], rowdesign[,i], sep = ":", lex.order = TRUE) }))
    coldesign=data.frame(lapply(1:ncol(coldesign), function(i){ interaction(blkdesign[,i], coldesign[,i], sep = ":", lex.order = TRUE) }))
    colnames(rowdesign)=unlist(lapply(1:ncol(rowdesign), function(j) {paste0("Stratum_",j,":Rows")}))
    colnames(coldesign)=unlist(lapply(1:ncol(coldesign), function(j) {paste0("Stratum_",j,":Cols")}))
    list(blkdesign=blkdesign,rowdesign=rowdesign,coldesign=coldesign)
  }
  # ******************************************************************************************************************************************************** 
  # Main body of rows design function which tests inputs, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************     
  if (missing(treatments)|is.null(treatments)) stop(" Treatments missing or not defined ") 
  if (anyNA(replicates)|any(is.nan(replicates))|any(!is.finite(replicates))|any(replicates%%1!=0)|any(replicates<1)|is.null(replicates)) stop(" replicates invalid") 
  
  if (!is.data.frame(treatments)){
    if (anyNA(treatments)|any(is.nan(treatments))|any(!is.finite(treatments))|any(treatments%%1!=0)|any(treatments<1)) stop(" treatments parameter invalid") 
    if (length(replicates)!=length(treatments)) stop("treatments and replicates parameters must both be the same length")
    hcf=HCF(replicates)
    treatments=data.frame( Treatments=as.factor(rep(rep(1:sum(treatments), rep(replicates/hcf,treatments)), hcf)))
  } else {
    hcf=1
    treatments=treatments[rep(seq_len(nrow(treatments)), replicates), ,drop=FALSE]
  }
  if (is.null(rows)) rows=hcf
  if (anyNA(rows)|any(is.nan(rows))|any(!is.finite(rows))|any(rows%%1!=0)|any(rows<1)|is.null(rows)) stop(" rows invalid") 
  if (is.null(columns)) columns=rep(1,length(rows))
  if (anyNA(columns)|any(is.nan(columns))|any(!is.finite(columns))|any(columns%%1!=0)|any(columns<1)) stop(" columns parameter invalid") 
  if (length(columns)!=length(rows)) stop("rows and columns vectors must be the same length ")
  if (max(rows*columns)==1) {rows=1; columns=1} else {
    rows=rows[rows*columns>1] 
    columns=columns[rows*columns>1] 
  }
  if (is.na(seed) | !is.finite(seed) | is.nan(seed) | seed%%1!=0 | seed<0 ) stop(" seed parameter invalid  ") 
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) stop(" number of jumps parameter is invalid (max is 10) ") 
  if (!is.null(searches) && ( is.na(searches) | !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 )) stop(" number of searches parameter is invalid") 
  if (max(rows*columns)==1) { rows=1; columns=1} else {index=rows*columns>1; rows=rows[index]; columns=columns[index]}

  for (i in 1:ncol(treatments))
    if (isTRUE(all.equal(treatments[,i], rep(treatments[1,i], length(treatments[,i]))))) stop("One or more treatment factors is a constant which is not valid")
  fnames=colnames(treatments)
  strata=length(rows)
  cumrows=cumprod(rows)
  cumcols=cumprod(columns)
  cumblocks=c(1,cumprod(rows*columns))
  isrowcol=(max(rows)>1 & max(columns)>1)
  # omit any single replicate treatments for unstructured factorial designs and find hcf for factor replicates 
  fulltreatments=treatments
  fullreplicates=replicates
  simpleTF=(ncol(treatments)==1 & is.factor(treatments[,1]))
  if (simpleTF && min(replicates)==1 && max(replicates)>1) { 
    replicates=as.numeric(table((treatments[,1])))
    treatments=droplevels(treatments[ replicates[treatments[,1]]>1 ,1,drop=FALSE])
    replicates=replicates[replicates>1]
    hcf=HCF(replicates)
  }
  nunits=nrow(treatments)
  regReps=isTRUE(all.equal(max(replicates), min(replicates))) 
  
  # default model formula
  modelnames=function (treatments) {unlist(lapply(1:ncol(treatments), function(i) {
    if (!is.factor(treatments[,i])) paste0("poly(",colnames(treatments)[i],",",length(unique(treatments[,i]))-1,")") else colnames(treatments)[i]})) }
  if (is.null(model)) model=paste0("~ ",paste0(modelnames(treatments), collapse="*"))
  # tests for viable design sizes
  
  if (cumrows[strata]*2>nunits) stop("Too many row blocks for the available plots  - every row block must contain at least two plots")
  if (cumcols[strata]*2>nunits) stop("Too many column blocks for the available plots  - every column block must contain at least two plots")
  if (cumblocks[strata+1]>nunits && cumrows[strata]>1 && cumcols[strata]>1) stop("Too many blocks - every row-by-column intersection must contain at least one plot")
  trtdf=ncol(model.matrix(as.formula(model),treatments))-1
  blockdf=sum( (rows+columns-2)*cumblocks[-length(cumblocks)] )
  if ( (trtdf+blockdf+2) > nunits) stop(paste("Too many model parameters: plots df = ", nunits-1," treatments df = ",trtdf," blocks df = ", blockdf ))
  if (is.null(searches)) searches=1+10000%/%nunits
  Plots=factor(1:nunits)
  set.seed(seed)

  blocksizes=nunits
  for (i in 1:strata) 
    blocksizes=Sizes(blocksizes,i)
  regBlocks=isTRUE(all.equal(max(blocksizes), min(blocksizes)))
  
  if ( simpleTF  && isrowcol && min(fullreplicates)==1 && max(fullreplicates)>1 && regBlocks==FALSE )  
    stop("The algorithm does not deal with irregular row-and-column designs containing single replicate treatments ")
  if (isrowcol) 
    df1=dataframesRowCol(cumrows,cumblocks,rows,columns,strata) else 
      df1=dataframesBlocks(cumrows,cumblocks,rows,columns,strata)
  blkdesign=df1$blkdesign
  rowdesign=df1$rowdesign
  coldesign=df1$coldesign
  Mean=factor(rep(1,sum(blocksizes)))
  rowDesign=data.frame(Mean, rowdesign[rep(1:length(blocksizes),blocksizes),]) 
  if (isrowcol) colDesign=data.frame(Mean, coldesign[rep(1:length(blocksizes),  blocksizes ),])
  if (isrowcol) blkDesign=data.frame(blkdesign[rep(1:length(blocksizes),  blocksizes ),])  
  TF=NULL
  orthoSize=nunits/hcf
  # check for algebraic solution
  if (simpleTF) {
    v=sqrt(nlevels(treatments[,1]))  # dimension of a lattice square
    k=nunits/cumblocks[strata+1]  # block size
    orthoMain=(regReps & (replicates[1]==rows[1]))
    # for s orthogonal Latin squares of dimension r x r there are r x kr Trojan designs for r replicates of kr treatments in blocks of size k where k<=s
    if (regReps && regBlocks && orthoMain && !isrowcol && identical(v,floor(v)) & identical(k,v) && length(rows)==2)
      TF=lattice(v,replicates[1])
    else if (regReps && regBlocks && orthoMain && isrowcol && identical(columns[1],replicates[1]) && length(rows)==1 && length(columns)==1 && (k<replicates[1])) #?? identical(rows[1],r)
      TF=trojan(replicates[1],k)
    if (all(TF==FALSE)) TF=NULL
  }
  attempts=0
  while (is.null(TF) & attempts<10) {
    attempts=attempts+1
      TF=data.frame(do.call(rbind,lapply(1:hcf,function(i) {treatments[sample((1+(i-1)*orthoSize):(i*orthoSize)), ,drop=FALSE]}))) # randomize
      colnames(TF)=fnames
      for ( i in 1:strata) {
        if (!isrowcol && rows[i]>1)    TF=blocksOpt(TF,rowDesign[,i],rowDesign[,i+1],rowDesign[,i])
        if ( isrowcol && rows[i]>1)    TF=blocksOpt(TF,blkDesign[,i],rowDesign[,i+1],blkDesign[,i])
        if ( isrowcol && columns[i]>1) TF=blocksOpt(TF,blkDesign[,i],colDesign[,i+1],rowDesign[,i+1])
      }
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  colnames(TF)=fnames
  # add back single rep treatments for nested stratum blocks only
  if (nrow(fulltreatments)!=nrow(treatments) ) {
    treatments=fulltreatments
    nunits=nrow(treatments)
    fullblocksizes=nunits
    for (i in 1:strata)
      fullblocksizes=Sizes(fullblocksizes,i)
    TrtsInBlocks= split( levels(TF[,1])[TF[,1]] , rep(1:length(blocksizes),blocksizes))
    fullreps=as.numeric(table((fulltreatments[,1])))
    singleTF=split( (1:length(fullreps))[fullreps==1], rep(1:length(fullblocksizes),(fullblocksizes-blocksizes)))
    
    for (i in names(singleTF)) 
      TrtsInBlocks[[i]]=sample(append( TrtsInBlocks[[i]] ,  singleTF[[i]]))
    TF=data.frame(unlist(TrtsInBlocks))
    blocksizes=fullblocksizes
    Plots=factor(1:nunits)
    Stratum_0=as.factor(rep(1,nunits))
    rowDesign=data.frame(Stratum_0,rowdesign[rep( 1:length(blocksizes),  blocksizes ),]) 
  }
  
  if (!isrowcol) {
    rowDesign=cbind(rowDesign,Plots)
    rowDesign=data.frame(lapply(1:ncol(rowDesign), function(r){sample(nlevels(rowDesign[,r]))[rowDesign[,r]]})) # Randomize labels - NB gives numeric columns
    rowDesign=cbind(rowDesign,TF)
    rowDesign=rowDesign[do.call(order, rowDesign), ] # re-order
    blocksizes=table(rowDesign[,ncol(rowDesign)-ncol(TF)-1])[unique(rowDesign[,ncol(rowDesign)-ncol(TF)-1])]
    TF=rowDesign[,c((ncol(rowDesign)-ncol(TF)+1):ncol(rowDesign)),drop=FALSE]
    for (i in 1 : ncol(treatments)) 
      if (is.factor(treatments[,i])) TF[,i]=as.factor(TF[,i])
    Design  = data.frame(rowdesign[rep(1:length(blocksizes),blocksizes ),],Plots, TF)
    rownames(Design)=NULL
    colnames(Design)=c(colnames(rowdesign),"Plots",fnames)
    V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))],Design[,(ncol(Design)-ncol(TF)-1)])
    V=lapply(V, function(x){ length(x) =max(blocksizes); x })
    PlotsInBlocks=rep("",length(V))
    if (ncol(treatments)==1) {
      Plan=data.frame(rowdesign,PlotsInBlocks,matrix( unlist(V), nrow=length(V),byrow=TRUE))
      colnames(Plan)= c(colnames(Plan[,1:(strata+1)]), c(1:(ncol(Plan)-strata-1)))
      rownames(Plan)=NULL
    }
    else Plan=NULL
  }
  
  if (isrowcol) {
    rdf = data.frame(rowdesign,coldesign,blkdesign[,2:ncol(blkdesign),drop=FALSE])[c(rbind(1:strata,(strata+1):(2*strata),(2*strata+1):(3*strata)))]
    rcDesign=data.frame(rdf[rep(1:length(blocksizes), blocksizes),],Plots)
    rcDesign=data.frame(lapply(1:ncol(rcDesign), function(r){ sample(nlevels(rcDesign[,r]))[rcDesign[,r]]})  ) # Randomize
    rcDesign=data.frame(rcDesign, TF)
    rcDesign=rcDesign[do.call(order,rcDesign), ] # re-order
    TF=rcDesign[,c((ncol(rcDesign)-ncol(TF)+1):ncol(rcDesign)),drop=FALSE] 
    blocksizes=table(rcDesign[,ncol(rcDesign)-ncol(TF)-1])[unique(rcDesign[,ncol(rcDesign)-ncol(TF)-1])]
    Design  = data.frame( rdf[rep(1:length(blocksizes), blocksizes),],Plots,TF)  # rebuild factor levels
    colnames(Design)=c(colnames(rdf),"Plots",fnames)
    rownames(Design)=NULL
    if (max(blocksizes)>1 && ncol(treatments)==1)   {
      V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))],Design[,(ncol(Design)-ncol(TF)-1)]) # split on blocks
      V=lapply(V, function(x){ length(x) =max(blocksizes); x })
      PlotsInBlocks=rep("",length(V))
      Plan=data.frame(rdf,PlotsInBlocks,matrix( unlist(V), nrow=length(V),byrow=TRUE))
      colnames(Plan)=c(colnames(rdf),"Plots_In_Blocks",1:max(blocksizes))
      Plan=Plan[,-(ncol(Plan)-max(blocksizes)-1)]
      sets=nlevels(Plan[,ncol(rdf)-3])
      sizes=nrow(Plan)/sets
    } else if (max(blocksizes)==1 && ncol(treatments)==1) {
      V=split(Design[,c((ncol(Design)-ncol(TF)+1):ncol(Design))], Design[,(ncol(Design)-ncol(TF)-3)] ) # split on rows
      plan = rdf[seq(1,nrow(rdf)-columns[strata]+1,columns[strata]),-c(ncol(rdf)-1,ncol(rdf)) ,drop=FALSE]
      Columns=rep("",nrow(plan))
      Plan=data.frame(plan,Columns,matrix( unlist(V), nrow=length(V),byrow=TRUE))
      colnames(Plan)=c(colnames(plan),paste0("Stratum_",strata,".Cols"),1:columns[strata])
      rownames(Plan)=NULL
      sets=nlevels(Plan[,ncol(plan)-1])
      sizes=nrow(Plan)/sets
    } else Plan=NULL
    
    if (!is.null(Plan) && (strata>1)) {
      reorder=order(c(1:nrow(Plan)+rep(0:(sets-1),each=sizes),(1:sets)*(sizes+1)))
      Plan=rbind(Plan, setNames(data.frame( matrix(" ", nrow=sets,ncol=ncol(Plan)) ), names(Plan)))[reorder,]
    }
  }
  
  # efficiencies

  if (simpleTF && isrowcol) Efficiencies=RowColEfficiencies(Design)
  else if (simpleTF && max(replicates)==1) Efficiencies=UnrepEfficiencies()
  else if (simpleTF)  Efficiencies=BlockEfficiencies(Design)
  else if ( isrowcol) Efficiencies=FactRowColEffics(Design)
  else if ( !isrowcol) Efficiencies=FactBlocksEffics(Design)
  row.names(Efficiencies)=NULL
  # treatment replications
  if(simpleTF) {
  TreatmentsTable=data.frame(table(Design[,ncol(Design)]))
  TreatmentsTable=TreatmentsTable[order( as.numeric(levels(TreatmentsTable[,1])) [TreatmentsTable[,1]])           ,]
  TreatmentsTable[]=lapply(TreatmentsTable, as.factor) 
  colnames(TreatmentsTable)=c("Treatments","Replicates")
  rownames(TreatmentsTable)=1:nrow(TreatmentsTable)
  } else {
    TreatmentsTable=data.frame(treatments,rep(replicates,nrow(treatments)))
    colnames(TreatmentsTable)=c(colnames(treatments),"Replicates")
  }
  
  list(Treatments=TreatmentsTable,model=model,Efficiencies=Efficiencies,Plan=Plan,Design=Design,seed=seed,searches=searches,jumps=jumps) 
} 