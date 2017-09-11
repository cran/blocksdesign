#' @title Block designs for factorial treatment sets
#'
#' @description
#'
#' Constructs randomized nested block designs for factorial or fractional factorial treatment designs with any
#' feasible depth of nesting and up to two crossed block structures in each level of nesting.
#'
#' @details
#'
#' \code{factblocks} generates blocked factorial designs for general factorial treatment structures possibly including
#' mixtures of qualitative and quantitative level factors. Qualitative level factors are
#' modelled factorially while quantitative level factors are modelled by polynomials of the required degree.  
#' Designs can be based on any multiple, not necessarily integral, of the complete factorial 
#' treatment design where the fractional part of the design, if any, is chosen by optimizing a 
#' D-optimal fraction of that size for that treatment design.
#'
#' The \code{treatments} parameter defines the treatment factors of the design and must be be a data frame with
#' a column for each factor and a row for each factorial combination (see examples). The treatment factors 
#' can be any mixture of qualitative or quantitative level factors and the treatment model can be any feasible model defined 
#' by the \code{models} formula of the \code{\link[stats]{model.matrix}} package (see examples). 
#' 
#' Quantitative factors can be modelled either by raw or by orthogonal polynomials. Orthogonal polynomials are numerically more stable
#' than raw polynomials and are usually the best choice at the design stage. Polynomial models can be fitted at the analysis stage either by raw or
#' by orthogonal polynomials regardless of the type of polynomial fitted at the design stage.
#'   
#' The \code{replicates} parameter defines the required replication for the treatments design and should be a single number, not necessarily integral, 
#' representing a required multiple or a required fraction of the \code{treatments} data frame. The algorithm will find a
#' D-optimal or near D-optimal fraction of the required size for any fractional part of replication number, assuming the required design is non-singular.
#' 
#' The \code{rows} parameter, if any, defines the nested row blocks for each level of nesting taken in order from the highest to the lowest. The
#' first number, if any, is the number of nested row blocks in the first-level of nesting, the second number, if any, is the number of nested row blocks in
#' the second-level of nesting and so on for any required feasible depth of nesting.
#' 
#' The \code{columns} parameter, if any, defines the nested column blocks for each level of nesting taken in order from the highest to the lowest.
#' The first number, if any, is the number of nested column blocks in the first-level of nesting, the second, if any, is the number of nested column blocks in
#' the second-level of nesting and so on for the same required depth of nesting as in the \code{rows} parameter.
#' 
#' The \code{rows} and \code{columns} parameters, if defined, must be of equal length. If a simple set of nested blocks is required for
#' any particular level of nesting, the number of columns for that level should be set to unity. Any required combination of simple or
#' crossed blocks can be obtained by appropriate choice of the levels of the \code{rows} and \code{columns} parameters.
#' If the \code{rows} parameter is defined but the \code{columns} parameter is null, the design will be a simple nested
#' blocks design with numbers of block levels defined by the \code{rows} parameter. If both the \code{rows} parameter and the \code{columns} parameter are null, 
#' the default block design will be a set of orthogonal main blocks equal in number to the highest common factor of the replication numbers. 
#'
#' Block sizes are always as nearly equal as possible and will never differ by more than a single plot for any particular block classification. 
#' Row blocks and column blocks must always contain at least two plots per block and this restriction will constrain the permitted numbers of 
#' rows and columns in the various nested levels of a block design.
#'
#' For any particular level of nesting, the algorithm first optimizes the row blocks conditional on any higher-level blocks
#' and then optimizes the columns blocks, if any, conditional on the rows blocks.
#' 
#' The efficiency factor of a fractional factorial design is the generalized variance of the complete factorial design divided by the generalized variance of
#' the fractional factorial design where the generalized variance of a design is the (1/p)th power of the determinant of the crossed-product of the p-dimensional
#' model matrix divided by the number of observations in the design. 
#'
#' Comment:
#'
#' Row-and-column designs may contain useful treatment information in the individual row-by-column intersection blocks but \code{blocksdesign} does not currently
#' optimize the efficiency of these blocks.
#'
#' Row-and-column design with 2 complete treatment replicates, 2 complete rows and 2 complete columns will always confound one treatment contrast in the
#' rows-by-columns interaction. For these designs, it is impossible to nest a non-singular block design in the rows-by-columns intersections and instead
#' we suggest a randomized nested blocks design with four incomplete main blocks.
#'
#' Outputs:
#'
#' The principle design outputs comprise:
#' \itemize{
#'  \item  A data frame showing the allocation of treatments to blocks with successive nested strata arranged in standard block order. \cr
#'  \item  A table showing the replication number of each treatment in the design. \cr
#'  \item  An efficiency factor for fractional factorial treatment designs. \cr
#'  \item  A table showing the block levels and the achieved D-efficiency factors for each stratum. \cr
#' }
#'
#' @param treatments  a data frame with columns for individual treatment factors and rows for individual treatment factor combinations.
#'
#' @param replicates  a single replication number, not necessarily integral.
#'
#' @param rows the number of rows nested in each preceding block for each level of nesting from the top-level block downwards. The top-level block is a
#' single super-block which need not be defined explicitly. 
#' 
#' @param columns the number of columns nested in each preceding block for each level of nesting from the top-level block downwards. The \code{rows} and 
#' \code{columns} parameters must be of equal length unless the \code{columns} parameter is null, in which case the design has a
#' single column block for each level of nesting and the design becomes a simple nested row blocks design.
#'
#' @param model a model equation for the treatment factors in the design where the equation is defined using the model.matrix notation
#' in the {\link[stats]{model.matrix}} package. If undefined, the model is a full factorial treatment design.
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
#' \item{BlocksEfficiency}{The D-efficiencies of the blocks in each stratum of the design.}
#' \item{DesignEfficiency}{The generalized variance of the complete factorial design divided by the generalized variance of the fractional factorial design.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#'
#'
#' @references
#'
#' Sailer, M. O. (2013). crossdes: Construction of Crossover Designs. R package version 1.1-1. https://CRAN.R-project.org/package=crossdes
#'
#' Edmondson R. N. (1998). Trojan square and incomplete Trojan square designs for crop research. Journal of Agricultural Science, Cambridge, 131, pp.135-142
#'
#' Cochran, W.G., and G.M. Cox. 1957. Experimental Designs, 2nd ed., Wiley, New York.
#'
#'
#' @examples
#' 
#' ## The number of searches in the examples have been limited for fast execution. 
#' ## For optimum results, the number of searches may need to be increased in practice.
#' ## Designs should be rebuilt repeatedly to check that a near-optimum design has been found.  
#' 
#' 
#' ## Factorial designs defined by a treatments data frame and a factorial model equation.
#' 
#' # Main effects of five 2-level factors in a half-fraction of a 4 x 4 row-and column design.
#' GF = expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' factblocks(treatments=GF,model="~ F1+F2+F3+F4+F5",replicates=.5,rows=4,columns=4,searches=20)
#' 
#' # Quadratic regression for one 6-level numeric factor in 2 randomized blocks assuming 10/6 fraction
#' factblocks(treatments=expand.grid(X=1:6),model=" ~ poly(X,2)",rows=2,searches=5,replicates=10/6) 
#' 
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2))
#' factblocks(treatments=GF,model=" ~ (F1+F2+F3+F4+F5)^2",rows=4,searches=5)
#' 
#' # First-order model for 1/3rd fraction of four qualitative 3-level factors in 3  blocks
#' GF=expand.grid(F1=factor(1:3),F2=factor(1:3),F3=factor(1:3),F4=factor(1:3))
#' factblocks(treatments=GF,model=" ~ F1+F2+F3+F4",replicates=(1/3),rows=3,searches=5)
#' 
#' # Second-order model for a 1/3rd fraction of five qualitative 3-level factors in 3 blocks
#' GF=expand.grid( F1=factor(1:3), F2=factor(1:3), F3=factor(1:3), F4=factor(1:3), F5=factor(1:3) )
#' factblocks(treatments=GF,model=" ~ (F1+F2+F3+F4+F5)^2",rows=3,replicates=(1/3),searches=5)
#' 
#' # Second-order model for two qualitative and two quantitative level factors in 4 randomized blocks
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:3),V1=1:3,V2=1:4)
#' modelform=" ~ F1 + F2 + poly(V1,2) +  poly(V2,2) + (poly(V1,1)+F1+F2):(poly(V2,1)+F1+F2) "
#' \dontrun{factblocks(treatments=GF,model=modelform,rows=4,searches=5)}
#' 
#' # Plackett and Burman design for eleven 2-level factors in 12 runs (needs large number of searches)
#' GF=expand.grid(F1=factor(1:2),F2=factor(1:2),F3=factor(1:2),F4=factor(1:2),F5=factor(1:2),
#' F6=factor(1:2),F7=factor(1:2),F8=factor(1:2),F9=factor(1:2),F10=factor(1:2),F11=factor(1:2))
#' \dontrun{factblocks(GF,model="~ F1+F2+F3+F4+F5+F6+F7+F8+F9+F10+F11",replicates=(12/2048))}
#' 
#'
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames
#'
factblocks = function(treatments,replicates=1,rows=NULL,columns=NULL,model=NULL,searches=NULL,seed=sample(10000,1),jumps=1) {

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
    return(newblocksizes)
  }
  # ********************************************************************************************************************************************************
  # Calculates D and A-efficiency factors for treatment factor TF assuming block factor BF
  # ********************************************************************************************************************************************************
  FactEstEffics=function(TF,MF,BF) {
    if (nlevels(MF)==nlevels(BF)) return(1)
    dd=data.frame(TF,BF)
    TM=model.matrix(as.formula(model),dd)[,-1,drop=FALSE] # drops mean contrast
    BM=model.matrix(~BF,dd)[,-1,drop=FALSE] # drops mean contrast
    TM=do.call(rbind,lapply(1:nlevels(MF),function(i) {scale(TM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    BM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(BM[MF==levels(MF)[i],] , center = TRUE, scale = FALSE)}))
    BM=BM[, -seq( nlevels(BF)/nlevels(MF), nlevels(BF) , by=nlevels(BF)/nlevels(MF) ) ,drop=FALSE]
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
    TF=Design[, c( (ncol(Design)-ncol(treatments)+1):ncol(Design)),drop=FALSE]
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata))
      for (j in 1:3)
        effics=c(effics,FactEstEffics(TF,Design[,3*(i-1)+1],Design[,3*(i-1)+1+j]))
    names =unlist(lapply(1:strata, function(j) {c(paste("Rows",j),paste("Columns",j),paste("Rows x Columns",j))}))
    blocklevs=unlist(lapply(1:strata, function(j) {c( nlevels(Design[,3*(j-1)+2]),nlevels(Design[,3*(j-1)+3]),nlevels(Design[,3*(j-1)+4]))}))
    efficiencies=data.frame(cbind(names,blocklevs,effics))
    colnames(efficiencies)=c("Stratum","Blocks","D-Efficiencies")
    return(efficiencies)
  }
  # ********************************************************************************************************************************************************
  # Calculates D and A-efficiency factors for treatment factors TF assuming block factor BF
  # ********************************************************************************************************************************************************
  FactBlocksEffics=function(Design) {
    TF=Design[, c( (ncol(Design)-ncol(treatments)+1):ncol(Design)),drop=FALSE]
    effics=NULL
    Design=data.frame(as.factor(rep(1,nrow(Design))),Design)
    for (i in seq_len(strata))
      effics=c(effics,FactEstEffics(TF,Design[,i],Design[,i+1]))
    names =unlist(lapply(1:strata, function(j) {paste0("Stratum_",j)}))
    blocklevs=unlist(lapply(2:(strata+1), function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocklevs,effics))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies")
    return(efficiencies)
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
    blocklevs=unlist(lapply(1:strata, function(j) {nlevels(Design[,j])}))
    efficiencies=data.frame(cbind(names,blocklevs,effics,bounds))
    colnames(efficiencies)=c("Strata","Blocks","D-Efficiencies","A-Efficiencies", "A-Bounds")
    return(efficiencies)
  }
     # ******************************************************************************************************************************************************** 
     # Maximises the blocks design matrix using the matrix function dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
     # Sampling is used initially when many feasible swaps are available but later a full search is used to ensure steepest ascent optimization.
     # ********************************************************************************************************************************************************
     DMax=function (VTT,VBB,VTB,TF,MF,TM,BM,restrict)  {
       locrelD=1
       mainSizes=tabulate(restrict)
       nSamp=pmin( rep(8,nlevels(restrict)), mainSizes)
       repeat {
         kmax=1
         for (k in 1: nlevels(restrict)) {
           s=sort(sample( seq_len(nrow(TF)) [restrict==levels(restrict)[k]], nSamp[k])) 
           TMB=crossprod(t(crossprod(t(TM[s,]),VTB)),t(BM[s,]))
           TMT=crossprod(t(crossprod(t(TM[s,]),VTT)),t(TM[s,]))
           BMB=crossprod(t(crossprod(t(BM[s,]),VBB)),t(BM[s,]))
           TMB=sweep(TMB,1,diag(TMB))
           TMT=sweep(TMT,1,diag(TMT))
           BMB=sweep(BMB,1,diag(BMB))
           dMat=(1+TMB+t(TMB))**2 - (TMT + t(TMT))*(BMB + t(BMB))
           sampn=which.max(dMat)
           i=1+(sampn-1)%%nrow(dMat)
           j=1+(sampn-1)%/%nrow(dMat)
           if (dMat[i,j]>kmax) {kmax=dMat[i,j]; pi=s[i]; pj=s[j]} 
         }
         if (kmax>(1+tol)) {
           locrelD=locrelD*kmax
           t=TM[pi,]-TM[pj,]
           b=BM[pj,]-BM[pi,]
           up=UpDate(VTT,VBB,VTB,t,b)
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
 
   # ********************************************************************************************************************************************************
   # Updates variance matrix for pairs of swapped treatments using standard matrix updating formula
   # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
   # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
   # ********************************************************************************************************************************************************
   UpDate=function(VTT,VBB,VTB,t,b) {
     VTTt=crossprod(VTT,t)
     VBBb=crossprod(VBB,b)
     VTBt=crossprod(VTB,t)
     VTBb=crossprod(t(VTB),b)
     tMt=as.numeric(crossprod(t,VTTt))
     bMb=as.numeric(crossprod(b,VBBb))
     tMb=as.numeric(crossprod(b,VTBt))
     f1=(VTTt+VTBb)/sqrt(2)
     f2=(VBBb+VTBt)/sqrt(2)
     g1=(VTBb-VTTt)/sqrt(2)
     g2=(VBBb-VTBt)/sqrt(2)
     a=(tMt+bMb+2*tMb)/2
     b=(tMt+bMb-2*tMb)/2
     c=(bMb-tMt)/2
     d=(1+a)*(1-b)+c*c
     VTT=VTT- (tcrossprod(f1)*(1-b) - tcrossprod(g1)*(1+a) + (tcrossprod(g1,f1)+tcrossprod(f1,g1))*c)/d
     VBB=VBB- (tcrossprod(f2)*(1-b) - tcrossprod(g2)*(1+a) + (tcrossprod(g2,f2)+tcrossprod(f2,g2))*c)/d
     VTB=VTB- (tcrossprod(f1,f2)*(1-b) - tcrossprod(g1,g2)*(1+a) + (tcrossprod(g1,f2)+tcrossprod(f1,g2))*c)/d
     list(VTT=VTT,VBB=VBB,VTB=VTB)
   }
  # ********************************************************************************************************************************************************
  #  Searches for an optimization with selected number of searches and selected number of junps to escape local optima
  # ********************************************************************************************************************************************************
  Optimise=function(TF,MF,fBF,restrict,VTT,VBB,VTB,TM,BM) {
    globrelD=0
    relD=1
    globTF=TF
    blocksizes=tabulate(fBF)
    if (regReps & max(blocksizes)==min(blocksizes) & ncol(treatments)==1 & is.factor(treatments[,1]))
      rowsEffBound=upper_bounds(nrow(TF),nlevels(TF[,1]),nlevels(fBF)) 
    else if ( ncol(treatments)==1 & nlevels(fBF)>1 & is.factor(treatments[,1])) 
      rowsEffBound=1
    else rowsEffBound=NA
    for (r in 1:searches) {
      dmax=DMax(VTT,VBB,VTB,TF,MF,TM,BM,restrict) 
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
          if (!is.na(rowsEffBound)) {
            reff=EstEffics(globTF[,1],fBF)[2]
            if (isTRUE(all.equal(rowsEffBound,reff))) return(globTF)
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
            z= seq_len(nunits)[MF==MF[s1] & restrict==restrict[s1] & fBF!=fBF[s1] & available]   
            if (length(z)==0) next
            if (length(z)>1) s=c(s1,sample(z,1))  else s=c(s1,z)
            TMT=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],VTT)),TM[s[2],]-TM[s[1],])
            BMB=crossprod(t(crossprod(BM[s[1],]-BM[s[2],],VBB)),BM[s[2],]-BM[s[1],])
            TMB=crossprod(t(crossprod(TM[s[1],]-TM[s[2],],VTB)),BM[s[2],]-BM[s[1],] )
            Dswap=(1+TMB)**2-TMT*BMB
            if (Dswap>.1| counter>1000) break
          }
          if (counter>1000) return(globTF) # finish with no non-singular swaps
          relD=relD*Dswap 
          up=UpDate(VTT,VBB,VTB, TM[s[1],]-TM[s[2],], BM[s[2],]-BM[s[1],] )
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
  NonSingular=function(TF,MF,BF,restrict,TM,BM) {
    fullrank=ncol(cbind(TM,BM))
    Q=qr(t(cbind(BM,TM)))
    rank=Q$rank
    pivot=Q$pivot
    for ( i in 1:1000) {
      if (rank==fullrank) return(list(TF=TF,TM=TM))
      s=Swaps(TF,MF,BF,restrict,pivot,rank,nunits)
      tindex=1:nrow(TF)
      tindex[c(s[1],s[2])]=tindex[c(s[2],s[1])]
      Q=qr(t(cbind(BM,TM[tindex,])))
      if (Q$rank>rank) {
        TF=TF[tindex,,drop=FALSE]
        TM=TM[tindex,,drop=FALSE]
        rank=Q$rank
        pivot=Q$pivot
      }
    }
    stop(" Unable to find an initial non-singular choice of treatment design for this choice of block design")
  }
  # *******************************************************************************************************************************************************
  # Optimize the nested Blocks assuming a possible set of Main block constraints Initial randomized starting design.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # ********************************************************************************************************************************************************
  blocksOpt=function(TF,MF,fBF,BF,restrict) {
      TM=model.matrix(as.formula(model),TF)[,-1,drop=FALSE] # drops mean contrast
      TM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(TM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)})) #centres treatments within main blocks
      if (nlevels(MF)==1) BM=model.matrix(as.formula(~BF))[,-1,drop=FALSE] else BM=model.matrix(as.formula(~MF+MF:BF))[,-c(1:nlevels(MF)),drop=FALSE]
      BM=do.call(rbind,lapply(1:length(levels(MF)),function(i) {scale(BM[MF==levels(MF)[i],], center = TRUE,scale = FALSE)})) #centres sub-blocks within main blocks
      BM=BM[ ,as.numeric(matrix(1:(nlevels(BF)*nlevels(MF)),nrow=nlevels(BF),ncol=nlevels(MF),byrow=TRUE)[-nlevels(BF),]) ,drop=FALSE] # reorder within main blocks
      if ( (ncol(cbind(TM,BM))+1) > nrow(TM)  ) stop( paste("Too many parameters: plots df = ",nrow(TM)-1," model:df = ",ncol(cbind(TM,BM)) )) 
      nonsing=NonSingular(TF,MF,fBF,restrict,TM,BM)
      TF=nonsing$TF
      TM=nonsing$TM
      V=chol2inv(chol(crossprod(cbind(TM,BM))))
      VTT = V[1:ncol(TM),1:ncol(TM),drop=FALSE]
      VBB = V[ (ncol(TM)+1):ncol(V) , (ncol(TM)+1):ncol(V),drop=FALSE]
      VTB = V[1:ncol(TM), (ncol(TM)+1):ncol(V), drop=FALSE]
      TF=Optimise(TF,MF,fBF,restrict,VTT,VBB,VTB,TM,BM)
    return(TF)
  }
  # ********************************************************************************************************************************************************
  # Updates variance matrix for swapped rows where mi is swapped out and mj is swapped in
  # ********************************************************************************************************************************************************
  fractUpDate=function(D,mi,mj) {
    f=crossprod(D,mj)
    g=crossprod(D,mi)
    a=as.numeric(crossprod(mj,f))
    b=as.numeric(crossprod(mi,g))
    c=as.numeric(crossprod(mi,f))
    W=g*sqrt(1+a)-f*c/sqrt(1+a)
    d=(1+a)*(1-b)+c*c
    V=f*sqrt(d/(1+a))
    D=D-(tcrossprod(V)-tcrossprod(W))/d
    return(D)
  }
  # ********************************************************************************************************************************************************
  # Fractional factorials
  # ********************************************************************************************************************************************************
  factorial=function(TF,model,replicates) {
    geff=1
    nunits=nrow(TF)*floor(replicates) + floor( nrow(TF)*(replicates%%1) ) # total design number of units 
    if (floor( nrow(TF)*(replicates%%1) )==0 ) {
      gTF=TF[unlist(lapply(1:replicates,function(j){sample(1:nrow(TF))})),,drop=FALSE] 
    } else { 
      allfactors=all(unlist(lapply(TF, class))=="factor") 
      bunits=nrow(TF)*floor(replicates) # base units
      FTF=TF[rep(seq_len(nrow(TF)),ceiling(replicates)), ,drop=FALSE] # complete orthogonal design
      FTM=model.matrix(as.formula(model),FTF) # complete orthogonal model matrix
      if ( ncol(FTM)>nunits) stop("Fractional factorial design too small to estimate all the required model parameters ")
      maxDeff=(det(crossprod(FTM))**(1/ncol(FTM)))/nrow(FTM) 
      gfracDeff=0
      if (is.null(searches)) searches=1+10000%/%nunits
      for (i in 1:searches) {
        counter=0
        repeat {
         rerand=unlist(lapply(1:ceiling(replicates),function(j){(j-1)*nrow(TF)+sample(1:nrow(TF))}))
         fTM=FTM[rerand,,drop=FALSE]
         if (counter>99 | qr(fTM[1:nunits,,drop=FALSE])$rank==ncol(FTM)) break else counter=counter+1
        }
       if (counter>99) stop("Unable to find a non-singular starting design of the required size by random search ")
       fTF=FTF[rerand,,drop=FALSE]
       Info=crossprod(fTM[1:nunits,,drop=FALSE]) # the information matrix of the the non-singular starting design
       V=chol2inv(chol(Info)) # V is the variance matrix of the the non-singular starting design
       repeat {
          M1VM2=       tcrossprod(tcrossprod(fTM[(bunits+1):nunits,,drop=FALSE],V),fTM[(nunits+1):nrow(fTM),,drop=FALSE])
          M1VM1=1-diag(tcrossprod(tcrossprod(fTM[(bunits+1):nunits,,drop=FALSE],V),fTM[(bunits+1):nunits,,drop=FALSE] ))
          M2VM2=1+diag(tcrossprod(tcrossprod(fTM[(nunits+1):nrow(fTM),,drop=FALSE],V),fTM[(nunits+1):nrow(fTM),,drop=FALSE]))
          Z=M1VM2**2 + tcrossprod(M1VM1,M2VM2)
          maxindex=which.max(Z)
          i=1+(maxindex-1)%%nrow(Z)
          j=1+(maxindex-1)%/%nrow(Z)
          if (Z[i,j]<(1+tol)) break
          V=fractUpDate(V ,as.numeric(fTM[bunits+i,]) ,  as.numeric(fTM[nunits+j,]) )# parameters(V,row_swappedout,row_swappedin) 
          fTF[ c(bunits+i,nunits+j), ] =fTF[ c(nunits+j,bunits+i), ]
          fTM[ c(bunits+i,nunits+j), ]= fTM[ c(nunits+j,bunits+i), ]
        }
      fracDeff=det(crossprod(fTM[(1:nunits),,drop=FALSE]))**(1/ncol(fTM))/nunits
      if (fracDeff>gfracDeff) {
        gfracDeff=fracDeff
        geff=gfracDeff/maxDeff
        gTF=fTF[(1:nunits),,drop=FALSE]
        if (geff>(1-tol) & allfactors) break
        }
      }
    }
    return(list(TF=gTF,eff=geff,fraction=nunits/nrow(TF)))
  }
  # ******************************************************************************************************************************************************** 
  # Design data frames for rows columns and row.column blocks
  # ********************************************************************************************************************************************************     
  dataframes=function(rows,columns) {
    nblkdesign=as.data.frame(lapply(1:(strata+1),function(i) {gl(nestblocks[i],cumblocks[strata+1]/cumblocks[i])}))
    for ( i in 1:length(nestblocks)) levels(nblkdesign[,i])=unlist(lapply(1:nestblocks[i],function(j) {paste0("Blocks_",j)})) 
    
    nrowdesign=data.frame(lapply(1:strata,function(i) {gl(rows[i],cumblocks[strata+1]/cumblocks[i]/rows[i],cumblocks[strata+1]) }))
    for ( i in 1:length(rows)) levels(nrowdesign[,i])=unlist(lapply(1:rows[i],function(j)  {paste0("Rows_",j)})) 

    ncoldesign=data.frame(lapply(1:strata,function(i) {gl(columns[i],cumblocks[strata+1]/cumblocks[i+1],cumblocks[strata+1]) }))
    for ( i in 1:length(columns)) levels(ncoldesign[,i])=unlist(lapply(1:columns[i],function(j) {paste0("Cols_",j)}))
    
    fblkdesign=as.data.frame(lapply(1:(strata+1),function(i) {gl(cumblocks[i],cumblocks[strata+1]/cumblocks[i],
                                                             labels=unlist(lapply(1:cumblocks[i], function(j) {paste0("Blocks_",j)})))}))
    frowdesign=data.frame(lapply(1:ncol(nrowdesign), function(i){ interaction(fblkdesign[,i], nrowdesign[,i], sep = ":", lex.order = TRUE) }))
    fcoldesign=data.frame(lapply(1:ncol(ncoldesign), function(i){ interaction(fblkdesign[,i], ncoldesign[,i], sep = ":", lex.order = TRUE) }))
    colnames(fblkdesign)=unlist(lapply(1:ncol(fblkdesign), function(j) {paste0("Level_",j-1)}))
    colnames(frowdesign)=unlist(lapply(1:ncol(frowdesign), function(j) {paste0("Level_",j,":Rows")}))
    colnames(fcoldesign)=unlist(lapply(1:ncol(fcoldesign), function(j) {paste0("Level_",j,":Cols")}))
    colnames(nblkdesign)=unlist(lapply(1:ncol(fblkdesign), function(j) {paste0("Level_",j-1)}))
    colnames(nrowdesign)=unlist(lapply(1:ncol(frowdesign), function(j) {paste0("Level_",j,":Rows")}))
    colnames(ncoldesign)=unlist(lapply(1:ncol(fcoldesign), function(j) {paste0("Level_",j,":Cols")}))
    list(nblkdesign=nblkdesign,nrowdesign=nrowdesign,ncoldesign=ncoldesign,fblkdesign=fblkdesign,frowdesign=frowdesign,fcoldesign=fcoldesign)
  }
  # ********************************************************************************************************************************************************
  # Main design program which tests input variables, omits any single replicate treatments, optimizes design, replaces single replicate
  # treatments, randomizes design and prints design outputs including design plans, incidence matrices and efficiency factors
  # ********************************************************************************************************************************************************
  options(contrasts=c('contr.SAS','contr.poly'))
  tol=.Machine$double.eps^0.5
  if (missing(treatments)|is.null(treatments)|!is.data.frame(treatments) ) stop(" Treatments data frame missing or not defined ")
  for (i in 1:ncol(treatments))
    if (isTRUE(all.equal(treatments[,i], rep(treatments[1,i], length(treatments[,i]))))) stop("One or more treatment factors is a constant which is not valid")
  if (is.null(replicates)|anyNA(replicates)|any(is.nan(replicates))|any(!is.finite(replicates))) stop(" replicates invalid")
  if (is.na(seed) | !is.finite(seed) | is.nan(seed) | seed%%1!=0 | seed<0 ) stop(" seed parameter invalid  ")
    set.seed(seed)
  if (is.na(jumps) | !is.finite(jumps) | is.nan(jumps) | jumps<1 | jumps%%1!=0 | jumps>10) stop(" number of jumps parameter is invalid (max is 10) ")
  if (ceiling(replicates)==floor(replicates)) hcf=replicates else hcf=1
  if (is.null(rows)) rows=hcf
  if (anyNA(rows)|any(is.nan(rows))|any(!is.finite(rows))|any(rows%%1!=0)|any(rows<1)|is.null(rows)) stop(" rows invalid")
  if (is.null(columns)) columns=rep(1,length(rows))
  if (anyNA(columns)|any(is.nan(columns))|any(!is.finite(columns))|any(columns%%1!=0)|any(columns<1)) stop(" columns parameter invalid")
  if (length(columns)!=length(rows)) stop("rows and columns vectors must be the same length ")
  if (max(rows*columns)==1) { rows=1; columns=1} else {index=rows*columns>1; rows=rows[index]; columns=columns[index]}
  if ( replicates==2 && length(rows)>1 && length(columns)>1 && any((rows==2 & columns==2)[-length(rows)]) && is.null(model) ) 
      stop(" 2 x 2 row-and-column designs with two treatment replicates confound one treatment contrast between the row-by-column blocks
        and are unsuitable for complete factorial designs: you could try replacing the 2 x 2 row-and-column design by 4 x 1 row-and-column design.") 
  if (is.null(model)) model=paste0("~ ",paste0( unlist(lapply( 1:ncol(treatments), 
               function(i) {if (!is.factor(treatments[,i])) paste0("poly(",colnames(treatments)[i],",",length(unique(treatments[,i]))-1,")") else 
                            colnames(treatments)[i]})), collapse="*"))
  # constructs treatment  data frame possibly a fractional factorial design 
  Z=factorial(treatments,model,replicates)
  treatments=Z$TF
  fractionalEff=data.frame(Z$fraction,Z$eff)
  nunits=nrow(treatments)
  colnames(fractionalEff)=c("Fraction","D-Efficiency")
  fnames=colnames(treatments)
  strata=length(rows)
  blocksizes=nunits
  for (i in 1:strata)
    blocksizes=Sizes(blocksizes,i)
  regBlocks=isTRUE(all.equal(max(blocksizes), min(blocksizes)))
  
  if (is.null(searches)) 
    if (nunits<1000) searches=10000%/%nunits else if (nunits<5000) searches=5000%/%nunits else searches=1
  if( !is.finite(searches) | is.nan(searches) | searches<1 | searches%%1!=0 ) stop(" searches parameter is invalid")
  regReps=isTRUE(all.equal(max(replicates), min(replicates)))
  # tests for viable design sizes
  nestblocks=c(1,rows*columns)
  cumrows=cumprod(rows)
  cumcols=cumprod(columns)
  cumblocks=c(1,cumprod(rows*columns))
  
  if (cumrows[strata]*2>nunits) stop("Too many row blocks for the available plots  - every row block must contain at least two plots")
  if (cumcols[strata]*2>nunits) stop("Too many column blocks for the available plots  - every column block must contain at least two plots")
  if (cumblocks[strata+1]>nunits & cumrows[strata]>1 & cumcols[strata]>1) stop("Too many blocks - every row-by-column intersection must contain at least one plot")
  
 # nested factor level data frames
  df1=dataframes(rows,columns)
  nblkDesign=df1$nblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  nrowDesign=df1$nrowdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  ncolDesign=df1$ncoldesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE] 
  fblkDesign=df1$fblkdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  frowDesign=df1$frowdesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE]
  fcolDesign=df1$fcoldesign[rep(1:length(blocksizes),blocksizes),,drop=FALSE] 
  attempts=0
  CRB= (max(columns)==1 & ncol(treatments)==1 & is.factor(treatments[,1]) & length(rows)==1 & hcf%%rows[1]==0)
  TF=NULL
  while (is.null(TF) & attempts<10) {
    attempts=attempts+1
    TF=treatments
    colnames(TF)=fnames
    if (!CRB) {
      for ( i in 1:strata) {
        if (hcf%%cumrows[i]!=0) TF=blocksOpt(TF,fblkDesign[,i],frowDesign[,i],nrowDesign[,i],fblkDesign[,i])
        if (columns[i]>1)       TF=blocksOpt(TF,fblkDesign[,i],fcolDesign[,i],ncolDesign[,i],frowDesign[,i])
      }
    }
  }
  if (is.null(TF)) stop("Unable to find a non-singular solution for this design - please try a simpler block or treatment design")
  if (max(columns)==1) {
    fblkDesign=cbind(fblkDesign,Plots=factor(1:nunits))
    fblkDesign=data.frame(lapply(1:ncol(fblkDesign), function(r){sample(nlevels(fblkDesign[,r]))[fblkDesign[,r]]})) # Randomize labels - NB gives numeric columns
    fblkDesign=cbind(fblkDesign,TF)
    fblkDesign=fblkDesign[do.call(order, fblkDesign), ] # re-order
    blocksizes=table(fblkDesign[,ncol(fblkDesign)-ncol(TF)-1])[unique(fblkDesign[,ncol(fblkDesign)-ncol(TF)-1])]
    TF=fblkDesign[,c((ncol(fblkDesign)-ncol(TF)+1):ncol(fblkDesign)),drop=FALSE]
    for (i in 1 : ncol(treatments))
      if (is.factor(treatments[,i])) TF[,i]=as.factor(TF[,i])
    Design  = data.frame(df1$fblkdesign[rep(1:length(blocksizes),blocksizes ),-1,drop=FALSE],Plots=factor(1:nunits), TF)
    Efficiencies=FactBlocksEffics(Design)
    colnames(Design)=c(colnames(df1$fblkdesign)[-1],"Plots",fnames)
  }
  if (max(columns)>1) {
    rdf = data.frame(df1$frowdesign,df1$fcoldesign,df1$fblkdesign[,-1,drop=FALSE])[c(t(matrix(1:(3*strata),nrow=strata)))] # reorder columns
    rcDesign=data.frame(rdf[rep(1:length(blocksizes), blocksizes),],Plots=factor(1:nunits))
    rcDesign=data.frame(lapply(1:ncol(rcDesign), function(r){ sample(nlevels(rcDesign[,r]))[rcDesign[,r]]})  ) # Randomize
    rcDesign=data.frame(rcDesign, TF)
    rcDesign=rcDesign[do.call(order,rcDesign), ] # re-order
    blocksizes=table(rcDesign[,ncol(rcDesign)-ncol(TF)-1])[unique(rcDesign[,ncol(rcDesign)-ncol(TF)-1])]
    TF=rcDesign[,c((ncol(rcDesign)-ncol(TF)+1):ncol(rcDesign)),drop=FALSE]
    Efficiencies=FactRowColEffics( data.frame( rdf[rep(1:length(blocksizes), blocksizes),],Plots=factor(1:nunits),TF)) 
    ndf = data.frame(df1$nrowdesign,df1$ncoldesign)[c(t(matrix(1:(2*strata),nrow=strata)))]
    Design  = data.frame( ndf[rep(1:length(blocksizes), blocksizes),],Plots=factor(1:nunits),TF)  # rebuild factor levels
    colnames(Design)=c(colnames(ndf),"Plots",fnames)
  }
  # treatment replications
  TF=data.frame(TF[do.call(order, TF), ])  # re-order
  index=which(duplicated(TF)==FALSE)
  augmented_index=c(index, (nrow(TF)+1))
  TreatmentsTable=data.frame(TF[index,], Replication=diff(augmented_index) )
  colnames(TreatmentsTable)=c(fnames,"Replication")
  row.names(Design)=NULL
  row.names(Efficiencies)=NULL
  row.names(TreatmentsTable)=NULL
  list(Treatments=TreatmentsTable,model=model,DesignEfficiency=fractionalEff,BlocksEfficiency=Efficiencies,Design=Design,seed=seed,searches=searches,jumps=jumps)
}
