#' @title Optimum treatment set from a candidate set of treatments
#' 
#' @description
#' Finds an optimum set of treatments from a candidate set of treatments for any arbitrary 
#' treatments design formula.
#' 
#' @seealso \code{\link{design}}
#'  
#' @param treatments is a data frame or a list containing a candidate set of factorial treatments.  
#' 
#' @param size is the required number of treatments in the fractional set of treatments.
#' 
#' @param treatments_model is a model formula for the required treatments design.
#' 
#' @param restriction_model is a model formula which is a subset of the \code{treatments_model}
#' formula and which fixes those treatment factors contained in the restriction model formula.
#' 
#' @param searches are the maximum number of searches for selecting the best optimization.  
#' 
#' @details The candidate set \code{treatments} will normally contain one or more complete sets of 
#' treatment replicates. The algorithm re-arranges the rows of the \code{treatments} set to ensure 
#' that the first \code{size} rows of the optimized \code{treatments} set contains the optimized treatment
#' fraction. The maximum replication of any treatment is the number of times that treatment occurs in the 
#' candidate treatment set and for a polynomial response surface design extra replication of the candidate
#' set may be necessary to allow for differential replication of the design points. The design is
#' optimized with respect to the \code{treatments_model} conditional on the treatment factors
#' in the \code{restriction_model} being held constant. The \code{restriction_model} must be a subset
#' of the full \code{treatments_model} otherwise the design will be fully fixed and no further optimization
#' will be possible. Fitting a non-null \code{restriction_model} allows sequential optimization
#' with each successively \code{treatments_model} optimized conditional on all previously optimized models.
#' The D-optimal efficiency of the design for the optimized treatment set is calculated relative to the 
#' D-optimal efficiency of the design for the candidate treatment set. 
#' 
#' The default \code{treatments_model} parameter is an additive model for all treatment factors.
#' 
#' @return
#' A list containing:
#'  \itemize{
#'  \item{"TF"} {An optimized treatment fraction of the required \code{size}}
#'  \item{"fullTF"} {The full candidate set of treatments but with the first \code{size} 
#'  rows containing the optimized treatment fraction}
#'  \item{"Efficiency"} {The D-efficiency of the optimized treatment fraction relative to the full candidate set of treatments}
#' }
#'  
#' @examples
#' #' ## Plackett and Burman (P&B) type design for eleven 2-level factors in 12 runs 
#' ## NB. The algorithmic method is unlikely to succeed for larger P&B type designs. 
#'
#' GF = list(F1 = factor(1:2,labels=c("a","b")), F2 = factor(1:2,labels=c("a","b")), 
#'                  F3 = factor(1:2,labels=c("a","b")), F4 = factor(1:2,labels=c("a","b")),
#'                  F5 = factor(1:2,labels=c("a","b")), F6 = factor(1:2,labels=c("a","b")),
#'                  F7 = factor(1:2,labels=c("a","b")), F8 = factor(1:2,labels=c("a","b")), 
#'                  F9 = factor(1:2,labels=c("a","b")), F10= factor(1:2,labels=c("a","b")), 
#'                  F11= factor(1:2,labels=c("a","b")) )
#' model = ~ F1 + F2 + F3 + F4 + F5 + F6 + F7 + F8 + F9 + F10 + F11
#' Z=fraction(GF,size=12,treatments_model=model,searches=100)
#' print(Z$TF)
#' print(Z$Efficiency)
#' round(crossprod(scale(data.matrix(Z$TF))),6)
#' 
#' ## Factorial treatment designs defined by sequentially fitted factorial treatment models
#' ## 4 varieties by 3 levels of N by 3 levels of K assuming degree-2 treatment model in 24 plots.
#' ## The single stage model gives an unequal split for the replication of the four varieties
#' ## whereas the two stage model forces an equal split of 6 plots per variety.
#' ## The single stage model is slightly more efficient overall (about 1.052045 versus 1.043662)
#' ## but unequal variety replication is undesirable if all varieties are equally important.
#' 
#' ## model terms
#' treatments = list(Variety = factor(1:4), N = 1:3, K = 1:3)
#' variety_model = ~ Variety
#' full_model = ~ (Variety + N + K)^2  + I(N^2) + I(K^2)
#' 
#' ## single stage model
#' opt_full_treatments = fraction(treatments,24,full_model,searches=10)
#' opt_full_treatments$Efficiency
#' table(opt_full_treatments$TF[,1]) # variety replication
#' 
#' ## two stage model
#' opt_var_treatments  = fraction(treatments,24,variety_model,searches=10)
#' opt_full_treatments = fraction(opt_var_treatments$fullTF,24,full_model,variety_model,searches=10)
#' opt_full_treatments$Efficiency
#' table(opt_full_treatments$TF[,1]) # variety replication
#' 
#' @export
 fraction = function(treatments,size,treatments_model=NULL,restriction_model=NULL,searches=50) {
  tol = .Machine$double.eps ^ 0.5
  options(contrasts = c('contr.treatment', 'contr.poly'),
          warn = 0)
  
  # *************************************************************************************************************
  # Updates variance matrix for swapped rows where mi is swapped out and qj is swapped in
  # *************************************************************************************************************
  fractUpDate = function(V, mi, qj) {
    f = crossprod(V, qj)
    g = crossprod(V, mi)
    a = as.numeric(crossprod(qj, f))
    b = as.numeric(crossprod(mi, g))
    c = as.numeric(crossprod(mi, f))
    W = g * sqrt(1 + a) - f * c / sqrt(1 + a)
    U = f * sqrt(1 - b + (c * c) / (1 + a))
    V = V - (tcrossprod(U) - tcrossprod(W)) / ((1 + a) * (1 - b) + c * c)
    return(V)
  }
  # *************************************************************************************************************
  # non singular starting design for the first nunits of TF given the design matrix TM
  # *************************************************************************************************************
  nonSingular = function(TF, TM, fixed, size, rank)  {
    
    for (z in seq_len(10000)) {
      k = sample(nlevels(fixed), 1)
      available = which(fixed == levels(fixed)[k])
      swapout = available[available <= size]
      swapin = available[available > size]
      if (length(swapin) == 0 | length(swapout) == 0)
        next
      if (length(swapin) > 1)
        p2 = sample(swapin, 1)
      else
        p2 = swapin
      if (length(swapout) > 1)
        p1 = sample(swapout, 1)
      else
        p1 = swapout
      TM[c(p1, p2), ] = TM[c(p2, p1), ]
      TF[c(p1, p2), ] = TF[c(p2, p1), ]
      newrank = qr(TM[1:size, ])$rank
      if (newrank < rank) {
        TM[c(p1, p2), ] = TM[c(p2, p1), ]
        TF[c(p1, p2), ] = TF[c(p2, p1), ]
      } else
        rank = newrank
      if (ncol(TM) == rank)
        break
    }
    if (rank < ncol(TM))
      stop ("Cannot find a starting design")
    return(list(TF = TF, TM = TM))
  }
  # *******************************************************************************************************************
  #  treatment information
  # *******************************************************************************************************************
  information = function(TF,size, treatments_model) {
    TM = model.matrix(treatments_model, TF)
    fullinf = exp(determinant(crossprod(TM), logarithm = TRUE)$modulus / ncol(TM)) / nrow(TM)
    TM=TM[1:size,]
    modinf = exp(determinant(crossprod(TM), logarithm = TRUE)$modulus / ncol(TM)) / nrow(TM)
    rel_info=modinf/fullinf
    return(rel_info)
  }
  # *************************************************************************************************************
  # main program
  # *************************************************************************************************************
  TF=treatments
  if (is.list(TF) & !is.data.frame(TF)) 
    TF=expand.grid(TF) else
      TF = data.frame(TF)
  if (is.null(treatments_model))
    treatments_model = as.formula( paste("~",paste0(colnames(TF), collapse = "+")))
  TFlive=colnames(TF) %in% all.vars(treatments_model)

  if (!is.null(restriction_model)) {
    tnames = names(treatments) %in%  all.vars(restriction_model)
    fixed = interaction(data.frame(lapply(TF[, tnames, drop = FALSE], factor)), drop = TRUE)
  } else
    fixed = factor(rep(1, nrow(TF)))
  TM = model.matrix(treatments_model, TF)
  if (ncol(TM) > size)
    stop("too many treatment parameters for the given number of experimental units")
  # information matrix for full treatment set before fractionation
  Dmax = exp((determinant(crossprod(TM), logarithm = TRUE)$modulus) /
               ncol(TM)) / nrow(TM)
  gDfrac = 0
  for (i in 1:searches) {
    # randomize treatments within the fixed levels of TF
    for (j in 1:nlevels(fixed))
      TF[seq(nrow(TF))[fixed == levels(fixed)[j]],] = TF[sample(seq(nrow(TF))[fixed ==
                                                                                levels(fixed)[j]]), , drop = FALSE]
    TM = model.matrix(as.formula(treatments_model), TF)
    rank = qr(TM[1:size,])$rank
    if (rank < ncol(TM)) {
      nonSing = nonSingular(TF, TM, fixed, size, rank)
      TF = nonSing$TF
      TM = nonSing$TM
    }
    V = chol2inv(chol(crossprod(TM[1:size, , drop = FALSE]))) # V is the variance matrix of the starting design
    counter = 0
    locrelD = 1
    repeat {
      kmax = 1
      counter = counter + 1
      for (k in 1:nlevels(fixed)) {
        available = which(fixed == levels(fixed)[k])
        swapout = available[available <= size]
        swapin = available[available > size]
        if ((length(swapin) == 0 |
             length(swapout) == 0) & k < nlevels(fixed))
          next
        if ((length(swapin) == 0 |
             length(swapout) == 0) & k == nlevels(fixed))
          break
        if (length(swapin) > 1)
          swap_in = sample(swapin, min(1024, length(swapin)))
        else
          swap_in = swapin
        if (length(swapout) > 1)
          swap_out = sample(swapout, min(1024, length(swapout)))
        else
          swap_out = swapout
        M1VM2 =        tcrossprod(tcrossprod(TM[swap_out, , drop =
                                                  FALSE], V), TM[swap_in, , drop = FALSE])
        M1VM1 = 1 - diag(tcrossprod(tcrossprod(TM[swap_out, , drop =
                                                    FALSE], V), TM[swap_out, , drop = FALSE]))
        M2VM2 = 1 + diag(tcrossprod(tcrossprod(TM[swap_in, , drop =
                                                    FALSE], V), TM[swap_in, , drop = FALSE]))
        Z = M1VM2 ** 2 + tcrossprod(M1VM1, M2VM2)
        z = which(Z == max(Z), arr.ind = TRUE)[1,]
        if (Z[z[1], z[2]] > kmax) {
          kmax = Z[z[1], z[2]]
          pi = swap_out[z[1]]
          pj = swap_in[z[2]]
        }
      }
      if (kmax > (1 + tol)) {
        counter = 0
        locrelD = locrelD * kmax
        V = fractUpDate(V, TM[pi,], TM[pj,])   # parameters(V,row_swappedout,row_swappedin)
        TM[c(pi, pj),] = TM[c(pj, pi),]
        TF[c(pi, pj),] = TF[c(pj, pi),]
      }
      if (counter == 5)
        break
      if (counter == 1 &
          length(swapin) <= 1024 & length(swapout) <= 1024)
        break
    }
    Dfrac = exp((determinant(crossprod(TM[1:size, , drop = FALSE]), logarithm = TRUE)$modulus) /
                  ncol(TM)) / size
    if (Dfrac > (gDfrac + tol)) {
      gTF = TF
      gDfrac = Dfrac
    }
    if ( all(sapply(TF[,TFlive,drop=FALSE], is.factor)) &
        isTRUE(all.equal(gDfrac, Dmax, tolerance = tol)))
      break
  
  }
  Effics = information(gTF,size,treatments_model)
 return(list(TF=gTF[1:size,,drop=FALSE],fullTF=gTF,Efficiency=Effics))
}
