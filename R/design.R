##' @title General block and treatment designs.
#'
#' @description
#' Constructs D-optimal block and treatment designs for any feasible linear treatment model and
#' any feasible combination of block factors.
#' 
#' @param treatments a candidate set of treatments for any combination of qualitative or quantitative level factors.
#' 
#' @param blocks a data frame of block factors.
#' 
#' @param treatments_model a list containing one or more nested treatment model formula. 
#' 
#' @param weighting a weighting factor between 0 and 1 for weighting the 2-factor interaction effects of
#' factorial blocks.
#'
#' @param seed an integer initializing the random number generator.
#'
#' @param searches the maximum number of local optima searched at each stage of an
#' optimization.
#'
#' @param jumps the number of pairwise random treatment swaps used to escape a local maxima. 
#' 
#' @details
#' 
#' \code{treatments} a factor or data frame or list of generators for the candidate set of treatments. 
#' The candidate treatment set can include both quantitative and qualitative level factors and the required design 
#' is selected from the candidate treatment set without replacement. The number of times 
#' a treatment can be included in a design depends on the replication of that treatment in the candidate set therefore
#' increasing the treatment replication in the candidate set allows increased replication for individual treatments. 
#' 
#' \code{blocks} a data frame of nested or crossed block factors. The length of the block factors defines the 
#' required size of the treatment design. the treatment design is optimized for the required blocks design by maximising
#' the determinant of the block-adjusted treatment information matrix (D-optimality. For crossed block designs,
#' the information matrix is based on a weighted combination of block main effects and two-factor interaction effects
#' where the information on the interaction effects is down-weighted by a coefficient 0 < w < 1 where w = 0 
#' gives blocks main effects only whereas w = 1 gives the full two-factor blocks interaction model.
#' This process ensures that factorial main effects can be given more importance than the interaction effects
#'  in the optimization of a crossed blocks design. 
#' 
#' \code{treatments_model} a character vector containing one or more treatments model formula. The 
#' treatment models are optimized sequentially for each model formula in turn assuming the treatments of
#' any previously fitted models are held constant. Sequential fitting provides improved flexibility 
#' for fitting treatment factors or variables of different status or importance (see examples below). 
#' The design criterion for each treatment model is maximization of the determinant of the information 
#' matrix of that treatment model conditional on the the information matrices of all previously 
#' fitted treatment models.
#' 
#' The D-optimal treatment design efficiency is the ratio of the generalized variance of the full 
#' candidate treatment model relative to the generalized variance of the optimized design. The efficiency
#' will be less than or equal to 1 for factorial models but may exceed 1 for polynomial models. 
#' 
#' The design outputs include a table of efficiency factors for first and second-order factorial block 
#' effects and comparison of the efficiency factors for different choices of w can be used to find a 
#' good compromise design for fitting both main block and interaction block effects.
#' 
#' For more details see \code{vignette(package = "blocksdesign")}  
#' 
#' @return
#' \item{Replication}{The treatments included in the design and the replication of each individual 
#'  treatment taken in de-randomized standard order.}
#' \item{Design}{The design layout showing the randomized allocation of treatments to blocks and plots
#' in standardized block order.}
#' \item{Treatments_model}{The fitted treatment model, the number of model parameters (DF)
#'   and the D-efficiency of each sequentially fitted treatment model.}
#' \item{Blocks_model}{The blocks sub-model design and the D- and A-efficiency factors of each successively
#'  fitted sub-blocks model.}
#' \item{seed}{Numerical seed for random number generator.}
#' \item{searches}{Maximum number of searches in each stratum.}
#' \item{jumps}{Number of random treatment swaps to escape a local maxima.}
#' 
#' @references
#' 
#' Cochran W. G. & Cox G. M. (1957) Experimental Designs 2nd Edition John Wiley & Sons
#' 
#' Gupta, S. C., and B. Jones. Equireplicate Balanced Block Designs with Unequal Block Sizes.
#' Biometrika, vol. 70, no. 2, 1983, pp. 433–440. 
#' 
#' @examples
#' ## For best results, the number of searches may need to be increased.
#' 
#' ## 4 replicates of 12 treatments with 16 nested blocks of size 3
#' ## giving a rectangular lattice: see Plan 10.10 Cochran and Cox 1957.
#' 
#' blocks = data.frame(Main = gl(4,12), Sub = gl(16,3))
#' treatments = data.frame(treatments = gl(12,1,48))
#' Z=design(treatments, blocks)
#' print(Z)
#' 
#' ## 3 replicates of 15 treatments in 3 main blocks with two sets of
#' ## nested blocks and one control treatment
#' \donttest{
#' blocks=data.frame( Main = gl(3,18,54),Sub1 = gl(9,6,54),Sub2 = gl(27,2,54))
#' treatments=factor(rep(c(1:15,rep("control",3)),3))
#' Z=design(treatments,blocks)
#' print(Z)
#' incid1=table(interaction(Z$Design$Main,Z$Design$Sub1,lex.order = TRUE),Z$Design$treatments)
#' crossprod(incid1) # print pairwise concurrences within Sub1 blocks
#' incid2=table(interaction(Z$Design$Main,Z$Design$Sub2,lex.order = TRUE),Z$Design$treatments)
#' crossprod(incid2) # print pairwise concurrences within Sub2 blocks
#' }
#'  
#' ## 4 x 12 design for 4 replicates of 12 treatments with 3 plots in each intersection block
#' ## and  3 sub-columns nested within each main column. The optimal row-and column design
#' ## is Trojan with A-efficiency = 22/31 for the intersection blocks
#'  \donttest{
#' blocks = data.frame(Rows = gl(4,12), Cols = gl(4,3,48), subCols = gl(12,1,48))
#' treatments = data.frame(treatments =gl(12,1,48))
#' Z=design(treatments,blocks)
#' print(Z)
#' incid=table(interaction(Z$Design$Rows,Z$Design$Cols,lex.order = TRUE),Z$Design$treatments)
#' crossprod(incid) # print pairwise concurrences within blocks
#' }
#' 
#' ## 4 x 13 Row-and-column design for 4 replicates of 13 treatments 
#' ## Youden design Plan 13.5 Cochran and Cox (1957).
#' \donttest{
#' blocks = data.frame(Rows = gl(4,13), Cols = gl(13,1,52))
#' treatments = data.frame(treatments= gl(13,1,52))
#' Z=design(treatments,blocks,searches = 700)
#' print(Z)
#' incid=table(Z$Design$Cols,Z$Design$treatments)
#' crossprod(incid) # print pairwise concurrences of treatments within column blocks (BIB's)
#' }
#' 
#' ## 48 treatments in 2 replicate blocks with 2 nested rows in each replicate and 3 main columns
#' 
#' blocks = data.frame(Reps = gl(2,48), Rows = gl(4,24,96), Cols = gl(3,8,96))
#' treatments = data.frame(treatments=gl(48,1,96))
#' Z=design(treatments,blocks,searches=5)
#' print(Z)
#' 
#' ## 48 treatments in 2 replicate blocks with 2 main columns
#' ## The default weighting gives non-estimable Reps:Cols effects due to inherent aliasing
#' ## Increased weighting gives estimable Reps:Cols effects but non-orthogonal main effects
#' \donttest{
#' blocks = data.frame(Rows = gl(2,48), Cols = gl(2,24,96))
#' Z1=design(treatments=gl(48,1,96),blocks,searches=5)
#' print(Z1)
#' Z2=design(treatments=gl(48,1,96),blocks,searches=5,weighting=.9)
#' print(Z2)
#' }
#' 
#' ## Equireplicate balanced designs with unequal block sizes. Gupta & Jones (1983). 
#' ## Check for equality of D and A-efficiency in optimized design
#' 
#' t=factor(c(rep(1:12,each=7)))
#' b=factor(c(rep(1:12,each=6),rep(13:18,each=2)))
#' Z=design(t,b,searches=100)$Blocks_model # max efficiency = 6/7
#' print(Z)
#' 
#' \donttest{
#' t=factor(c(rep(1:14,each=8)))
#' b=factor(c(rep(1:14,each=4),rep(15:21,each=8)))
#' Z=design(t,b,searches=500)$Blocks_model # max efficiency = 7/8
#' print(Z)
#' }
#' \donttest{
#' t=factor(c(rep(1:16,each=7)))
#' b=factor(c(rep(1:16,each=4),rep(17:22,each=8)))
#' Z=design(t,b,searches=1000)$Blocks_model # max efficiency = 6/7
#' print(Z)
#' }
#' \donttest{
#' t=factor(c(rep(1:18,each=7)))
#' b=factor(c(rep(1:18,each=6),rep(19:24,each=3)))
#' Z=design(t,b,searches=500)$Blocks_model # max efficiency = 6/7
#' print(Z)
#' }
#' 
#' ## Factorial treatment designs defined by a factorial treatment model
#' 
#' ## Main effects of five 2-level factors in a half-fraction in 2/2/2 nested blocks design 
#' ## The algorithmic search method is likely to require a very long search time to reach the
#' ## known orthogonal solution for this regular fractional factorial design 
#' \donttest{
#' treatments = list(F1 = gl(2,1), F2 = gl(2,1),F3 = gl(2,1), F4 = gl(2,1), F5 = gl(2,1))
#' blocks = data.frame(b1 = gl(2,8),b2 = gl(4,4),b3 = gl(8,2))
#' model= ~ F1 + F2 + F3 + F4 + F5
#' repeat {Z = design(treatments,blocks,treatments_model=model,searches=50)
#' if ( isTRUE(all.equal(Z$Blocks_model[3,3],1) ) ) break }
#' print(Z)
#' }
#'  
#' # Second-order model for five qualitative 2-level factors in 4 randomized blocks with two
#' # nested sub-blocks in each main plot
#' treatments = list(F1 = gl(2,1), F2 = gl(2,1),F3 = gl(2,1), F4 = gl(2,1), F5 = gl(2,1))
#' blocks = data.frame(main = gl(4,8),sub=gl(8,4))
#' model = ~ (F1 + F2 + F3 + F4 + F5)^2
#' Z=design(treatments,blocks,treatments_model=model,searches = 10)
#' print(Z)
#' 
#' # Main effects of five 2-level factors in a half-fraction of a 4 x 4 row-and column design
#' \donttest{
#' treatments = list(F1 = gl(2,1), F2 = gl(2,1),F3 = gl(2,1), F4 = gl(2,1), F5 = gl(2,1))
#' blocks = data.frame(rows = gl(4,4), cols = gl(4,1,16))
#' model = ~ F1 + F2 + F3 + F4 + F5
#' repeat {Z = design(treatments,blocks,treatments_model=model,searches=50)
#' if ( isTRUE(all.equal(Z$Blocks_model[2,3],1) ) ) break}
#' print(Z)
#' }
#' 
#' # Quadratic regression for three 3-level numeric factor assuming a 10/27 fraction
#' treatments = list(A = 1:3, B = 1:3, C = 1:3)
#' blocks=factor(rep(1,10))
#' model = ~ ( A + B + C)^2 + I(A^2) + I(B^2) + I(C^2)
#' Z=design(treatments,blocks,treatments_model=model,searches=10) 
#' print(Z)
#' 
#' ## response surface designs with doubled treatment candidate sets 
#' ## allowing duplicated design points
#' 
#' ## two treatment factor design showing double replication on the 
#' ## four corner points and double replication on the centre point
#' treatments = list( V1 = 1:3, V2 = rep(1:3,2))
#' blocks=factor(rep(1,14))
#' model =  ~ (V1 + V2)^2 + I(V1^2) + I(V2^2)
#' design(treatments,blocks,treatments_model=model,searches=10)
#' 
#' ## three treatment factor design with eight corner points, 12 edge points and 1 centre 
#' ## point. All 21 design points are accounted for leaving nothing for extra replication.
#' treatments = list( V1 = 1:3, V2 = 1:3, V3 = rep(1:3,2))
#' blocks=factor(rep(1,21))
#' model =  ~ (V1 + V2 + V3)^2 + I(V1^2) + I(V2^2) + I(V3^2)
#' design(treatments,blocks,treatments_model=model,searches=20)
#' 
#' ## as above but with extra replication on each corner point and on the centre point
#' treatments = list( V1 = 1:3, V2 = 1:3, V3 = rep(1:3,2))
#' blocks=factor(rep(1,30))
#' model =  ~ (V1 + V2 + V3)^2 + I(V1^2) + I(V2^2) + I(V3^2)
#' design(treatments,blocks,treatments_model=model,searches=20)
#' 
#' ## Factorial treatment designs defined by sequentially fitted factorial treatment models
#' 
#' ## 4 varieties by 3 levels of N by 3 levels of K assuming degree-2 treatment
#' ## interaction effects and two blocks of 12 plots
#' ## the single stage model gives an unequal split for the replication of the four varieties
#' ## which may be undesirable whereas the two stage model forces an equal split of 6 plots
#' ## per variety. The single stage model is slightly more efficient than the two stage model
#' ## (1.052045 versus 0.9761135 standardized relative to the full factorial design)  but
#' ## in this example global D-optimality does not give the most suitable design structure. 
#' \donttest{
#' treatments = list(Variety = factor(1:4), N = 1:3, K = 1:3)
#' blocks = data.frame(main=gl(2,12))
#' treatments_model = ~ (Variety + N + K)^2  + I(N^2) + I(K^2)
#' Z1 = design(treatments,blocks,treatments_model=treatments_model,searches=10)
#' print(Z1)
#' treatments_model = list( ~ Variety , ~ Variety + (Variety + N + K)^2 + I(N^2) + I(K^2))
#' Z2 = design(treatments,blocks,treatments_model=treatments_model,searches=10)
#' print(Z2)
#' }
#' 
#' @export
#' @importFrom stats anova lm model.matrix as.formula setNames 
#' @importFrom plyr match_df
#' @importFrom plyr count
#' @importFrom plyr is.formula
#' 
design = function(treatments,
                  blocks,
                  treatments_model = NULL,
                  weighting = 0.5,
                  searches = NULL,
                  seed = NULL,
                  jumps = 1) {
  
  tol = .Machine$double.eps ^ 0.5
  options(contrasts = c('contr.treatment', 'contr.poly'),
          warn = 0)
  
  # ***********************************************************************************************************
  # dMat is a matrix where the ith, jth matrix element is the change in the D-max criterion of the block
  # design due to swapping ith and jth treatments of a sample s of the set of plots
  # ************************************************************************************************************
  dMat = function(TM, BM, V, s) {
    sTM = TM[s, , drop = FALSE]
    sBM = BM[s, , drop = FALSE]
    TT = tcrossprod(tcrossprod(sTM, V[1:ncol(TM), 1:ncol(TM), drop = FALSE]), sTM)
    BB = tcrossprod(tcrossprod(sBM, V[(ncol(TM) + 1):ncol(V), (ncol(TM) +
                                                                 1):ncol(V), drop = FALSE]), sBM)
    TB = tcrossprod(tcrossprod(sTM, V[(ncol(TM) + 1):ncol(V), 1:ncol(TM), drop =
                                        FALSE]), sBM)
    TT = sweep(sweep(2 * TT, 1, diag(TT)), 2, diag(TT))
    BB = sweep(sweep(2 * BB, 1, diag(BB)), 2, diag(BB))
    TBBT = sweep(sweep(TB + t(TB), 1, diag(TB)), 2, diag(TB))
    dMat = (1 + TBBT) ** 2 - TT * BB
    return(dMat)
  }
  # *************************************************************************************************************
  # Updates variance matrix for pairs of swapped plot treatments using standard matrix updating formula
  # mtb**2-mtt*mbb is > 0 because the swap is a positive element of dMat=(TB+t(TB)+1)**2-TT*BB
  # 2*mtb+mtt+mbb > mtt + mbb + 2*(mtt*mbb)**.5 > 0 because mtb**2 > mtt*mbb
  # general designs with factorial treatment structures or crossed blocks designs
  # *************************************************************************************************************
  UpDate = function(V, t, b, v) {
    VTTt = crossprod(V[1:v, 1:v, drop = FALSE], t)
    VBBb = crossprod(V[(1 + v):ncol(V), (1 + v):ncol(V), drop = FALSE], b)
    VTBt = crossprod(V[1:v, (1 + v):ncol(V), drop = FALSE], t)
    VTBb = crossprod(V[(1 + v):ncol(V), 1:v, drop = FALSE], b)
    tMt = as.numeric(crossprod(t, VTTt))
    bMb = as.numeric(crossprod(b, VBBb))
    tMb = as.numeric(crossprod(b, VTBt))
    f1 = (VTTt + VTBb) / sqrt(2)
    f2 = (VBBb + VTBt) / sqrt(2)
    g1 = (VTBb - VTTt) / sqrt(2)
    g2 = (VBBb - VTBt) / sqrt(2)
    a = (tMt + bMb + 2 * tMb) / 2
    b = (tMt + bMb - 2 * tMb) / 2
    c = (bMb - tMt) / 2
    div = (1 + a) * (1 - b) + c * c
    V[1:v, 1:v] = V[1:v, 1:v, drop = FALSE] - (tcrossprod(f1) * (1 - b) - tcrossprod(g1) *
                                                 (1 + a) +
                                                 (tcrossprod(g1, f1) + tcrossprod(f1, g1)) *
                                                 c) / div
    V[(1 + v):ncol(V), (1 + v):ncol(V)] = V[(1 + v):ncol(V), (1 + v):ncol(V), drop =
                                              FALSE] -
      (tcrossprod(f2) * (1 - b) - tcrossprod(g2) * (1 + a) + (tcrossprod(g2, f2) +
                                                                tcrossprod(f2, g2)) * c) / div
    V[1:v, (1 + v):ncol(V)] = V[1:v, (1 + v):ncol(V), drop = FALSE] -
      (
        tcrossprod(f1, f2) * (1 - b) - tcrossprod(g1, g2) * (1 + a) + (tcrossprod(g1, f2) +
                                                                         tcrossprod(f1, g2)) * c
      ) / div
    V[(1 + v):ncol(V), 1:v] = t(V[1:v, (1 + v):ncol(V), drop = FALSE])
    return(V)
  }
  # **************************************************************************************************************
  # Maximises the matrix dMat=TB**2-TT*BB to compare and choose the best swap for D-efficiency improvement.
  # Sampling is used initially but later a full search is used to ensure steepest ascent optimization.
  # **************************************************************************************************************
  DMax = function(TF, TM, BM, V, mainBlocks) {
    locrelD = 1
    mainSets = tabulate(mainBlocks)
    nSamp = pmin(rep(8, nlevels(mainBlocks)), mainSets)
    repeat {
      kmax = 1
      for (k in 1:nlevels(mainBlocks)) {
        s = sort(sample((1:nrow(TF))[mainBlocks == levels(mainBlocks)[k]], nSamp[k]))
        dMat = dMat(TM, BM, V, s)
        z = which(dMat == max(dMat, na.rm = TRUE), arr.ind = TRUE)[1, ]
        if (dMat[z[1], z[2]] > kmax) {
          kmax = dMat[z[1], z[2], drop = FALSE]
          pi = s[z[1]]
          pj = s[z[2]]
        }
      }
      if (kmax > (1 + tol)) {
        locrelD = locrelD * kmax
        V = UpDate(V, (TM[pi, , drop = TRUE] - TM[pj, , drop = TRUE]), (BM[pj, , drop =
                                                                             TRUE] - BM[pi, , drop = TRUE]), ncol(TM))
        TM[c(pi, pj), ] = TM[c(pj, pi), , drop = FALSE]
        TF[c(pi, pj), ] = TF[c(pj, pi), , drop = FALSE]
      }  else if (sum(nSamp) == nrow(TF))
        break
      else
        nSamp = pmin(mainSets, 2 * nSamp)
    }
    return(list(
      V = V,
      locrelD = locrelD,
      TF = TF,
      TM = TM
    ))
  }
  # *****************************************************************************************************************
  # Random swaps for blocks optimization
  # *****************************************************************************************************************
  Swaps = function(TF, BF, pivot, rank, restrict, blocks) {
    candidates = NULL
    counter = 0
    while (length(candidates) == 0 & counter < 500) {
      counter = counter + 1
      swapOut = pivot[1:rank]
      swapIn = pivot[(1 + rank):nrow(TF)]
      if (length(swapIn) > 1)
        s1 = sample(swapIn, 1)
      else
        s1 = swapIn
      rownames(TF) = NULL
      sameTrts = suppressMessages(as.integer(rownames(match_df(TF, TF[s1, , drop =
                                                                        FALSE])))) # row replicates of TF[s1,]
      candidates = seq_len(nrow(TF))[restrict == restrict[s1] &
                                       blocks != blocks[s1] &
                                       !(seq_len(nrow(TF)) %in% sameTrts) &
                                       (seq_len(nrow(TF)) %in% swapOut)]
    }
    if (length(candidates) == 0)
      stop(" Unable to find an initial non-singular starting design ")
    if (length(candidates) > 1)
      s2 = sample(candidates, 1)
    else
      s2 = candidates
    return(c(s1, s2))
  }
  # *******************************************************************************************************************
  # Initial randomized starting design. If the initial design is rank deficient,
  # random swaps with positive selection are used to to increase design rank
  # *******************************************************************************************************************
  blocksNonSingular = function(TF, TM, BF, BM, restrict, degree2) {
    fullrank = FALSE
    D=cbind(TM, BM)
    if (ncol(D)<=nrow(D)) {
      Q = qr(t(D))
      rank = Q$rank
      pivot = Q$pivot
      for (i in 1:500) {
        if (rank < ncol(D)) {
          swap = Swaps(TF, BF, pivot, rank, restrict, BF)
          TM[c(swap[1], swap[2]), ] = TM[c(swap[2], swap[1]), , drop =
                                           FALSE]
          Q = qr(t(cbind(TM, BM)))
          if (Q$rank > rank) {
            rank = Q$rank
            pivot = Q$pivot
            TF[c(swap[1], swap[2]), ] = TF[c(swap[2], swap[1]), , drop =
                                             FALSE]
          } else {
            TM[c(swap[1], swap[2]), ] = TM[c(swap[2], swap[1]), , drop = FALSE]
          }
        } else {
          fullrank = TRUE
          break
        }
      }
    }
    return(list(
      TF = TF,
      TM = TM,
      fullrank = fullrank,
      degree2 = degree2
    ))
  }
 
  # *****************************************************************************************************************
  # Blocks design efficiency factors
  # ******************************************************************************************************************
  effBlocks  = function(TF, TM, BF, deg1_BM, deg2_BM) {
    D_Effic  = rep(NA, length = ncol(BF))
    A_Effic  = rep(NA, length = ncol(BF))
    D_IntEff = rep(NA, length = ncol(BF))
    A_IntEff = rep(NA, length = ncol(BF))
    deg1_levs = sapply(1:ncol(BF), function(u)
      ncol(deg1_BM[[u]]))
    deg2_levs = sapply(1:ncol(BF), function(u)
      ncol(deg2_BM[[u]]))
    for (u in 1:ncol(BF)) {
      E1 = eigen(
        diag(ncol(TM)) - tcrossprod(crossprod(TM, deg1_BM[[u]])),
        symmetric = TRUE,
        only.values = TRUE
      )
      D_Effic[u] = 0
      A_Effic[u] = 0
      if (all(E1$values > tol)) {
        D_Effic[u] = exp(sum(log(E1$values)) / ncol(TM))
        A_Effic[u] = ncol(TM) / sum(1 / E1$values)
      }
      if (ncol(deg2_BM[[u]]) > ncol(deg1_BM[[u]]) &
          weighting > 0) {
        E2 = eigen(
          diag(ncol(TM)) - tcrossprod(crossprod(TM, deg2_BM[[u]])),
          symmetric = TRUE,
          only.values = TRUE
        )
        D_IntEff[u] = 0
        A_IntEff[u] = 0
        if (all(E2$values > tol)) {
          D_IntEff[u] = exp(sum(log(E2$values)) / ncol(TM))
          A_IntEff[u] = ncol(TM) / sum(1 / E2$values)
        }
      }
    }
    if (isTRUE(all.equal(deg1_levs, deg2_levs))) {
      Effics = data.frame(
        Level = 1:ncol(BF),
        Blocks = sapply(BF[, 1:ncol(BF), drop = FALSE], nlevels),
        "D-Efficiency" = round(D_Effic, 6),
        "A-Efficiency" = round(A_Effic, 6),
        Bound = 1
      )
    } else {
      Effics1 = data.frame(
        First_order = format(first_order),
        effects = deg1_levs,
        "D-Effic" = round(D_Effic, 6),
        "A-Effic" = round(A_Effic, 6)
      )
      Effics2 = data.frame(
        Second_order = format(second_order),
        effects = deg2_levs,
        "D-Effic" = round(D_IntEff, 6),
        "A-Effic" = round(A_IntEff, 6)
      )
      Effics = cbind(Effics1, Effics2)
    }
    return(Effics)
  }
  
  # *****************************************************************************************************************
  # Orthogonal basis for M
  # ******************************************************************************************************************
  orthogM = function(M) {
    QR = qr(M) # qr transformation
    qr.Q(QR)[, 1:QR$rank, drop = FALSE] # orthogonal basis for M where Q'Q=I
  }
  
  # *****************************************************************************************************************
  # TF and BF are factor levels and TM and BM are indicator matrices for each factor level. Centering eliminates the
  # the block and treatment means. Optimizes the  blocks assuming a possible set of Main block constraints.
  # If the initial design is rank deficient, random swaps with positive selection are used to to increase design rank
  # Each factor in BF may include levels from preceding factors hence pivoting sweeps out factor dependencies
  # ******************************************************************************************************************
  optBlocks = function(TF, BF) {
    TM = orthogM(model.matrix(treatments_model[[length(treatments_model)]], TF))[,-1, drop = FALSE] #omit mean
    if (nrow(TM) < (ncol(TM) + nlevels(BF[, ncol(BF)])))
      stop("Too many model parameters for the size of the design")
    Effic=NULL
    
    if (nlevels(BF[,1])>1) {
    deg1_BM = lapply(1:ncol(BF),function(u) orthogM(model.matrix(first_order[[u]],BF))[,-1,drop = FALSE])
    deg2_BM = lapply(1:ncol(BF),function(u) orthogM(model.matrix(second_order[[u]],BF))[,-1,drop = FALSE])

    for (u in 1:ncol(BF)) {
      if (u>1) Restrict = droplevels(interaction(BF[,1:(u-1)])) else
        Restrict = factor(rep(1, nrow(BF)))
      NS = blocksNonSingular(TF, TM, BF[, u], deg2_BM[[u]], Restrict, TRUE)
      if (!NS$fullrank)
        NS = blocksNonSingular(TF, TM, BF[, u], deg1_BM[[u]], Restrict, FALSE)
      if (!NS$fullrank)
        stop("Cannot find a non-singular starting block design ") 
      TM = NS$TM
      TF = NS$TF
      if (NS$degree2 & weighting > 0)
        BM = deg2_BM[[u]]
      else
        BM = deg1_BM[[u]]
      Info = crossprod(cbind(TM, BM))
      # down weights the two factor information, if present
      if (ncol(deg2_BM[[u]]) > ncol(deg1_BM[[u]]) &
          NS$degree2 & weighting > 0)
        Info[(ncol(TM) + ncol(deg1_BM[[u]]) + 1):ncol(Info), (ncol(TM) +
                                                                ncol(deg1_BM[[u]]) + 1):ncol(Info)] =
        Info[(ncol(TM) + ncol(deg1_BM[[u]]) + 1):ncol(Info), (ncol(TM) +
                                                                ncol(deg1_BM[[u]]) + 1):ncol(Info), drop = FALSE] / (weighting ^ 2)
      globrelD = 0
      relD = 1
      globTM = TM
      globTF = TF
      V = chol2inv(chol(Info))

      for (r in 1:searches) {
        dmax = DMax(TF, TM, BM, V, Restrict)
        if (dmax$locrelD > (1 + tol)) {
          relD = relD * dmax$locrelD
          TM = dmax$TM
          TF = dmax$TF
          V = dmax$V
          if (relD > globrelD) {
            globTM = TM
            globTF = TF
            globrelD = relD
          }
        }
        if (r == searches)
          break
        for (iswap in 1:jumps) {
          counter = 0
          repeat {
            counter = counter + 1
            s1 = sample(seq_len(nrow(TF)), 1)
            rownames(TF) = NULL
            unavailable = suppressMessages(as.integer(rownames(match_df(TF, TF[s1, , drop =
                                                                                 FALSE]))))
            z = seq_len(nrow(TF))[Restrict == Restrict[s1] &
                                    BF[, u] != BF[s1, u] &
                                    !(seq_len(nrow(TF)) %in% unavailable)]
            if (length(z) == 0 & counter < 501)
              next
            else if (length(z) == 0 & counter > 500)
              break
            if (length(z) > 1)
              s = c(s1, sample(z, 1))
            else
              s = c(s1, z)
            testDswap = dMat(TM, BM, V, s)[2, 1, drop = FALSE]
            if (testDswap < tol & counter < 501)
              next
            else if (testDswap < tol & counter > 500)
              break
            Dswap = testDswap
            break
          }
          if (counter > 500)
            break
          relD = relD * Dswap
          V = UpDate(V, (TM[s[1],] - TM[s[2],]), (BM[s[2],] - BM[s[1],]), ncol(TM))
          TM[c(s[1], s[2]),] = TM[c(s[2], s[1]), , drop = FALSE]
          TF[c(s[1], s[2]),] = TF[c(s[2], s[1]), , drop = FALSE]
        } #jumps
      } #searches
      TF = globTF
      TM = globTM
    } # levels
    Effic = effBlocks(TF, TM, BF, deg1_BM, deg2_BM)
    }
    Design=data.frame(BF,plots=factor(1:nrow(BF)),treatments=TF)
    return(list(TF = TF, Design = Design, Effic = Effic))
  }
  # *******************************************************************************************************************
  #  main program
  # *******************************************************************************************************************
  if (is.factor(blocks)) blocks=data.frame(blocks)
  if (is.data.frame((blocks))) blocks=list(blocks)
  BF=do.call(cbind, blocks) #bind listed data frames into a single data frame 
  if (!all(sapply(BF, is.factor))) stop("blocks must be factors")
  if (!all(sapply(BF, nlevels)==1))
    BF=BF[,sapply(BF, nlevels)>1,drop=FALSE] else
      BF=data.frame(Main=factor(rep(1,nrow(BF))))
  BF = BF[do.call(order, BF), , drop = FALSE] # makes sure nested blocks are in standard order
  
  if (!is.null(seed))
    set.seed(seed)
  if (!is.numeric(weighting) || weighting < 0 || weighting > 1)
    stop("weighting must be a number between 0 and 1")
  if (!is.numeric(jumps) || jumps < 1 || jumps %% 1 != 0 )
    stop("jumps must be an integer greater than 0")
  if (is.null(searches))
    searches = 1 + 5000 %/% nrow(BF)
  if (!is.numeric(searches) || searches < 1 || searches %% 1 != 0 )
    stop("searches must be an integer greater than 0")
  
  if (is.list(treatments) & !is.data.frame(treatments)) 
    TF=expand.grid(treatments) else
      TF = data.frame(treatments)
  TF = TF[do.call(order, TF), , drop = FALSE]  # makes sure the treatments are in standard order 
  fullrand = unlist(lapply(1:ceiling(nrow(BF) / nrow(TF)), function(j) {
    sample(seq_len(nrow(TF)))
  }))
  TF = TF[fullrand, ,drop = FALSE]
  
  if (!is.null(treatments_model) & !is.list(treatments_model)) 
    treatments_model=list(treatments_model)
  
  if (is.list(treatments_model))   
    treatments_model=sapply(treatments_model,as.formula)
  
  if (is.null(treatments_model))
    treatments_model =list(as.formula( paste("~",paste0(colnames(TF), collapse = "+"))))
  
  # is each block factor nested within the levels of all preceding blocks factors
  nested=
    all(sapply(1:ncol(BF),function(i){nlevels(droplevels(interaction(BF[,1:i])))==nlevels(droplevels(BF[,i]))}))
  
  first_order  = sapply(1:ncol(BF), function(j) {
    as.formula(paste0("~",paste0(colnames(BF)[1:j], collapse = "+")))
  })
  
  second_order = sapply(1:ncol(BF), function(j) {
    as.formula( paste0("~ (", paste0(colnames(BF)[1:j], collapse = "+"), ")^2"))
  })
  
  restriction_model=append(list(NULL),treatments_model)
  EfficsTrts=NULL
  
  for (u in 1:length(treatments_model)) {
    W=fraction(TF,nrow(BF),treatments_model[[u]],restriction_model[[u]])
    TF=W$fullTF
    EfficsTrts=c(EfficsTrts,W$Efficiency)
  }
  TF=W$TF
  
  DF = sapply(1:length(treatments_model), function(i) {
    ncol(model.matrix(as.formula(treatments_model[[i]]), TF)) - 1
  })

  if (nested & ncol(TF) == 1 & is.factor(TF[, 1]) & length(blocks) == 1) {
    Z = nestedBlocks(TF, BF, searches, seed, jumps)
  } else 
    Z = optBlocks(TF, BF) 
  
  blocksModel = Z$Effic
  Design=Z$Design
  names=c(colnames(BF),"plots",colnames(TF))    
  colnames(Design)=names
  TF=Design[, (ncol(Design)-ncol(TF)+1):ncol(Design),drop=FALSE]
  rownames(TF) = NULL
  Replicates = count(TF)
  rownames(Design) = NULL
  rownames(Replicates) = NULL
  # Table for Treatments_model DF and efficiency
  
  treatsModel = data.frame(
    cbind(
      "Treatment model" = treatments_model,
      "Model DF" = DF ,
      "D-Efficiency" = EfficsTrts
    )
  )
  list(
    Replication = Replicates,
    Design = Design,
    Treatments_model = treatsModel,
    Blocks_model = blocksModel,
    weighting = weighting,
    seed = seed,
    searches = searches,
    jumps = jumps
  ) 
  
}
