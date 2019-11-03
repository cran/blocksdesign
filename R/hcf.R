#' @title Finds hcf of any set of positive integers
#' 
#' @description
#' Finds the highest common factor (hcf) of a set of integer numbers greater than zero (Euclidean algorithm).
#' 
#' @details
#' Finds the hcf of any set of positive integers which can be in any order. 
#' 
#' @param ... any set of positive integers, in any order, for which the hcf is required.
#' 
#' @return
#' hcf
#'  
#' @examples
#' 
#' # hcf of vectors of integers
#' HCF(56,77,616)
#' HCF(3,56,77,616) 
#' 
#' @export
HCF=function(...)  {
  v = sort( c(...))
  for (i in  seq_len(length(v))) 
    while ( v[i]%%v[1]!=0  )
      v[c(1,i)] = c(v[i]%%v[1], v[1])
    return(v[1])
}
