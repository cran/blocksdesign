#' @title Find hcf
#' 
#' @description
#' Finds the highest common factor (hcf) of a set of integer numbers greater than zero (Euclidean algorithm)
#' 
#' @details
#' Finds the hcf of a vector of positive integers which can be in any order 
#' 
#' @param v is the vector of integers for which the hcf is required (must be integers)
#' 
#' @return
#' hcf
#'  
#' @examples
#' 
#' # hcf of vectors of integers
#' HCF(c(56,77,616))
#' HCF(c(3,56,77,616))  
#' 
#' @export
 HCF=function(v)  {
  v=sort(v)
  for (i in  seq_len(length(v))) 
    while ( v[i]%%v[1]!=0  )
      v[c(1,i)] = c(v[i]%%v[1], v[1])
  return(v[1])
}

