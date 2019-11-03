#' @title Finds a prime power solution for N, if available.
#' 
#' @description
#' Tests if a given number N is a prime power and returns either the base prime p and power q 
#' or p = 0 and q = 0.
#' 
#' @details
#' Finds the smallest integral solution for s = N**(1/i), which gives the smallest s such that
#' s**i = N. Then, if s is a prime, the number N is a prime power with p = s and q = i.  
#' 
#' @param 
#' N the number to be tested for primality
#' 
#' @return
#' Returns the base prime p and the power q if N is a prime power; otherwise returns p = 0  and q = 0. 
#'  
#' @examples
#' 
#' isPrimePower(10000)
#'  
#' @export
#' 
 isPrimePower=function(N) {
   if (N<1) return(list(base=0,power=0))
  for (newpower in 1:floor(log2(N))) {
    newbase = N**(1/newpower)
    if (abs( newbase - round(newbase)) < .Machine$double.eps^0.75 ) {
      power=newpower
      base=newbase
    }
 } 
  if (isPrime(base)) return(list(base=base,power=power))
  else return(list(base=0,power=0))
}
