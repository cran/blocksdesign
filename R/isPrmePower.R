#' @title Prime number test
#' 
#' @description
#' Tests if a given number v is prime power and returns the base prime and the power or -1 
#' 
#' @details
#' tests if a given number N gives a prime Q when raised to the power Q=N**1/i for all i such that Q>1 
#' 
#' @param N is the number to be tested for primality
#' 
#' @return
#' logical TRUE or FALSE
#'  
#' @examples
#' 
#' isPrimePower(10000)
#'  
#' @export
#' 
 isPrimePower=function(N) {
   if (N==1) return(list(base=1,power=1))
  for (i in 1:floor(log2(N)) ) 
    if ( (N **(1/i))%%1 == 0) power=i
  if (power>0) base=N**(1/power)
  if (isPrime(base)) return(list(base=base,power=power))
  else return(list(base=0,power=0))
}
