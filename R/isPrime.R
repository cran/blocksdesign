#' @title Prime number test
#' 
#' @description
#' Tests if a given number is prime and returns TRUE or FALSE
#' 
#' @details
#' Tests for the primality of any positive integer using the fact that all primes except 2 and 3 can be
#' expressed as 6k-1 or 6k+1 for integer k.
#' 
#' @param 
#' v the number to be tested for primality
#' 
#' @return
#' logical TRUE or FALSE
#'  
#' @examples
#' 
#' isPrime(731563)
#' isPrime(7315631)
#' isPrime(31**2)
#'  
#' @export
 isPrime=function(v) {
  if (abs( v - round(v)) > .Machine$double.eps^.75) return(FALSE) 
  else v=round(v) # must be integer
  if (v < 4) return(TRUE)
  if ((v %% 2==0) | (v %% 3==0))  return(FALSE)
  if (v<25) return(TRUE)
  for (i in  6*seq_len(floor((sqrt(v)+1)/6)))
    if ( (v %% (i-1) == 0) | (v %% (i+1) == 0) ) return(FALSE)
  return(TRUE)
}

