#' @title Prime number test
#' 
#' @description
#' Tests if a given number v is prime and returns TRUE or FALSE
#' 
#' @details
#' primality of any positive integer based on the fact that all primes except 2 and 3 can be expressed as 6k-1 or 6k+1
#' 
#' @param v is the number to be tested for primality
#' 
#' @return
#' logical TRUE or FALSE
#'  
#' @examples
#' 
#' isPrime(731563)
#' isPrime(7315631)
#'  
#' @export
 isPrime=function(v) {
  if (v < 4) return(TRUE)
  if ((v %% 2==0) | (v %% 3==0))  return(FALSE)
  if (v<25) return(TRUE)
  for (i in  6*seq_len(floor((sqrt(v)+1)/6)))
    if ( (v %% (i-1) == 0) | (v %% (i+1) == 0) ) return(FALSE)
  return(TRUE)
}
