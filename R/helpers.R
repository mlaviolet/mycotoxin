# helper functions

# Wilson binomial confidence interval
# https://rpubs.com/brouwern/binomialCI2
# retrieved 2025-12-28

binom.CI <- function(events, #events = outcomes
                     trials, #number of individuals, test, etc
                     alpha = 0.05)
  {
  n <- trials
  x <- events
  p.hat <- x/n
  
  # Calculate upper and lower limit
  upper.lim <- (p.hat + (qnorm(1-(alpha/2))^2/(2*n)) + qnorm(1-(alpha/2)) * sqrt(((p.hat*(1-p.hat))/n) + (qnorm(1-(alpha/2))^2/(4*n^2))))/(1 + (qnorm(1-(alpha/2))^2/(n)))
  
  lower.lim <- (p.hat + (qnorm(alpha/2)^2/(2*n)) + qnorm(alpha/2) * sqrt(((p.hat*(1-p.hat))/n) + (qnorm(alpha/2)^2/(4*n^2))))/(1 + (qnorm(alpha/2)^2/(n)))
  
  # Modification for probabilities close to boundaries
  
  if ((n <= 50 & x %in% c(1, 2)) | (n >= 51 & n <= 100 & x %in% c(1:3))) {
    lower.lim <- 0.5 * qchisq(alpha, 2 * x)/n
  }
  
  if ((n <= 50 & x %in% c(n - 1, n - 2)) | (n >= 51 & n <= 100 & x %in% c(n - (1:3)))) {
    upper.lim <- 1 - 0.5 * qchisq(alpha, 2 * (n - x))/n
  }
  
  out <- c(lower.lim,upper.lim)
  names(out) <- c("lower.CI","upper.CI")
  return(out)
  
  # The core code for this function is from
  #   stats.stackexchange.com/questions/59733/can-agresti-coull-binomial#-confidence-intervals-be-negative
  #   and was written by stats.stackexchange.com/users/21054/coolserdash
  
  # for more info see 
  #   wikipedia.org/wiki/Binomial_proportion_confidence_interval
  }
