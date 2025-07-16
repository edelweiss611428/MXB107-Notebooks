skewness = function(x) {
  
  if (na.rm){
    x = x[!is.na(x)]
    }
  
  n = length(x)
  if (n < 3){
    return(NA_real_)
  }
  
  m = mean(x)
  s = sd(x)
  
  if (s == 0){
    return(0)  
  }
  
  sum(((x - m)/s)^3) * (n / ((n - 1) * (n - 2)))
}


kurtosis = function(x) {
  
  if (na.rm) x = x[!is.na(x)]
  
  n = length(x)
  if (n < 4) return(NA_real_)
  
  m = mean(x)
  s = sd(x)
  
  if (s == 0) return(0)
  
  # 4th tandardised moment
  mu4 = sum(((x - m)/s)^4) / n
  
  # bias-corrected kurtosis
  ((n*(n+1)) / ((n-1)*(n-2)*(n-3))) * sum(((x - m)/s)^4) -
    (3 * ((n-1)^2) / ((n-2)*(n-3)))
}