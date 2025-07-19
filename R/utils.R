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


Mode = function(x) {
  
  if (length(x) == 0) return(NA)
  
  # Remove NA values
  x = x[!is.na(x)]
  if (length(x) == 0) return(NA)
  
  # Compute frequency table
  freq = table(x)
  
  # Find max frequency
  max_freq = max(freq)
  
  # Return all values with max frequency 
  modes = names(freq)[freq == max_freq]
  
  # Convert to the original data type (names(freq) return a character vector)
  if (is.numeric(x)) {
    modes = as.numeric(modes)
  } else if (is.factor(x)) {
    modes = factor(modes, levels = levels(x))
  } else if (is.character(x)) {
    modes = as.character(modes)
  }
  return(modes)
}

FDbinning = function(x) {
  
  if (!is.numeric(x)) stop("x must be numeric.")
  
  if (length(x) == 0) return(NA)
  
  # Remove NA values
  x = x[!is.na(x)]
  if (length(x) == 0) return(NA)
  
  n = length(x)
  
  # Freedman–Diaconis rule
  iqr = IQR(x)
  binwidth = 2 * iqr / n^(1/3)
  
  if (binwidth <= 0){
    stop("IQR is 0!")
  }
  
  # Define breaks
  range_x = range(x)
  breaks = seq(floor(range_x[1]), ceiling(range_x[2]) + binwidth, by = binwidth)
  
  # Cut into labeled intervals
  binned = cut(
    x,
    breaks = breaks,
    right = TRUE,
    include.lowest = TRUE
  )
  
  return(binned)
}

ModeBinMidpoint = function(x) {
  # Get bins using your FDbinning function
  bins = FDbinning(x)
  
  # Count frequency per bin
  tab = table(bins)
  
  if(length(tab) == 0) return(NA)  # no data or no bins
  
  # Find the mode bin (highest frequency)
  mode_bin = names(tab)[which.max(tab)]
  
  # mode_bin is something like "(a,b]"
  # Extract numeric limits from interval string
  # Remove parentheses/brackets and split by comma
  limits_str = gsub("\\(|\\)|\\[|\\]", "", mode_bin)
  limits = as.numeric(strsplit(limits_str, ",")[[1]])
  
  # Calculate midpoint
  midpoint = mean(limits)
  
  return(midpoint)
}