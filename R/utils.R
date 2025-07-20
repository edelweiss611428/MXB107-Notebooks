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


empiricalRuleGaussian = function(data,   xlim = c(min(data), max(data))) {
  
  # Empirical mean and SD
  empMean = mean(data)
  empSd = sd(data)
  
  # Define intervals: ±1 SD, ±2 SD, ±3 SD
  intervals = list(
    "±1 SD" = c(empMean - empSd, empMean + empSd),
    "±2 SD" = c(empMean - 2 * empSd, empMean + 2 * empSd),
    "±3 SD" = c(empMean - 3 * empSd, empMean + 3 * empSd)
  )
  
  # Plot histogram with density curve
  hist(data, breaks = 25, probability = TRUE,
       main = "Histogram with ±1 SD, ±2 SD, and ±3 SD Intervals",
       xlab = "Value", col = "lightgray", border = "white", xlim = xlim)
  curve(dnorm(x, mean = empMean, sd = empSd), add = TRUE, col = "blue", lwd = 2)
  
  # Add vertical lines for intervals
  cols = c("red", "green", "purple")
  ltys = c(2, 3, 4)
  
  for (i in seq_along(intervals)) {
    abline(v = intervals[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  }
  
  legend("topright", legend = names(intervals), col = cols, lty = ltys, lwd = 2, bty = "n")
  
  # Compute and print empirical coverage
  coverage = sapply(intervals, function(bounds) {
    mean(data >= bounds[1] & data <= bounds[2])
  })
  
  # Print results
  cat("Empirical coverage:\n")
  for (name in names(coverage)) {
    cat(sprintf("%s: %.2f%%\n", name, 100 * coverage[[name]]))
  }
}


chebyshevRule = function(data, xlim = c(min(data), max(data))) {
  n = length(data)
  empMean = mean(data)
  empSd = sd(data)
  
  ks = c(1, 2, 3)
  intervals = lapply(ks, function(k) c(empMean - k * empSd, empMean + k * empSd))
  names(intervals) = paste0("±", ks, " SD")
  
  # Plot histogram
  hist(data, breaks = 25, probability = TRUE,
       main = "Histogram with Chebyshev ±k SD Intervals",
       xlab = "Value", col = "lightgray", border = "white", xlim = xlim)
  
  cols = c("red", "green", "purple")
  ltys = c(2, 3, 4)
  
  for (i in seq_along(intervals)) {
    abline(v = intervals[[i]], col = cols[i], lwd = 2, lty = ltys[i])
  }
  
  legend("topright", legend = names(intervals), col = cols, lty = ltys, lwd = 2, bty = "n")
  
  # Calculate empirical coverage
  empirical_coverage = sapply(intervals, function(bounds) {
    mean(data >= bounds[1] & data <= bounds[2])
  })
  
  # Chebyshev lower bounds (valid for k >= 1)
  chebyshev_bounds = sapply(ks, function(k) 1 - 1 / k^2)
  
  # Print results
  cat("Coverage vs. Chebyshev Lower Bound:\n")
  for (i in seq_along(ks)) {
    cat(sprintf("k = %d: Empirical = %.2f%%, Chebyshev lower bound = %.2f%%\n",
                ks[i], 100 * empirical_coverage[i], 100 * chebyshev_bounds[i]))
  }
}

