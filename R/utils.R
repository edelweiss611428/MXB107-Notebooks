skewness = function(x) {
  
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

empiricalRuleGaussian = function(data, xlim = c(min(data), max(data)), ks = c(1, 2, 3)) {
  emp_mean = mean(data)
  emp_sd = sd(data)
  
  # Create ±k*SD intervals
  intervals = lapply(ks, function(k) c(emp_mean - k * emp_sd, emp_mean + k * emp_sd))
  names(intervals) = paste0("±", ks, " SD")
  
  # Plot histogram and Gaussian curve
  hist(data, breaks = 25, probability = TRUE,
       main = "Histogram with Gaussian ±k SD Intervals",
       xlab = "Value", col = "lightgray", border = "white", xlim = xlim)
  
  curve(dnorm(x, mean = emp_mean, sd = emp_sd),
        add = TRUE, col = "blue", lwd = 2)
  
  # Colors and line types
  cols = c("red", "green", "purple", "orange", "brown")[seq_along(ks)]
  ltys = 2:(1 + length(ks))
  
  # Add interval lines and labels
  for (i in seq_along(intervals)) {
    bounds = intervals[[i]]
    abline(v = bounds, col = cols[i], lwd = 2, lty = ltys[i])
    
    # Add numeric labels at base of histogram
    text(x = bounds[1], y = 0, labels = sprintf("%.2f", bounds[1]), 
         col = cols[i], pos = 4, cex = 0.8, offset = 0.2)
    text(x = bounds[2], y = 0, labels = sprintf("%.2f", bounds[2]), 
         col = cols[i], pos = 2, cex = 0.8, offset = 0.2)
  }
  
  # Add legend with intervals and ranges
  legend("topright",
         legend = paste0(names(intervals),
                         "\n[", sapply(intervals, function(b) sprintf("%.2f – %.2f", b[1], b[2])), "]"),
         col = cols, lty = ltys, lwd = 2, bty = "n", text.col = cols)
  
  # Calculate empirical coverage
  empirical_coverage = sapply(intervals, function(bounds) {
    mean(data >= bounds[1] & data <= bounds[2])
  })
  
  # Print coverage results
  cat("Empirical coverage:\n")
  for (i in seq_along(ks)) {
    cat(sprintf("k = %d (±%.0f SD): %.2f%% of data\n",
                ks[i], ks[i], 100 * empirical_coverage[i]))
  }
  
  invisible(list(
    mean = emp_mean,
    sd = emp_sd,
    ks = ks,
    intervals = intervals,
    empirical_coverage = empirical_coverage
  ))
}

chebyshevRule = function(data, xlim = c(min(data), max(data)), ks = c(1, 2, 3)) {
  n = length(data)
  emp_mean = mean(data)
  emp_sd = sd(data)
  
  # Construct intervals
  intervals = lapply(ks, function(k) c(emp_mean - k * emp_sd, emp_mean + k * emp_sd))
  names(intervals) = paste0("±", ks, " SD")
  
  # Plot histogram
  hist(data, breaks = 25, probability = TRUE,
       main = "Histogram with Chebyshev ±k SD Intervals",
       xlab = "Value", col = "lightgray", border = "white", xlim = xlim)
  
  # Colors and line types for overlays
  cols = c("red", "green", "blue", "purple", "orange")[seq_along(ks)]
  ltys = 2:(1 + length(ks))
  
  # Draw vertical lines and annotate intervals
  for (i in seq_along(intervals)) {
    bounds = intervals[[i]]
    abline(v = bounds, col = cols[i], lwd = 2, lty = ltys[i])
    
    # Add text showing interval bounds
    text(x = bounds[1], y = 0, labels = sprintf("%.2f", bounds[1]), 
         col = cols[i], pos = 4, cex = 0.8, offset = 0.2)
    text(x = bounds[2], y = 0, labels = sprintf("%.2f", bounds[2]), 
         col = cols[i], pos = 2, cex = 0.8, offset = 0.2)
  }
  
  legend("topright",
         legend = paste0(names(intervals), 
                         "\n[", sapply(intervals, function(b) sprintf("%.2f – %.2f", b[1], b[2])), "]"),
         col = cols, lty = ltys, lwd = 2, bty = "n", text.col = cols)
  
  # Calculate empirical coverage
  empirical_coverage = sapply(intervals, function(bounds) {
    mean(data >= bounds[1] & data <= bounds[2])
  })
  
  # Chebyshev lower bounds
  chebyshev_bounds = sapply(ks, function(k) if (k >= 1) 1 - 1 / k^2 else NA)
  
  # Print results
  cat("Coverage vs. Chebyshev Lower Bound:\n")
  for (i in seq_along(ks)) {
    cat(sprintf("k = %d: Empirical = %.2f%%, Chebyshev bound = %.2f%%\n",
                ks[i], 100 * empirical_coverage[i], 100 * chebyshev_bounds[i]))
  }
  
  invisible(list(
    mean = emp_mean,
    sd = emp_sd,
    ks = ks,
    intervals = intervals,
    empirical_coverage = empirical_coverage,
    chebyshev_bounds = chebyshev_bounds
  ))
}



rangeBasedSD = function(x){
  
  range = max(x) - min(x)
  return(range/4)
  
}

IQRBasedSD = function(x){
  
  iqr = quantile(x, 0.75) - quantile(x, 0.25) 
  return(iqr/1.349)
  
}

boxPlotDescribe = function(){

  set.seed(123)
  data = c(rnorm(100, mean = 10, sd = 2), 20, 22)  # two minor outliers
  

  boxplot(data,
          main = "Boxplot Example: Distribution with Minor Outliers",
          ylab = "Value",
          col = "lightblue")
  

  text(x = 1.2, y = median(data), labels = "Median (Q2)", pos = 4, col = "blue")
  text(x = 0.8, y = quantile(data, 0.25), labels = "Q1 (25%)", pos = 2, col = "darkgreen")
  text(x = 0.8, y = quantile(data, 0.75), labels = "Q3 (75%)", pos = 2, col = "darkgreen")
  text(x = 1, y = min(data[data > quantile(data, 0.25) - 1.5*IQR(data)]),
       labels = "Lower whisker", pos = 3, cex = 0.8)
  text(x = 1, y = max(data[data < quantile(data, 0.75) + 1.5*IQR(data)]),
       labels = "Upper whisker", pos = 3, cex = 0.8)
  text(x = 1.1, y = 20, labels = "Outlier", pos = 3, col = "red", cex = 0.8)
  text(x = 1.1, y = 22, labels = "Outlier", pos = 3, col = "red", cex = 0.8)
  
  
}


coinToss = function(nTimes, nRepeats) {
  results = numeric(nRepeats)
  
  for (i in 1:nRepeats) {
    tosses = sample(c("H", "T"), size = nTimes, replace = TRUE)
    heads = sum(tosses == "H")
    results[i] = heads / nTimes
  }
  
  return(results)
}

generateZTable = function() {

  row_vals = seq(-3.4, 3.4, by = 0.1)  # rows
  col_vals = seq(0, 0.09, by = 0.01)   # columns
  
  z_table = matrix(NA, nrow = length(row_vals), ncol = length(col_vals))

  for (i in seq_along(row_vals)) {
    for (j in seq_along(col_vals)) {
      z = row_vals[i] + col_vals[j]
      z_table[i, j] = pnorm(z)  
    }
  }
  
  rownames(z_table) = sprintf("%.1f", row_vals)
  colnames(z_table) = sprintf("%.2f", col_vals)
  
  return(round(z_table, 4) %>% kable())
}

CI_visualiser = function(seed = NULL){
  
  set.seed(seed)
  
  mu_true  = 2.0
  sigma    = 1.5
  n        = 10
  reps     = 100
  conf_lvl = 0.95
  
  t_ci = function(x, conf = 0.95) {
    n = length(x)
    xbar = mean(x)
    s    = sd(x)
    alpha = 1 - conf
    tcrit = qt(1 - alpha/2, df = n - 1)
    se = s / sqrt(n)
    c(lower = xbar - tcrit * se, upper = xbar + tcrit * se, mean = xbar)
  }
  
  # simulate
  cis = t(sapply(1:reps, function(i) {
    x = rnorm(n, mean = mu_true, sd = sigma)
    t_ci(x, conf_lvl)
  }))
  
  covered = (cis[,1] <= mu_true & mu_true <= cis[,2])
  coverage = mean(covered)
  
  # plot
  plot(NA, xlim = range(cis[,1:2]), ylim = c(1, reps),
       xlab = "Value", ylab = "Iteration",
       main = sprintf("95%% CI Coverage = %.2f (%d/%d)", coverage, sum(covered), reps))
  abline(v = mu_true, lty = 2, col = "blue")
  
  for (i in 1:reps) {
    col = if (covered[i]) "black" else "red"
    segments(cis[i,1], i, cis[i,2], i, col = col, lwd = 2)
    points(cis[i,3], i, col = col, pch = 16, cex = 0.6)  # sample mean
  }
}

generateTTable = function() {

  alpha = c(0.10, 0.05, 0.025, 0.01, 0.005)  # common alpha levels
  prob = 1 - alpha / 2  # for two-sided tests
  
  df = 1:30  # example df from 1 to 30
  
  t_table = sapply(prob, function(p) qt(p, df))
  
  colnames(t_table) = paste0("alpha=", alpha)
  rownames(t_table) = paste0("df=", df)
  
  return(round(t_table, 4) %>% kable())
}



