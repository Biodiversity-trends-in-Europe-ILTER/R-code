#####################################################################################################
# Modified version of the function mmkh (R package: modifiedmk; Patakamuri & O’Brien, 2019. Modified versions of Mann Kendall and Spearman’s Rho trend tests.)
# The only difference from the original function is at line 87, to include S and n in the output
#####################################################################################################

My.mmkh=function (x, ci = 0.95) 
{
  x = x
  z = NULL
  z0 = NULL
  pval = NULL
  pval0 = NULL
  S = 0
  Tau = NULL
  essf = NULL
  ci = ci
  if (is.vector(x) == FALSE) {
  stop("Input data must be a vector")
  }
  if (any(is.finite(x) == FALSE)) {
  x <- x[-c(which(is.finite(x) == FALSE))]
  warning("The input vector contains non-finite numbers. An attempt was made to remove them")
  }
  n <- length(x)
  V <- rep(NA, n * (n - 1)/2)
  k = 0
  for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
  k = k + 1
  V[k] = (x[j] - x[i])/(j - i)
  }
  }
  slp <- median(V, na.rm = TRUE)
  t = 1:length(x)
  xn <- (x[1:n]) - ((slp) * (t))
  for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
  S = S + sign(x[j] - x[i])
  }
  }
  ro <- acf(rank(xn), lag.max = (n - 1), plot = FALSE)$acf[-1]
  sig <- qnorm((1 + ci)/2)/sqrt(n)
  rof <- rep(NA, length(ro))
  for (i in 1:(length(ro))) {
  if (ro[i] > sig || ro[i] < -sig) {
  rof[i] <- ro[i]
  }
  else {
  rof[i] = 0
  }
  }
  cte <- 2/(n * (n - 1) * (n - 2))
  ess = 0
  for (i in 1:(n - 1)) {
  ess = ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
  }
  essf = 1 + ess * cte
  var.S = n * (n - 1) * (2 * n + 5) * (1/18)
  if (length(unique(x)) < n) {
  aux <- unique(x)
  for (i in 1:length(aux)) {
  tie <- length(which(x == aux[i]))
  if (tie > 1) {
  var.S = var.S - tie * (tie - 1) * (2 * tie + 
  5) * (1/18)
  }
  }
  }
  VS = var.S * essf
  if (S == 0) {
  z = 0
  z0 = 0
  }
  if (S > 0) {
  z = (S - 1)/sqrt(VS)
  z0 = (S - 1)/sqrt(var.S)
  }
  else {
  z = (S + 1)/sqrt(VS)
  z0 = (S + 1)/sqrt(var.S)
  }
  pval = 2 * pnorm(-abs(z))
  pval0 = 2 * pnorm(-abs(z0))
  Tau = S/(0.5 * n * (n - 1))
  return(c("Corrected Zc" = z, "new P-value" = pval, "N/N*" = essf, 
  "Original Z" = z0, "old P.value" = pval0, "Tau" = Tau, 
  "Sen s slope" = slp, "old.variance" = var.S, "new.variance" = VS, "S statistic" = S, "n" = n))
}
