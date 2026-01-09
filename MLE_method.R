# Load MLE function ####
#
# This custom function to estimate N-M slopes through maximum-likelihood method 
# was created by Ignasi Arranz using original functions and materials from package sizeSpectra.
# See: Edwards, A. M. sizeSpectra (github, 2019); https://github.com/andrew-edwards/sizeSpectra
#
MLE.method <- function (x) {
  #MLE
  xmin = min(x)
  xmax = max(x)
  log.x = log(x)
  sum.log.x = sum(log.x)
  PL.bMLE = 1/(log(min(x)) - sum.log.x/length(x)) - 1
  PLB.minLL = nlm(negLL.PLB, p = PL.bMLE, x = x, n = length(x), 
                  xmin = xmin, xmax = xmax, sumlogx = sum.log.x)
  PLB.bMLE = PLB.minLL$estimate
  PLB.minNegLL = PLB.minLL$minimum
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 1e-05)
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB(bvec[i], x = x, n = length(x), 
                              xmin = xmin, xmax = xmax, sumlogx = sum.log.x)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  hMLE.list = list(b = PLB.bMLE, confVals = PLB.MLE.bConf)
  
  return(hMLE.list = hMLE.list)
  
}