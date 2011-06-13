# example
# chr     pos     gene    qual    c1      c0      n1      n0      g1      g0      hweChisq1       hweChisq0       function        rare
# 1       900559  KLHL17  214.45  13      7       117     211     94/13/0 204/7/0 0.447587930991807       0.0600319349687908      frameshift      0
 

indelFisher <- function(data = NULL)
{
  nvar = dim(data)[1]
  pValue = rep(0,nvar)
  hweP1 = rep(0,nvar)
  hweP0 = rep(0, nvar)
  for (i in 1:nvar) {
      x = c(data[i,5], data[i,7]*2 - data[i,5])
      y = c(data[i,6], data[i,8]*2 - data[i,6])
      pValue[i] = fisher.test(cbind(x,y))$p 
      if (is.na(data[i,11])) {
	hweP1[i] = 1
      } else {
        hweP1[i] = pchisq(data[i,11], df = 1, lower.tail=F)
      }

      if (is.na(data[i,12])) { 
	hweP0[i]  = 1
      } else {
         hweP0[i] = pchisq(data[i,12], df = 1, lower.tail=F)
      } 
  }

  cbind(data,pValue, hweP1, hweP0)
}
