#

dhits <- function(data = NULL)
{
  ngene = dim(data)[1]
  pValue = rep(0,ngene)

  for (i in 1:ngene) {
      x = c(data[i,2], data[i,4] - data[i,2])
      y = c(data[i,3], data[i,5] - data[i,3])
      pValue[i] = fisher.test(cbind(x,y))$p
  }

  cbind(data,pValue)
}
