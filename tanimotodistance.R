#Reference:
#http://statmodeling.hatenablog.com/entry/r-tanimoto-distance
#the function calculates tanimoto coefficient, in order to define the distance you have to use inverse.

tanimoto <- function(x)@{
  if(is.data.frame(x)) x <- as.matrix(x)
  x[is.na(x)] <- 0
  n <- ncol(x)
  ab <- t(x) %*% x
  aa <- diag(ab)
  aa.m <- matrix(rep(aa,n),n)
  y <- ab/(aa.m + t(aa.m) - ab)
  return(y)
}

View(1-tanimoto(t(multi.regression.compounds[c(1:100),])))
