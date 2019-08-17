#Reference:
#http://statmodeling.hatenablog.com/entry/r-tanimoto-distance
#the function calculates tanimoto coefficient, in order to define the distance you have to use inverse.

tanimoto <- function(x)　{
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

hist <- hclust(dist(multi.regression.compounds),method = "ward.D2")
plot(hist)

# 階層型クラスタリングからクラスタ数を指定して分類する（コサイン類似度、ウォード法）。
clusters <- cutree(hist, k = 6)

# クラスタiの中心点を算出
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

# 全クラスタに適用
centers <- sapply(unique(clusters), clust.centroid, multi.regression.compounds, clusters)

# k-meansの実行。sapply()の結果だと行と列が逆なので、転置して引数に与える
km <- kmeans(multi.regression.compounds, centers=t(centers)) 

