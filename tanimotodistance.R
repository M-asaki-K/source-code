#Reference:
#http://statmodeling.hatenablog.com/entry/r-tanimoto-distance
#the function calculates tanimoto coefficient, in order to define the distance you have to use inverse.

tanimoto <- function(x)�@{
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

# �K�w�^�N���X�^�����O����N���X�^�����w�肵�ĕ��ނ���i�R�T�C���ގ��x�A�E�H�[�h�@�j�B
clusters <- cutree(hist, k = 6)

# �N���X�^i�̒��S�_���Z�o
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

# �S�N���X�^�ɓK�p
centers <- sapply(unique(clusters), clust.centroid, multi.regression.compounds, clusters)

# k-means�̎��s�Bsapply()�̌��ʂ��ƍs�Ɨ񂪋t�Ȃ̂ŁA�]�u���Ĉ����ɗ^����
km <- kmeans(multi.regression.compounds, centers=t(centers)) 
