
library(GPfit)
x = preprocessed.x1;
y = preprocessed.y1;
GPmodel = GP_fit(x,y);
print(GPmodel)

GPprediction = predict.GP(GPmodel,preprocessed.x2);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;
yhat

plot(yhat,preprocessed.y2)

sigma <- 2*(mse**0.5)

f <- cbind(yhat, preprocessed.y2, sigma)
f <- as.data.frame(f)

library(reshape2)
library(ggplot2)
library(ggsci)

g <- ggplot(f, aes(x = preprocessed.y2, y = yhat))
g <- g + geom_point()
g <- g + geom_errorbar(aes(ymin = yhat - sigma, ymax = yhat + sigma, width = 0.03))
g <- g + scale_fill_nejm()
g <- g + geom_abline(intercept = 0, slope = 1) + xlim(-0.2,1.2) + ylim(-0.75,1.5)

plot(g)

#Partial dependence plot

m <- c(-0.5,-0.2,0.2,0.4,0.6,0.8,1.2,1.5) # valiable
n <- c(0,0,0,0,0,0,0,0) # constant
M.data <-cbind.data.frame(n,n,m,n) #you can choose the valiable manually

GPprediction = predict.GP(GPmodel,xnew=M.data);
yhat = GPprediction$Y_hat;
mse = GPprediction$MSE;
completedata = GPprediction$complete_data;
completedata;
plot(m,yhat)

sigma <- 2*(mse**0.5) # 2 sigma
f <- cbind(yhat, m, sigma)
f <- as.data.frame(f)

g <- ggplot(f, aes(x = m, y = yhat))
g <- g + geom_point()
g <- g + geom_errorbar(aes(ymin = yhat - sigma, ymax = yhat + sigma, width = 0.03))
g <- g + scale_fill_nejm()
g <- g + xlim(-1,2) + ylim(-3,3)

plot(g)
