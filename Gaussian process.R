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
g <- g + geom_abline(intercept = 0, slope = 1) + xlim(-0.75,1) + ylim(-0.75,1)

plot(g)

## 1D Example 1
n = 5; d = 1;
computer_simulator <- function(x){
  x = 2*x+0.5;
  y = sin(10*pi*x)/(2*x) + (x-1)^4;
  return(y)
}
set.seed(3);
library(lhs);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y)
plot(GPmodel)
##