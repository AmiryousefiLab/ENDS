

am_spline = function(x, C, a, lambda, K){
  y_est = c()
  if(length(x)==1){
    y_est = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )
  }else{
    for(i in  1:length(x))
      y_est[i] = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )  
  }
  return(y_est)
}

make_plots_example = function(x,y, C_post, sigma2_post, a_post, K, lambda, logplotx=T, title='example'){
  
  n_knots = length(K)
  
  n = length(C_post)
  obs = seq(1,n, by = 100)
  
  
  library(randomcoloR)
  file_plot = paste('images/npb_chains',title,length(K), lambda, '.png', sep = '_')
  png(file=file_plot, res=200, width=10, height=15,units = "cm" )
    par(mfrow = c(3,1))
    plot(C_post[obs], type='l',col='red',main = 'C (right limit)') 
    plot(sigma2_post[obs], type='l',col='blue',main='Sigma (std)')
    plot( a_post[obs,1], type='l',col='blue',main='a (weights)', ylim=c(0,0.5) )
    for( i in 2:(n_knots-1))
      lines(a_post[obs,i], type='l',col=randomColor(count=1))
  dev.off()
  
  # Final mean estimates, same  procedure as done in paper to construct point process
  a_est = colMeans(a_post)
  C_est = mean(C_post)
  sigma2_est = mean(sigma2_post)
  
  x_est = seq(min(x),max(x), length.out = 200)
  y_est = am_spline(x_est, C_est, a_est, lambda, K)
  
  file_plot = paste('images/npb',title,length(K),lambda, '.png', sep = '_')
  png(file=file_plot,res=400, width=20, height=16, units='cm')
  par(mfrow=c(1,1))
  main_plot = paste('NPB fit K=',length(K),' Lambda=',lambda, sep='')
  if(logplotx==T){plot(x,y, ylim = c(C_est-0.1,1.1), main=main_plot, log='x')
  }else{
    plot(x,y, ylim = c(C_est-0.1,1.1), main=main_plot)
  }
  lines(x_est,y_est, type='l', col='red')
  dev.off()
}


# Lets simulate data and try to fit 
x = 1:100/10
a = rev(c(1/2, 1/4, 1/4))
K = c(1,3, 6)
C= 0.5
lambda=0.1
set.seed(41)
y = am_spline(x,C,a, lambda, K) + rnorm(100,0,0.01)
# plot(x,y)

options(mc.cores = parallel::detectCores())
library(rstan)

list_am <- list(
  K = K,
  lambda = lambda,
  n_knots = length(K),
  x=x,
  y=y,
  n = length(x),
  alpha=1
)
setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/extras')
# fit <- stan('AMSpline.stan', data = list_am, iter = 1e4)
fit1 <- stan('StanExample.stan', data = list_am, iter = 1e4)
posterior1 <- extract(fit1)
C_post = posterior1$C
sigma2_post = posterior1$sigma2
a_post = posterior1$a
make_plots_example(x,y,C_post, sigma2_post, a_post,K,lambda, logplotx=F, title='stick')

fit2 <- stan('StanExampleDirUniPrior.stan', data = list_am, iter = 1e4)
posterior2 <- extract(fit2)
C_post = posterior2$C
sigma2_post = posterior2$sigma2
a_post = posterior2$a
make_plots_example(x,y,C_post, sigma2_post, a_post,K,lambda, logplotx=F,'dir')


# I would like to try also my manual implementation of MH and with stick prior, see if it converges, also increase number of iterations
fit3 <- stan('', data = list_am, iter = 1e4)
posterior2 <- extract(fit2)
make_plots_example(x,y,C_post, sigma2_post, a_post,K,lambda, logplotx=F,'dir')
