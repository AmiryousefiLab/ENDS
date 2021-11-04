# When we only have one sample replicate
logistic = function(x){
  return(  (1+exp(-x))^{-1} )
}

dlogistic  = function(x){
  return( exp(-x)/((1+exp(-x))^2) )
}
invlogistic  = function(y){
  return( log(y / (1-y)) )
}

library(knitr)
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

dsoftmax = function(s){
  # Return the Jacobian matrix of  s.
  sm = softmax(s)
  dsm = diag(sm) - outer(sm, sm)
  return(dsm)
}


likelihood = function(param){
  C_phi = param[1]
  sigma_phi = param[2]
  a_phi = param[3:(3+n_knots-1)]
  
  C = logistic(C_phi)
  sigma2 = exp(sigma_phi)
  a = softmax(a_phi)
  
  sl = rep(0, length(x))
  for(i  in 1:length(x)){
    pred = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )
    singlelikelihoods = dnorm(as.numeric(y[i]) ,
                              mean = pred,
                              sd = sqrt(sigma2), 
                              log = T)
    sl[i] = singlelikelihoods
  }
  
  sumll = sum(sl)
  return(sumll)
}


# Prior distribution
library(DirichletReg)
prior = function(param){
  # work with transformation of parameters
  C_phi = param[1]
  sigma_phi = param[2]
  a_phi = param[3:(3+n_knots-1)]
  
  C = logistic(C_phi)
  sigma2 = exp(sigma_phi)
  a = softmax(a_phi)
  
  C_prior = dbeta(C, shape1=1, shape2=1, log = T)
  sigma2_prior = dchisq(sigma2, df = 2, log = T)
  a_prior = DirichletReg::ddirichlet( matrix(a+1,ncol=n_knots), alpha = rep(1, n_knots) , log=T)    
  
  # p(a_{i+1}|a_{i}) ~ Unif(0, a_{i}) so p(a_1,..,a_n) = p(a_n|a_{n-1})..p(a_2|a_{1})p(a_1)
  # Not only this but sampling needs to be changed as well, for example with the beta contruction would work, 
  # this is, sample from a normal to reject a beta and then construct stick-breaking process
  # if(stick_breaking == TRUE){
  #   n_knots = length(a)
  #   a_prior_vec = rep(0, n_knots)
  #   a_prior_vec[1] = dunif(a[1],0,1, log = T)
  #   for(i in 2:n_knots){
  #     a_prior_vec[i] = dunif(a[i], 0, a[i-1], log=T )
  #   }
  #   a_prior = sum(a_prior_vec)
  # }
  
  
  dC_prior = log(dlogistic(C_phi)  )
  dSigma2_prior = log(exp(sigma_phi)  )
  da_prior = log(abs(det(dsoftmax(a_phi))))
  
  return(C_prior+sigma2_prior+a_prior+dC_prior+dSigma2_prior+da_prior)
}

posterior = function(param){
  return (likelihood(param) + prior(param))
}


######## Metropolis algorithm ################
proposalfunction = function(param){
  temp = rnorm(2+n_knots,
               mean = param,
               sd= c(0.1,0.1,rep(0.05, n_knots)))
  return(temp)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2+n_knots))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}


chain_NPB = function(x,y,K,lambda, iter = 10000 ){
  n_knots <<- length(K)
  
  startvalue = c(0.5,1, rep(1, n_knots)/n_knots)
  
  chain = run_metropolis_MCMC(startvalue, iter)
  burnIn = round(iter/2)
  acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
  print(paste('Acceptance rate:', acceptance))
  return(chain)
}

# Posterior predictive
posterior_predictive = function(x, C, a){
  y_est = c()
  if(length(x)==1){
    y_est = C  + (1-C)*sum(a*(1-pnorm(x, K, sqrt(lambda)) ) )
  }else{
    for(i in  1:length(x))
      y_est[i] = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )  
  }
  return(y_est)
}


lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

parameter_mean_estimation = function(chain){
  
  n_knots = length(K)
  
  chain_transformed = chain 
  chain_transformed[,1] = logistic(chain[,1])
  chain_transformed[,2] = exp(chain[,1])
  chain_transformed[,3:(3+n_knots-1)] = t(apply(chain[,3:(3+n_knots-1)], softmax, MARGIN = 1))
  
  C_post = chain_transformed[,1]
  sigma2_post = chain_transformed[,2]
  a_post = chain_transformed[,3:(3+n_knots-1)]
  
  n = length(C_post)
  obs = seq(1,n, by = 100)
  n_obs = length(obs)
  obs_drop = obs[floor(n_obs/2):n_obs]
  
  a_est = colMeans(a_post[obs_drop,])
  C_est = mean(C_post[obs_drop])
  sigma2_est = mean(sigma2_post[obs_drop])
  
  x_est = seq(min(x),max(x), length.out = 100)
  y_est = posterior_predictive(x_est, C_est, a_est)
  
  param_est = list(C_est, sigma2_est, a_est, x_est, y_est)
  return(param_est)
}



make_plots_example = function(x,y, chain, K, lambda, logplotx=T, title='example'){
  
  n_knots = length(K)
  
  chain_transformed = chain 
  chain_transformed[,1] = logistic(chain[,1])
  chain_transformed[,2] = exp(chain[,1])
  chain_transformed[,3:(3+n_knots-1)] = t(apply(chain[,3:(3+n_knots-1)], softmax, MARGIN = 1))
  
  C_post = chain_transformed[,1]
  sigma2_post = chain_transformed[,2]
  a_post = chain_transformed[,3:(3+n_knots-1)]
  
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

am_spline = function(x, C, a, lambda, K){
  y_est = c()
  if(length(x)==1){
    y_est = C  + (1-C)*sum(a*(1-pnorm(x, K, sqrt(lambda)) ) )
  }else{
    for(i in  1:length(x))
      y_est[i] = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )  
  }
  return(y_est)
}


# Run chains
# Different K, different Lambdas
setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/extras')
set.seed(42)

# Lets simulate data and try to fit 
x = 1:100/10
a = rev(c(1/2, 1/4, 1/4))
K = c(1,3, 6)
C= 0.5
lambda=0.1
set.seed(41)
y = am_spline(x,C,a, lambda, K) + rnorm(100,0,0.01)


chain = chain_NPB(x,y,K,lambda, iter = 10000)
make_plots_example(x,y,chain,K,lambda, logplotx=F,'dir_mh')

param_est = parameter_mean_estimation(chain)
C_est = param_est[[1]]
sigma2_est = param_est[[2]]
a_est = param_est[[3]]
x_est = param_est[[4]]
y_est = param_est[[5]]

# Functions to obtain AUC and IC50
posterior_predictive_integrate = function(x) posterior_predictive(x, C_est, a_est)
posterior_predictive_inverse = function(x){ posterior_predictive(x, C_est, a_est)- (C_est + (1-C_est)*0.5) }
auc = integrate(posterior_predictive_integrate, lower = min(x), max(x))$value
uniroot = uniroot(posterior_predictive_inverse, interval = c(min(x), max(x)), tol=1e-9)
x_ic50 = uniroot$root
mse = sample_meansquarederror(posterior_predictive_integrate(x), y)

plot(x,y)
curve(posterior_predictive_integrate, from=min(x), to=max(x), add=T, col='red')
abline(v=x_ic50, col='red')
abline(h=C_est + (1-C_est)*0.5 , col='red')

