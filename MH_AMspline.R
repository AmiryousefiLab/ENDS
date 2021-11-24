# Multiple sample replicates per dose
am_spline = function(x, C, a, lambda, K, scale=1){
  y_est = c()
  if(length(x)==1){
    y_est = C  + (1-C)*sum(a*(1-pnorm(x, K, sqrt(lambda)) ) )
  }else{
    for(i in  1:length(x))
      y_est[i] = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )  
  }
  return(scale*y_est)
}

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


likelihood = function(param, x,y,K,lambda){
  C_phi = param[1]
  sigma_phi = param[2]
  a_phi = param[3:(3+n_knots-1)]
  
  C = logistic(C_phi)
  sigma2 = exp(sigma_phi)
  a = softmax(a_phi)
  
  sl = matrix(0, nrow = length(x), ncol =  n_samples )
  for(i  in 1:length(x)){
    pred = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )
    singlelikelihoods = dnorm(as.numeric(y[i,]) ,
                              mean = pred,
                              sd = sqrt(sigma2), 
                              log = T)
    sl[i,] = singlelikelihoods
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
  
  dC_prior = log(dlogistic(C_phi)  )
  dSigma2_prior = log(exp(sigma_phi)  )
  da_prior = log(abs(det(dsoftmax(a_phi))))
  
  return(C_prior+sigma2_prior+a_prior+dC_prior+dSigma2_prior+da_prior)
}

posterior = function(param,x,y,K,lambda){
  return (likelihood(param,x,y,K,lambda) + prior(param))
}


######## Metropolis algorithm ################
proposalfunction = function(param){
  temp = rnorm(2+n_knots,
               mean = param,
               sd= c(0.1,0.1,rep(0.05, n_knots)))
  return(temp)
}

run_metropolis_MCMC = function(startvalue, iter, x,y,K,lambda){
  chain = array(dim = c(iter+1,2+n_knots))
  chain[1,] = startvalue
  for (i in 1:iter){
    proposal = proposalfunction(chain[i,])
    
    probab = exp(posterior(proposal,x,y,K,lambda) - posterior(chain[i,],x,y,K,lambda))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}


chain_NPB = function(x,y,K,lambda, iter = 10000){
  n_samples <<- length(y)
  n_knots <<- length(K)
  
  startvalue = c(0.5,1, rep(1, n_knots)/n_knots)
  
  chain = run_metropolis_MCMC(startvalue, iter, x,y,K,lambda)
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

parameter_mean_estimation = function(chain,K, burnin=0.5, dropout=100){
  
  n_knots = length(K)
  
  chain_transformed = chain 
  chain_transformed[,1] = logistic(chain[,1])
  chain_transformed[,2] = exp(chain[,1])
  chain_transformed[,3:(3+n_knots-1)] = t(apply(chain[,3:(3+n_knots-1)], softmax, MARGIN = 1))
  
  C_post = chain_transformed[,1]
  sigma2_post = chain_transformed[,2]
  a_post = chain_transformed[,3:(3+n_knots-1)]
  
  n = length(C_post)
  obs = seq(1,n, by = dropout)
  n_obs = length(obs)
  obs_drop = obs[floor(n_obs*burnin):n_obs]
  
  a_est = colMeans(a_post[obs_drop,])
  C_est = mean(C_post[obs_drop])
  sigma2_est = mean(sigma2_post[obs_drop])
  
  param_est = list(C_est, sigma2_est, a_est)
  return(param_est)
}

meansd_lambda = function(y){
  if(!is.vector(y) | dim(y)[2]>1){
    lambda_ml = mean(apply(y, MARGIN = 1, sd))  
  }else{
    lambda_ml = sd(y)
  }
  return(lambda_ml)
}

find_optimal_lambda = function(x,y,K,y_max, iter = 100, lambdas = c(0.01, 0.1, 2, 5)){
  # Find optimal lambda parameter in fit in terms of squared errors
  # Idea run small chains for several values of lambda, choose one with smallest squared error. 
  # Lambda options 0.01, 0.1, 2, 5 and mean(sd(doses))
  
  lambda_ml = meansd_lambda(y)
  
  MSEs = c()
  lambdas = c(lambdas, lambda_ml)
  start = Sys.time()
  for(lambda in lambdas){
    
    
    chain = chain_NPB(x,y,K,lambda, iter = iter)  
    param_est = parameter_mean_estimation(chain,K, dropout=1)
    
    C_est = param_est[[1]]
    sigma2_est = param_est[[2]]
    a_est = param_est[[3]]
    
    posterior_predictive_integrate = function(x) am_spline(x, C_est, a_est, lambda, K, y_max)
    mse = sample_meansquarederror(posterior_predictive_integrate(x), y)
    MSEs = c(MSEs, mse)
    end = Sys.time()
    
    # par(mfrow=c(1,1))
    # main_plot = paste('NPB fit K=',length(K),' Lambda=',lambda, sep='')
    # plot(x,y_og[,1], ylim = c(min(y_og)-10,max(y_og)), main=main_plot, log='x')
    # points(x, y_og[,2])
    # points(x,y_og[,3])
    # points(x,y_og[,4])
    # lines(x, rowMeans(y_og), col='blue')
    # curve(posterior_predictive_integrate(x), col='red', from = min(x), to=max(x) ,add=T)
    
    
  }
  print(end-start)
  
  min_id =which.min(MSEs)
  lambda_selected = lambdas[min_id]
  return(lambda_selected)
}

npb_fit = function(block2, dose_dependent_auc=TRUE, p_ic=50){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  
  q_ic = 100-p_ic
  set.seed(42)
  y_max <-  max(as.matrix(block2[,2:m]))
  x <- block2$doses
  y <- block2[,2:m]/y_max
  y_og <-  block2[,2:m]
  K <- block2$doses
  # lambda = 2 
  lambda  <- find_optimal_lambda(x,y,K,y_max, iter = 300, lambdas = c(0.01, 0.1, 2, 5))
  # lambda= meansd_lambda(y)
  chain = chain_NPB(x,y,K,lambda, iter = 1000)
  
  param_est = parameter_mean_estimation(chain,K, burnin = 0.5, dropout = 10)
  C_est = param_est[[1]]
  sigma2_est = param_est[[2]]
  a_est = param_est[[3]]
  
  
  posterior_predictive_integrate = function(x) am_spline(x, C_est, a_est, lambda, K, y_max)
  am_min = am_spline(max(x), C_est, a_est,lambda, K )
  am_max =  am_spline(min(x), C_est, a_est,lambda, K )
  y_ic_unomr = (am_min+(am_max-am_min)*(q_ic/100))
  posterior_predictive_inverse = function(x){ am_spline(x, C_est, a_est,lambda, K )- y_ic_unomr }
  
  
  if(q_ic<1){
    x_ic = max(x)
    y_ic = posterior_predictive_integrate(x_ic)
  }
  if(q_ic>99){
    x_ic = min(x)
    y_ic = posterior_predictive_integrate(x_ic)
  }
  if(1<=q_ic & q_ic<=99){
    uniroot = uniroot(posterior_predictive_inverse, interval = c(min(x), max(x)), tol=1e-9)
    x_ic = uniroot$root
    y_ic = posterior_predictive_integrate(x_ic)
  }
  
  mse = sample_meansquarederror(posterior_predictive_integrate(x), y_og)
  
  if(dose_dependent_auc==TRUE)
    auc = integrate(posterior_predictive_integrate, lower = min(x), max(x))$value
  if(dose_dependent_auc==FALSE){
    x_ = 1:length(x)
    chain = chain_NPB(x_,y,K,lambda, iter = 1000)
    param_est = parameter_mean_estimation(chain,K, burnin = 0.5, dropout = 10)
    C_est = param_est[[1]]
    sigma2_est = param_est[[2]]
    a_est = param_est[[3]]
    posterior_predictive_integrate = function(x) am_spline(x, C_est, a_est, lambda, K, y_max)
    auc = integrate(posterior_predictive_integrate, lower = min(x), max(x))$value
  }
  
  list_stats = list( ic50 = x_ic, 
                     y_ic = y_ic,
                     mse = mse,
                     auc = auc,
                     lambda=lambda,
                     param_est = param_est,
                     posterior_predictive_integrate=posterior_predictive_integrate
  )
  
  return(list_stats)  
}

plot_npbFit = function(block2, dose_dependent_auc=TRUE, p_ic=50, title = ''){
  
  m = dim(block2)[2]-2 
  if(m==0) m <- m+1
  
  # Fit on means
  # y = unlist(block2['y_mean'])
  # x = unlist(block2['doses'])
  
  # NPB fit over whole data
  list_npb = npb_fit(block2, dose_dependent_auc, p_ic)
  
  x_ic = list_npb[['ic50']] 
  posterior_predictive_integrate = list_npb[['posterior_predictive_integrate']] 
  y_ic = posterior_predictive_integrate(x_ic)
  mse = list_npb[['mse']] 
  auc = list_npb[['auc']] 
  lambda=list_npb[['lambda']] 
  param_est = list_npb[['param_est']]
  
  text0 = round(p_ic)
  text1 = max(round(x_ic,2), signif(x_ic,3) )
  text2 = round(mse,2)
  text3 = round(auc)
  
  text = sprintf("atop(atop(IC[%s] == %s, AUC == %s),atop( MSE == %s, \t) )",text0,text1, text3, text2)
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  if(title == '') title = 'Nonparemetric Bayesian'
  p = plot_initialize(block2)
  p = plot_point_samples(p, block2)
  
  colls <<- c(colls, "NPB"="purple", "IC"="red")
  linetypes <<- c(linetypes, "solid", "dotted")
  shapes <<- c(shapes, NA, NA)
  
  
  # We can add scalecolormanual since there is no other layer  to add on top
  p <- p +  
    ggtitle( title ) +
    geom_function(fun = posterior_predictive_integrate, aes(colour='NPB')) + 
    geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic, colour="IC"),  linetype="dotted", show.legend = F) + 
    annotate(geom = 'text', y= y_lim_right, x =max(block2$doses), 
             hjust=1,
             vjust=1,
             label = text, parse =T, size = 7, 
             color='purple') +
    scale_colour_manual(name="Labels",values=colls,
                        guide = guide_legend(
                          override.aes =
                            list(
                              linetype = linetypes,
                              shape = shapes
                            )
                        )
    )
  
  return(p)
}


make_plots = function(chain, K, lambda, logplotx=T, title='example'){
  
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
  if(logplotx==T){
    plot(x,y[,1], ylim = c(C_est-0.1,1.1), main=main_plot, log='x')
  }else{
    plot(x,y[,1], ylim = c(C_est-0.1,1.1), main=main_plot)
  }
  points(x,y[,2])
  points(x,y[,3])
  points(x,y[,4])
  lines(x, rowMeans(y), col='blue')
  lines(x_est,y_est, type='l', col='red')
  dev.off()
  
}



### Load data
# setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS')
# source('PreliminaryFunctions.R')
# source('NPDS.R')
# source('IsotonicRegressionFit.R')
# source('SigmoidFit.R') 
# 
# input = list(mean_switch=T, outlier_switch=T, onehunda_switch=T, dosedep_auc=T,
#              checkgroup1 = c("Point samples" ,
#                              "Spline fit" ,
#                              "Min-max bands",
#                              "Empirical viability bands",
#                              "Drug span gradient",
#                              # "Absolute doses" ,
#                              "Relative doses")
# )
# 
# df = openxlsx::read.xlsx('data/Drug_response_S8.xlsx', sheet = 1)
# df_list  = read_excel_allsheets('data/Drug_response_S8.xlsx')
# 
# block = extract_dose_block(df_list, '5FU', 'P1', 'T1', 1)
# block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch)



# # Run chains
# # Different K, different Lambdas
# setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/extras')
# set.seed(42)
# y_max = max(as.matrix(block2[,2:5]))
# x = block2$doses
# y = block2[,2:5]/y_max
# y_og =  block2[,2:5]
# K = block2$doses
# # lambda = 2 
# lambda  = find_optimal_lambda(x,y,K,y_max, iter = 100, lambdas = c(0.01, 0.1, 2, 5))
# chain = chain_NPB(x,y,K,lambda, iter = 1000)
# make_plots(chain, K, lambda, logplotx=T, title='dir_fit')
# 
# param_est = parameter_mean_estimation(chain,K, burnin = 0.5, dropout = 10)
# C_est = param_est[[1]]
# sigma2_est = param_est[[2]]
# a_est = param_est[[3]]
# x_est = param_est[[4]]
# y_est = param_est[[5]]
# 
# 
# # Functions to obtain AUC and IC50
# posterior_predictive_integrate = function(x) am_spline(x, C_est, a_est, lambda, K, y_max)
# posterior_predictive_inverse = function(x){ posterior_predictive(x, C_est, a_est)- (C_est + (1-C_est)*0.5) }
# uniroot = uniroot(posterior_predictive_inverse, interval = c(min(x), max(x)), tol=1e-9)
# x_ic50 = uniroot$root
# # x_ic50 is invariant to scale, for auc and mse bring back to original scale
# auc = integrate(posterior_predictive_integrate, lower = min(x), max(x))$value
# mse = sample_meansquarederror(posterior_predictive_integrate(x), y)
# 
# x_est = seq(min(x),max(x), length.out = 200)
# y_est = am_spline(x_est, C_est, a_est, lambda, K, scale=y_max)
# 
# par(mfrow=c(1,1))
# main_plot = paste('NPB fit K=',length(K),' Lambda=',lambda, sep='')
# plot(x,y_og[,1], main=main_plot, log='x')
# points(x,y_og[,2])
# points(x,y_og[,3])
# points(x,y_og[,4])
# lines(x, rowMeans(y_og), col='blue')
# lines(x_est,y_est, type='l', col='red')
# 


# 9 Plots
# fit is not good when knots are not on doses
# K1 =  lseq(from=min(x), to=max(x), length.out = 3)
# K2 = lseq(from=min(x), to=max(x), length.out = 10)
# K3= lseq(from=min(x), to=max(x), length.out = 30)
# Knots = list(K1, K2, K3)
# lambdas = c(0.1, 2, 20)
# 
# y_max = max(as.matrix(block2[,2: 5]))
# x = block2$doses
# y = block2[,2:5]/y_max
# for(i in 1:3){
#   for(j in 1:3){
#     set.seed(42)
#     K = Knots[[i]]
#     lambda = lambdas[[j]]
#     chain = chain_NPB(x,y,K,lambda) 
#     make_plots(chain,K  ,lambda )
#   }
# }

# source('PreliminaryFunctions.R')
# source('NPDS.R')
# input = list(mean_switch=T, outlier_switch=F, onehunda_switch=T, dosedep_auc=T,
#              checkgroup1 = c("Point samples" ,
#                              "Spline fit" ,
#                              "Min-max bands",
#                              "Empirical viability bands",
#                              "Drug span gradient",
#                              "Relative doses")
# )
# df = openxlsx::read.xlsx('data/Drug_response_S8.xlsx', sheet = 1)
# df_list  = read_excel_allsheets('data/Drug_response_S8.xlsx')
# df_example = read.csv('data/Example1.csv')
# block =extract_dose_block(df_list, '5FU', 'P1', 'T1', 1)
# block = df_example
# block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch)
# 
# p = plot_npbFit(block2, T, 50)
# p
