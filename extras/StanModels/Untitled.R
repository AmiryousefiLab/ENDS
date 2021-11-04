

find_optimal_lambda = function(x,y,K, iter = 100, lambdas = c(0.01, 0.1, 2, 5)){
  # Find optimal lambda parameter in fit in terms of squared errors
  # Idea run small chains for several values of lambda, choose one with smallest squared error. 
  # Lambda options 0.01, 0.1, 2, 5 and mean(sd(doses))
  
  if(!is.vector(y) | dim(y)[2]>1){
    labda_ml = mean(apply(y, MARGIN = 1, sd))  
  }else{
    labda_ml = sd(y)
  }
  
  MSEs = c()
  lambdas = c(lambdas, labda_ml)
  for(lambda in lambdas){
    print(lambda)
    start = Sys.time()
    chain = chain_NPB(x,y,K,lambda, iter = iter)  
    parameter_mean_estimation(chain, dropout=1)
    C_est = param_est[[1]]
    sigma2_est = param_est[[2]]
    a_est = param_est[[3]]
    x_est = param_est[[4]]
    y_est = param_est[[5]]
    posterior_predictive_integrate = function(x) am_spline(x, C_est, a_est, lambda, K, y_max)
    mse = sample_meansquarederror(posterior_predictive_integrate(x), y)
    MSEs = c(MSEs, mse)
    end = Sys.time()
    
    par(mfrow=c(1,1))
    main_plot = paste('NPB fit K=',length(K),' Lambda=',lambda, sep='')
    plot(x,y_og[,1], ylim = c(min(y_og)-10,max(y_og)), main=main_plot, log='x')
    points(x, y_og[,2])
    points(x,y_og[,3])
    points(x,y_og[,4])
    lines(x, rowMeans(y_og), col='blue')
    curve(posterior_predictive_integrate(x), col='red', from = min(x), to=max(x) ,add=T)
    
    print(end-start)
  }
  
  min_id =which.min(MSEs)
  lambda_selected = lambdas[min_id]
  return(lambda_selected)
}


