### AM-Splines
setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS')
source('PreliminaryFunctions.R')
source('NPDS.R')
source('IsotonicRegressionFit.R')
source('SigmoidFit.R') 


lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

posterior_predictive = function(x, C, a){
  y_est = c()
  if(length(x)==1){
    y_est = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )
  }else{
    for(i in  1:length(x))
      y_est[i] = C  + (1-C)*sum(a*(1-pnorm(x[i], K, sqrt(lambda)) ) )  
  }
  return(y_est)
}


make_plots = function(posterior,K,lambda, logplotx=T){
  file_plot = paste('images/npb_stanchains',length(K), lambda, '.png', sep = '_')
  
  C_post = posterior$C
  sigma2_post = posterior$sigma2
  a_post = posterior$a
  n_knots = ncol(a_post)
  
  n = length(C_post)
  obs = seq(1,n, by = 100)
  
  
  library(randomcoloR)
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
  
  x_est = seq(min(x),max(x), length.out = 100)
  y_est = posterior_predictive(x_est, C_est, a_est)
  
  file_plot = paste('images/npb_stanfit',length(K),lambda, '.png', sep = '_')
  png(file=file_plot,res=400, width=20, height=16, units='cm')
    par(mfrow=c(1,1))
    main_plot = paste('NPB fit K=',length(K),' Lambda=',lambda, sep='')
    if(logplotx==T){plot(x,y[,1], ylim = c(0.4,1.1), main=main_plot, log='x')
    }else{
      plot(x,y[,1], ylim = c(0.4,1.1), main=main_plot)
      }
    points(x,y[,2])
    points(x,y[,3])
    points(x,y[,4])
    lines(x_est,y_est, type='l', col='red')
    lines(x, rowMeans(y), col='blue')
  dev.off()
}


input = list(mean_switch=T, outlier_switch=T, onehunda_switch=T, dosedep_auc=T,
             checkgroup1 = c("Point samples" ,
                             "Spline fit" ,
                             "Min-max bands",
                             "Empirical viability bands",
                             "Drug span gradient",
                             # "Absolute doses" ,
                             "Relative doses")
)

df = openxlsx::read.xlsx('data/Drug_response_S8.xlsx', sheet = 1)
df_list  = read_excel_allsheets('data/Drug_response_S8.xlsx')

block = extract_dose_block(df_list, '5FU', 'P1', 'T1', 1)
block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch)



options(mc.cores = parallel::detectCores())
library(rstan)
# Log values for dir(1,..,2)
# K = log10(block2$doses)
# x = log10(block2$doses)
K = block2$doses
x = block2$doses
lambda = 0.1
y_max = max(as.matrix(block2[,2:5]))
y = block2[2:6]/y_max
list_am <- list(
                K = K,
                lambda = lambda,
                n_knots = length(K),
                x=x,
                y=y,
                n = length(x),
                m = ncol(y),
                alpha=1
                )
setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/extras')
# fit <- stan('AMSpline.stan', data = list_am, iter = 1e4)
fit <- stan('AMSpline_stickprior.stan', data = list_am, iter = 1e4)
posterior <- extract(fit)
make_plots(posterior,K,lambda, logplotx=F)


# 9 Plots
K1 =  lseq(from=min(x), to=max(x), length.out = 3)
K2 = lseq(from=min(x), to=max(x), length.out = 10)
K3= lseq(from=min(x), to=max(x), length.out = 30)
Knots = list(K1, K2, K3)
lambdas = c(0.1, 2, 20)


y_max = max(as.matrix(block2[,2: 5]))
x = block2$doses
y = block2[,2:5]/y_max
for(i in 1:3){
  for(j in 1:3){
    set.seed(42)
    K = Knots[[i]]
    lambda = lambdas[[j]]
    list_am <- list(
      K = K,
      lambda = lambda,
      n_knots = length(K),
      x=x,
      y=y,
      n = length(x),
      m = ncol(y),
      alpha=1
    )
    # fit <- stan('AMSpline.stan', data = list_am, iter = 1e4)
    fit <- stan('AMSpline_stickprior.stan', data = list_am, iter = 1e4)
    posterior <- extract(fit)
    
    make_plots(posterior, K  ,lambda )
  }
}

