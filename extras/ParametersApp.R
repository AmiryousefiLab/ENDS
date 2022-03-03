
# Example for shinny app

block = df_example

input = list(mean_switch=T, outlier_switch=T, onehunda_switch=T, dosedep_auc=T,
             checkgroup1 = c("Point samples" ,
                             "Spline fit" ,
                             "Min-max bands",
                             "Empirical viability bands",
                             "Drug span gradient",
                             # "Absolute doses" ,
                             "Relative doses"),
             SplinePlot=T, SigmoidPlot=T, MonotonePlot=T, NPBPlot=F, p_ic=50,
             NPS_title='', file_type='pdf', resolution=300
)


con = 'ModelEstimations.csv'
# Testing out multiple input  functions
mean_switch=T;outlier_switch=T; onehunda_switch=T; dosedependent_auc=T; p_ic=50; NPS_title='title'
check_boxes = c("Point Samples","Spline")

# tbl <- read.csv('/Users/bwilliams/GoogleDrive (syncing)/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/data/Example4.csv', header=T)
# tbl <- read.csv('/Users/bwilliams/GoogleDrive (syncing)/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/data/Example1.csv', header=T)
# block = create_blocks(tbl)
# block2 = preprocess_data_mult(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch)
# p1 = PlotOverlay_mult(block2, check_boxes, dosedependent_auc, p_ic, NPS_title)
# p1 = plot_npbFit_mult(block2, dosedependent_auc, p_ic, NPS_title)


#### Testing MH-npB
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


# Correct increasing fit for npS
# tbl <- read.csv('/Users/bwilliams/GoogleDrive (syncing)/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/data/d1.csv', header=T)
# block = create_blocks(tbl)
# block2 = preprocess_data_mult(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch)
# p1 = PlotOverlay_mult(block2, check_boxes, dosedependent_auc, p_ic, NPS_title)

