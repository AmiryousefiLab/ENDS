
# setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS')
source('PreliminaryFunctions.R')
source('NPDS.R')
source('IsotonicRegressionFit.R')
source('SigmoidFit.R') 
source('MH_AMspline.R')



nonparaametric_fit = function(block2){
  y_fit = block2$y_mean
  x_fit = block2$doses
  y_ic = 0.5*(max(y_fit)+min(y_fit))
  
  # now the value on the xaxis for which we get that value on the yaxis
  N = length(x_fit)
  x_ic_multiple = c()
  for( i in 1:(N-1)){
    if( (y_fit[i] <= y_ic & y_fit[i+1] >= y_ic ) |
        (y_fit[i+1] <= y_ic & y_fit[i] >= y_ic ) )
    {
      idx1 = i
      idx2 = i+1
      if(idx1!=idx2){
        y1 = y_fit[idx1]
        y2 = y_fit[idx2]
        x1 = x_fit[idx1]
        x2 = x_fit[idx2]
        x_ic_temp = (y_ic-y1)*((x2-x1)/(y2-y1)) + x1
      } else{
        print('IC_50 is not unique')
        x_ic_temp = c(x_fit[idx1], x_fit[idx2])
      }
      x_ic_multiple = c(x_ic_multiple, x_ic_temp)
    }
  }
  
  # Now find the ic_50 closest to ic_50 by monotone regression
  library(fdrtool)
  mono1 = fdrtool::monoreg(x = log10(block2$doses), y = block2$y_mean, type = 'antitonic')
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  y_fit = mono1$yf
  x_fit = (block2$doses)
  
  xy_fit = ic_50_monotonefit(x_fit, y_fit)
  x_ic_mono = xy_fit[[1]]
  y_ic_mono = xy_fit[[2]]
  
  # ic50 is defined as the closest one to ic_50_mono
  x_ic = x_ic_multiple[which.min(abs(x_ic_multiple-x_ic_mono))]
  
  if(dose_dependent_auc==TRUE) 
    auc = line_integral(block2$doses,block2$y_mean)
  
  if(dose_dependent_auc==FALSE) 
    auc = line_integral(1:length(block2$doses), block2$y_mean)
  
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  samples = block2[2:(m+1)]
  y = block2$y_mean
  mse = sample_meansquarederror(y, samples)
  
  list_stats = list( ic50 = x_ic, 
                     mse = mse,
                     auc = auc
                     )
  
  return(list_stats)  
}
monotone_fit = function(block2){
  mono1 = fdrtool::monoreg(x = log10(block2$doses), y = block2$y_mean, type = 'antitonic')
  
  # Gives exatly same fit
  # mono1 = fdrtool::monoreg(x = block2$doses, y = block2$y_mean, type = 'antitonic')
  
  # Only fits on the means, not on bulks of data
  # mono2 = fdrtool::monoreg(x = block2$doses, y = block2[,1:4], type = 'antitonic')
  # plot(mono1)
  # Calculate IC_{50}
  # It would be the middle point of the fit
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  y_fit = mono1$yf
  x_fit = (block2$doses)
  
  
  xy_fit = ic_50_monotonefit(x_fit, y_fit)
  x_ic = xy_fit[[1]]
  y_ic = xy_fit[[2]]
  
  block2_yf = mono1$yf
  
  # Mean squared error for samples is
  samples = block2[2:(m+1)]
  mse = sample_meansquarederror(mono1$yf, samples)
  
  
  if(dose_dependent_auc==TRUE) 
    auc = line_integral(x_fit,y_fit)
  
  if(dose_dependent_auc==FALSE) 
    auc =  line_integral(1:length(x_fit), y_fit)
  
  
  list_stats = list( ic50 = x_ic, 
                     mse = mse,
                     auc = auc
  )
  
  
  return(list_stats)  
}
sigmoid_fit = function(block2){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  # Fit on means
  y = unlist(block2['y_mean'])
  x = unlist(block2['doses'])
  
  # 4 parameter logistic fit,
  m0 <- drm(y ~ x,  fct = LL.4())
  
  
  # Mean squared error for samples is
  samples = block2[2:(m+1)]
  mse = sample_meansquarederror(m0$curve[[1]](x), samples)
  
  # Equation of the fit
  # f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f}         
  
  x_ic = m0$coefficients[4]
  y_ic = 0.5*(m0$coefficients['d:(Intercept)']+ m0$coefficients['c:(Intercept)'])
  if(dose_dependent_auc==TRUE)
    auc = integrate(m0$curve[[1]], lower = min(x), max(x))$value
  if(dose_dependent_auc==FALSE){
    x_ = 1:length(x)
    m1 = drm(y ~ x_,  fct = LL.4())
    auc = integrate(m1$curve[[1]], lower = min(x), max(x))$value
  }
  
  list_stats = list( ic50 = x_ic[[1]], 
                     mse = mse,
                     auc = auc
  )
  
  return(list_stats)  
}
simple_plot = function(block2, list_np, list_mf, list_sf, list_npb, d,p,t,sample){
  y_max = max(as.matrix(block2[,2:5]))
  x = block2$doses
  y = block2[,2:5]/y_max
  y_og =  block2[,2:5]
  K = block2$doses
  lambda = list_npb[['lambda']]
  C_est = list_npb$param_est[[1]]
  a_est = list_npb$param_est[[3]]
  posterior_predictive_integrate = function(x) am_spline(x, C_est, a_est, lambda, K, y_max)  
  
  # Fit on means
  y_sigm = unlist(block2['y_mean'])
  x_sigm = unlist(block2['doses'])
  # 4 parameter logistic fit,
  m0 <- drm(y_sigm ~ x_sigm,  fct = LL.4())
  
  file_plot = paste('images_test/',d,p,t,sample, '.png')
  png(file=file_plot, res=200, width=10, height=15,units = "cm" )
    par(mfrow=c(1,1))
    main_plot = paste('NPB fit K=',length(K),' Lambda=',lambda, sep='')
    plot(x,y_og[,1], main=main_plot, log='x', ylim = c(min(y_og), max(y_og)))
    points(x,y_og[,2])
    points(x,y_og[,3])
    points(x,y_og[,4])
    lines(x, rowMeans(y_og), col='blue')
    curve(posterior_predictive_integrate(x), col='red', from = min(x), to=max(x) ,add=T)
    curve(m0$curve[[1]](x), col='green', add=T)
  dev.off() 
}

# artificial input
input = list(mean_switch=T, outlier_switch=T, onehunda_switch=T, dosedep_auc=T,
             checkgroup1 = c("Point samples" ,
                             "Spline fit" ,
                             "Min-max bands",
                             "Empirical viability bands",
                             "Drug span gradient",
                             "Relative doses")
)


df = openxlsx::read.xlsx('data/Drug_response_S8.xlsx', sheet = 1)
df_list  = read_excel_allsheets('data/Drug_response_S8.xlsx')
dose_dependent_auc <- input$dosedep_auc


drugs = c('5FU','afatinib','akt','doxorubicin','irinotecan','MEK12','nutlin3a')
patients = c('P1','P2', 'P3')
treatments = c('T1','T2','T3', 'T4', 'T5')
samples = 1:5

df_stats = data.frame(matrix(ncol=14, nrow=0))
colnames(df_stats) = c('id','ic50_np','mse_np','auc_np','ic50_mf','mse_mf','auc_mf','ic50_sf','mse_sf','auc_sf','ic50_npb','mse_npb','auc_npb','lambda')

# for( d in drugs){
#   for(p in patients){
#     for(t in treatments){
#       for(sample in samples){
#             print(paste(d,p,t,sample, sep='_'))
#             block = extract_dose_block(df_list, d, p, t, sample)
#             if(!is.null(block)){
#               block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch)
#               list_np = nonparaametric_fit(block2)
#               list_mf = monotone_fit(block2)
#               list_sf = sigmoid_fit(block2)
#               list_npb = npb_fit(block2)
#               id = paste(d,p,t,sample, sep='_')
#               df_stats = rbind( df_stats,
#                                 cbind(id,
#                                       as.data.frame(list_np),
#                                       as.data.frame(list_mf),
#                                       as.data.frame(list_sf),
#                                       as.data.frame(list_npb[1:4]) ) )
#               # make a plot with all fits to check convergence
#               simple_plot(block2, list_np, list_mf, list_sf, list_npb,d,p,t,sample)
#             }
#       }
#     }
#   }
# }

# colnames(df_stats) = c('id','ic50_np','mse_np','auc_np','ic50_mf','mse_mf','auc_mf','ic50_sf','mse_sf','auc_sf','ic50_npb','mse_npb','auc_npb','lambda')
# date <- Sys.Date()
# write.csv(df_stats , paste0('data/Stats_models',date, '.csv'))

# df_stats = read.csv('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS/data/Stats_models2021-11-04.csv')

# Take out outlier (degenerate  fit)
df_stats2 = df_stats %>% 
  mutate(ic50_sf = ifelse(df_stats$ic50_sf>300, NA, ic50_sf ))

# Boxplot
library(tidyr)
library(forcats)

temp = df_stats2[,c('ic50_np','ic50_mf','ic50_sf','ic50_npb')]
names(temp) = c('NPS', 'NPM', 'PL', 'NPB')
df = gather(temp)
p1 = df %>% 
    mutate(key = factor(key,  levels = c('PL', 'NPS', 'NPM', 'NPB'))) %>% 
    ggplot( aes(x=key, y=value, fill=key)) +
    geom_violin( ) +
    geom_boxplot(width=0.1, show.legend = F)+
    scale_fill_brewer(palette="BuPu") +
    coord_cartesian(ylim =  c(0, 30))+ 
    ylab(expression(IC[50])) +
    xlab('Model') +
    labs(fill = "Model") + 
    ggtitle('IC50 per model')
    
temp = df_stats[,c('mse_np','mse_mf','mse_sf','mse_npb')]
names(temp) = c('NPS', 'NPM', 'PL', 'NPB')
df = gather(temp)
p2 = df %>%
  mutate(key = factor(key,  levels = c('PL', 'NPS', 'NPM', 'NPB'))) %>%
  ggplot( aes(x=key, y=value, fill=key))+
  geom_violin( )+
  geom_boxplot(width=0.1, show.legend = F)+
  scale_fill_brewer(palette="BuPu") +
  coord_cartesian(ylim =  c(0, 300))+ 
  ylab('Mean Square Error') +
  xlab('Model') +
  labs(fill = "Model")+
  ggtitle('MSE per model')

temp = df_stats[,c('auc_np','auc_mf','auc_sf','auc_npb')]
names(temp) = c('NPS', 'NPM', 'PL', 'NPB')
df = gather(temp)
p3 = df %>%  
  mutate(key = factor(key,  levels = c('PL', 'NPS', 'NPM', 'NPB'))) %>% 
  ggplot(aes(x=key, y=value, fill=key))+
  geom_violin( )+
  geom_boxplot(width=0.1, show.legend = F)+
  scale_fill_brewer(palette="BuPu") +
  coord_cartesian(ylim =  c(0, 2500))+ 
  ylab('AUC') +
  xlab('Model') +
  labs(fill = "Model")+
  ggtitle('AUC per model')

library(gridExtra)
p = arrangeGrob(p1+theme(plot.title = element_blank()), p2+theme(plot.title = element_blank()),p3+theme(plot.title = element_blank()), ncol=3, top = 'Statistics by model')

# setwd("/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS")
ggsave(paste('extras/images/boxplots_stats_p1_',Sys.Date(),'.png', sep=''), plot=p1, device = 'png',dpi = 400, width = 2*5, height = 6, unit = 'cm')
ggsave(paste('extras/images/boxplots_stats_p2_',Sys.Date(),'.png', sep=''), plot=p2, device = 'png',dpi = 400, width = 2*5, height = 6, unit = 'cm')
ggsave(paste('extras/images/boxplots_stats_p3_',Sys.Date(),'.png', sep=''), plot=p3, device = 'png',dpi = 400, width = 2*5, height = 6, unit = 'cm')
ggsave(paste('extras/images/boxplots_stats_',Sys.Date(),'.png', sep=''), plot=p, device = 'png',dpi = 400, width = 6*5, height = 6, unit = 'cm')

# Test differences parametric and  non parametric
# no normality, so Wilcoxon rank
shapiro.test(df_stats[,10])
# Wilcoxon test for  Montone, Nonparametric vs Logistic
wilcox.test(df_stats[,'ic50_np'], df_stats[,'ic50_sf'], alternative = "two.sided")
# p-value = 0.02585
wilcox.test(df_stats[,'ic50_mf'], df_stats[,'ic50_sf'], alternative = "two.sided")
# p-value = 0.04467
wilcox.test(df_stats[,'ic50_npb'], df_stats[,'ic50_sf'], alternative = "two.sided")
# p-value =  0.4172

wilcox.test(df_stats[,'mse_np'], df_stats[,'mse_sf'], alternative = "two.sided")
# p-value = 3.213e-10
wilcox.test(df_stats[,'mse_mf'], df_stats[,'mse_sf'], alternative = "two.sided")
# p-value = 0.001849
wilcox.test(df_stats[,'mse_npb'], df_stats[,'mse_sf'], alternative = "two.sided")
# p-value = 2.2e-16

wilcox.test(df_stats[,'auc_np'], df_stats[,'auc_sf'], alternative = "two.sided")
# p-value = 0.8674
wilcox.test(df_stats[,'auc_mf'], df_stats[,'auc_sf'], alternative = "two.sided")
# p-value = 0.965
wilcox.test(df_stats[,'auc_npb'], df_stats[,'auc_sf'], alternative = "two.sided")
# p-value = 0.8901



