
setwd('/Users/bwilliams/GoogleDrive/UniversityOfHelsinki/Summer2021/Network Pharmacology Group/ENDS/ENDS')

source('PreliminaryFunctions.R')
source('NPDS.R')
source('IsotonicRegressionFit.R')
source('SigmoidFit.R') 
source('MH_AMspline.R')

plot_sigmodiMonotoneNpbFit = function(block2, dose_dependent_auc=TRUE){
  
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
  
  x_ic_sigm = m0$coefficients[4]
  y_ic_sigm = 0.5*(m0$coefficients['d:(Intercept)']+ m0$coefficients['c:(Intercept)'])
  if(dose_dependent_auc)
    auc = integrate(m0$curve[[1]], lower = min(x), max(x))$value
  if(dose_dependent_auc==FALSE){
    x_ = 1:length(x)
    m1 = drm(y ~ x_,  fct = LL.4())
    auc = integrate(m1$curve[[1]], lower = min(x), max(x))$value
  }
  
  
  text1 = max(round(x_ic_sigm,2), signif(x_ic_sigm,3) )
  text2 = round(mse,2)
  text3 = round(auc)
  
  text = sprintf("atop(atop(IC[50] == %s, AUC == %s),atop( MSE == %s, \t) )",text1, text3, text2)
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  p = plot_initialize(block2)
  p = plot_point_samples(p, block2)
  
  colls <<- c(colls, "PL"="red", "IC50PL"="red")
  linetypes <<- c(linetypes, "solid", "dotted")
  shapes <<- c(shapes, NA, NA)
  
  
  # We can add scalecolormanual since there is no other layer  to add on top
  p <- p +  
    geom_function(fun = m0$curve[[1]], aes(colour='PL')) + 
    geom_hline( yintercept =  y_ic_sigm, color='red',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic_sigm, colour="IC50PL"),  linetype="dotted", show.legend = F) + 
    annotate(geom = 'text', y= y_lim_right, x =max(block2$doses), 
             hjust=1,
             vjust=1,
             label = text, parse =T, size = 7, 
             color='red') 
    

  mono1 = fdrtool::monoreg(x = log10(block2$doses), y = block2$y_mean, type = 'antitonic')
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  y_fit = mono1$yf
  x_fit = (block2$doses)
  
  
  xy_fit = ic_50_monotonefit(x_fit, y_fit)
  x_ic_mono = xy_fit[[1]]
  y_ic_mono = xy_fit[[2]]
  
  block2_yf = mono1$yf
  
  # Mean squared error for samples is
  samples = block2[2:(m+1)]
  mse = sample_meansquarederror(mono1$yf, samples)
  
  
  if(dose_dependent_auc==TRUE) 
    auc = line_integral(x_fit,y_fit)
  
  if(dose_dependent_auc==FALSE) 
    auc =  line_integral(1:length(x_fit), y_fit)
  
  text1 = max(round(x_ic_mono,2), signif(x_ic_mono,3) )
  text2 = round(mse,2)
  text3 = round(auc)
  text = sprintf("atop(atop(IC[50] == %s, AUC == %s),atop( MSE == %s, \t) )",text1, text3, text2)
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = min(block2[,2:(m+1)], na.rm=T)
  }
  
  
  colls <<- c(colls, "NPM"="darkgreen", "IC50NPM"="darkgreen")
  linetypes <<- c(linetypes, "solid","dotted")
  shapes <<- c(shapes, NA, NA)
  
  options(warn=-1)
  p = p +
    geom_line(aes(block2$doses, block2_yf, colour ='NPM')) + 
    geom_hline( yintercept =  y_ic_mono, color='darkgreen',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic_mono, colour="IC50NPM"),  linetype="dotted", show.legend = F) + 
    annotate(geom = 'text', y= y_lim_right, x =min(block2$doses), 
             hjust=-0.1,
             vjust=-0.1,
             label = text, parse =T, size = 7,
             color = 'darkgreen')
    
  
  
  list_npb = npb_fit(block2, dose_dependent_auc)
  
  x_ic = list_npb[['ic50']] 
  posterior_predictive_integrate = list_npb[['posterior_predictive_integrate']] 
  y_ic = posterior_predictive_integrate(x_ic)
  mse = list_npb[['mse']] 
  auc = list_npb[['auc']] 
  lambda=list_npb[['lambda']] 
  param_est = list_npb[['param_est']]
  
  
  text1 = max(round(x_ic,2), signif(x_ic,3) )
  text2 = round(mse,2)
  text3 = round(auc)
  
  text = sprintf("atop(atop(IC[50] == %s, AUC == %s),atop( MSE == %s, \t) )",text1, text3, text2)
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  colls <<- c(colls, "NPB"="purple", "IC50NPB"="purple")
  linetypes <<- c(linetypes, "solid", "dotted")
  shapes <<- c(shapes, NA, NA)
  
  
  # We can add scalecolormanual since there is no other layer  to add on top
  p <- p +  
    
    geom_function(fun = posterior_predictive_integrate, aes(colour='NPB')) + 
    geom_hline( yintercept =  y_ic, color='purple',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic, colour="IC50NPB"),  linetype="dotted", show.legend = F) + 
    annotate(geom = 'text', y= y_lim_right, x =min(block2$doses), 
             hjust=-0.1,
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

# Generate plot for paper

# artificial input
input = list(mean_switch=T, outlier_switch=T, onehunda_switch=T, dosedep_auc=T,
             checkgroup1 = c("Point Samples" ,
                             "Spline" ,
                             "Min-Max Bands",
                             "Empirical Viability Bands",
                             "Drug Span Gradient",
                             # "Absolute doses" ,
                             "Relative Doses")
)

df = openxlsx::read.xlsx('data/Drug_response_S8.xlsx', sheet = 1)
df_list  = read_excel_allsheets('data/Drug_response_S8.xlsx')
df_example = read.csv('data/Example1.csv')

block = df_example
block2 = preprocess_data(block, mean_samples = input$mean_switch, keep_outliers = input$outlier_switch, over_viability = input$onehunda_switch)

p1 = PlotOverlay(block2, input$checkgroup1, input$dosedep_auc)
p1_mod = p1  + theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank(),
                     plot.title = element_blank())
  
p_t = plot_sigmodiMonotoneNpbFit(block2) +
  labs(title = NULL)



library(gridExtra)
p_paper = gridExtra::grid.arrange(rbind(ggplotGrob(p1_mod), ggplotGrob(p_t), size = "first"))

wid  = 6*4
hei = 4*4
ggsave('extras/ends.pdf',plot = p_paper, device = 'pdf',dpi = 400, width = wid, height = hei*2, unit = 'cm' )

# Seems that ic50 for NPB is not correct, change it and recompute statistics and make plot