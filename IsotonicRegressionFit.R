
### Isotonic Regression Spline Fit
# By blocks adjust linear monotone decreasing function and save value, return plot with fit and full set of points


library(dplyr)
library(ggplot2)
library(fdrtool)



# plot_monotoneFit = function(p, block2){
#   
#   mono1 = fdrtool::monoreg(x = log10(block2$doses), y = block2$y_mean, type = 'antitonic')
#   
#   # Gives exatly same fit
#   # mono1 = fdrtool::monoreg(x = block2$doses, y = block2$y_mean, type = 'antitonic')
#   
#   # Only fits on the means, not on bulks of data
#   # mono2 = fdrtool::monoreg(x = block2$doses, y = block2[,1:4], type = 'antitonic')
#   # plot(mono1)
#   # Calculate IC_{50}
#   # It would be the middle point of the fit
#   
#   
#   y_fit = mono1$yf
#   x_fit = (block2$doses)
#   y_ic = 0.5*(max(y_fit)+min(y_fit))
#   # now the value on the xaxis for which we get that value on the yaxis
#   idx1 = max(which(y_fit>=y_ic))
#   idx2 = min(which(y_fit<=y_ic))
#   
#   if(idx1!=idx2){
#     y1 = y_fit[idx1]
#     y2 = y_fit[idx2]
#     x1 = x_fit[idx1]
#     x2 = x_fit[idx2]
#     x_ic = (y_ic-y1)*((x2-x1)/(y2-y1)) + x1
#   } else{
#     print('IC_50 is not unique')
#     x_ic = x_fit
#   }
#   
#   block2_yf = mono1$yf
#   
#   sqrd_error = sum( (block2$y_mean -  mono1$yf)^2)
#   
#   
#   text1 = max(round(x_ic,2), signif(x_ic,3) )
#   text2 = round(sqrd_error,2)
#   
#   text = sprintf("atop(IC[50] == %s, MSE == %s)",text1, text2 )
#   
#   p = p +
#     geom_line(aes(block2$doses, block2_yf), color='blue') + 
#     geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
#     geom_vline( xintercept =  x_ic, color='red',  linetype="dotted") + 
#     annotate(geom = 'text', y= max(block2$y_mean), x =max(block2$doses),hjust=1, label = text, parse=T) 
#   return(p)
# }
# 


monotone_fit = function(block2, dose_dependent_auc=TRUE, p_ic = 50, viability_switch=TRUE){
  
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  if(viability_switch==TRUE){
    mono1 = fdrtool::monoreg(x = log10(block2$doses), y = block2$y_mean, type = 'antitonic')
    y_fit = mono1$yf
    x_fit = (block2$doses)
    xy_fit = ic_50_monotonefit(x_fit, y_fit, p_ic)
    x_ic = xy_fit[[1]]
    y_ic = xy_fit[[2]]  
    block2_yf = mono1$yf
  }
  if(viability_switch==FALSE){
    mono1 = fdrtool::monoreg(x = log10(block2$doses), y = -block2$y_mean, type = 'antitonic')
    y_fit = mono1$yf
    x_fit = (block2$doses)
    xy_fit = ic_50_monotonefit(x_fit, y_fit, p_ic)
    x_ic = xy_fit[[1]]
    y_ic = -xy_fit[[2]]  
    block2_yf = -mono1$yf
  }
  

  
  # Mean squared error for samples is
  samples = block2[2:(m+1)]
  mse = sample_meansquarederror(block2_yf, samples)
  
  
  if(dose_dependent_auc==TRUE) 
    auc = line_integral(x_fit,block2_yf)
  
  if(dose_dependent_auc==FALSE) 
    auc =  line_integral(1:length(x_fit), block2_yf)
  
  
  
  list_stats = list( ic50 = x_ic, 
                     y_ic = y_ic,
                     mse = mse,
                     auc = auc,
                     y_fit = block2_yf
  )
  
  
  return(list_stats)  
}

plot_monotoneFit = function( block2, dose_dependent_auc=TRUE, p_ic=50, title = '', viability_switch=TRUE, stat_info=T, x_ticks=T){
  # Gives exatly same fit
  # mono1 = fdrtool::monoreg(x = block2$doses, y = block2$y_mean, type = 'antitonic')
  
  # Only fits on the means, not on bulks of data
  # mono2 = fdrtool::monoreg(x = block2$doses, y = block2[,1:4], type = 'antitonic')
  # plot(mono1)
  # Calculate IC_{50}
  # It would be the middle point of the fit
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  list_npM = monotone_fit(block2, dose_dependent_auc, p_ic, viability_switch)
  x_ic = list_npM$ic50
  y_ic = list_npM$y_ic
  mse = list_npM$mse
  auc = list_npM$auc
  y_fit = list_npM$y_fit
  
  if(title=='') title = 'Nonparametric Monotonic'
  
  text0 = p_ic
  text1 = max(round(x_ic,2), signif(x_ic,3) )
  text2 = round(mse,2)
  text3 = round(auc)
  text = sprintf("atop(atop(IC[%s] == %s, AUC == %s),atop( MSE == %s, \t) )",text0,text1, text3, text2)
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  p = plot_initialize(block2, x_ticks=x_ticks)
  p = plot_point_samples(p, block2)
  
  colls <<- c(colls, "npM"="darkgreen","IC"="red")
  linetypes <<- c(linetypes, "solid","dotted")
  shapes <<- c(shapes, NA, NA)
  
  options(warn=-1)
  p = p +
    ggtitle( title ) +
    geom_line(aes(block2$doses, y_fit, colour ='npM')) + 
    geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic, colour="IC" ),  linetype="dotted", show.legend = F)
  
  if(viability_switch==TRUE & stat_info==T){
     p = p + annotate(geom = 'text', y= y_lim_right, x =max(block2$doses), 
               hjust=1,
               vjust=1,
               label = text, parse =T, size = 7,
               color = 'darkgreen')
  }
  if(viability_switch==FALSE & stat_info==T){
    p = p + annotate(geom = 'text', y= y_lim_right, x =min(block2$doses), 
               hjust=0,
               vjust=1,
               label = text, parse =T, size = 7,
               color = 'darkgreen')
      
  }
  p <- p + scale_colour_manual(name="Labels",values=colls,
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

###########################
# Multiple input functions
plot_monotoneFit_mult = function( block2, dose_dependent_auc=TRUE, p_ic=50, title = '', viability_switch=TRUE, stat_info=T, x_ticks=T){
  # if(title=='') title = 'npM'
  n = length(block2)-1
  drugs = block2[[n+1]]
  plots = list()
  for(i in 1:n){
    plots[[i]] =  plot_monotoneFit(block2[[i]],  dose_dependent_auc=TRUE, p_ic=50, title = drugs[i], viability_switch, stat_info, x_ticks)
  }
  # Create row of plots with given title 
  wid  = 6*4
  hei = 4*4 
  p = gridExtra::grid.arrange(grobs=plots, ncol=n, nrow=1, widths = rep(wid, n), heights=hei, top=textGrob(title, gp=gpar(fontsize=15,font=8)) )
  return(p)
}

