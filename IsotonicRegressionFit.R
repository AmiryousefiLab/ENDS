
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

ic_50_monotonefit <- function(x_fit, y_fit){
  y_ic = 0.5*(max(y_fit)+min(y_fit))
  # now the value on the xaxis for which we get that value on the yaxis
  idx1 = max(which(y_fit>=y_ic))
  idx2 = min(which(y_fit<=y_ic))
  
  if(idx1!=idx2){
    y1 = y_fit[idx1]
    y2 = y_fit[idx2]
    x1 = x_fit[idx1]
    x2 = x_fit[idx2]
    x_ic = (y_ic-y1)*((x2-x1)/(y2-y1)) + x1
  } else{
    print('IC_50 is not unique')
    x_ic = x_fit
  }
  return(list(x_ic, y_ic))
}

plot_monotoneFit = function( block2, dose_dependent_auc=TRUE){
  
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
  
  text1 = max(round(x_ic,2), signif(x_ic,3) )
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
  
  colls <<- c(colls, "NP Monotone"="darkgreen","IC50"="red")
  linetypes <<- c(linetypes, "solid","dotted")
  shapes <<- c(shapes, NA, NA)
  
  options(warn=-1)
  p = p +
    ggtitle('Nonparametric Monotone') +
    geom_line(aes(block2$doses, block2_yf, colour ='NP Monotone')) + 
    geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic, colour="IC50" ),  linetype="dotted", show.legend = F) + 
    annotate(geom = 'text', y= y_lim_right, x =max(block2$doses), 
             hjust=1,
             vjust=1,
             label = text, parse =T, size = 7,
             color = 'darkgreen')+
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

# p = plot_monotoneFit(block2)
