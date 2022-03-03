
# Sigmoid Fit (Parametric Adjustment)
# By blocks adjust linear piecewise function and save value, return plot with means and points
# Sigmoid fit is modified to be over the means of the data, return also mean square error. 

library(dplyr)
library(ggplot2)
library(nplr)
library(drc)

ic_50_sigmoidfit <- function(x, m0, p_ic=50){
  
  q_ic = 100-p_ic
  c = m0$coefficients['c:(Intercept)']
  d = m0$coefficients['d:(Intercept)']
  function_ajusted = m0$curve[[1]]
  
  if(q_ic==50){
    y_ic = (c+(d-c)*(q_ic/100))
    x_ic = m0$coefficients[4] 
  }
  if(q_ic<=1){
    y_ic = (c+(d-c)*(1/100))
    function_inverse = function(x){ function_ajusted(x) - y_ic }
    uniroot = uniroot(function_inverse,interval=c(min(x),max(x)), tol=1e-9)
    x_ic = uniroot$root
  }
  if(q_ic>=99){
    y_ic = (c+(d-c)*(99/100))
    function_inverse = function(x){ function_ajusted(x) - y_ic }
    uniroot = uniroot(function_inverse,interval=c(min(x),max(x)), tol=1e-9)
    x_ic = uniroot$root
  }
  if(q_ic>1 & q_ic<99 & (q_ic!=50) ){
    y_ic = (c+(d-c)*(q_ic/100))
    function_inverse = function(x){ function_ajusted(x) - (c+(d-c)*(q_ic/100)) }
    uniroot = uniroot(function_inverse,interval=c(min(x),max(x)), tol=1e-9)
    x_ic = uniroot$root
    # Add a clause for errors
  }
  
  
  return(list(x_ic, y_ic))
}


sigmoid_fit = function(block2, dose_dependent_auc=TRUE, p_ic = 50, viability_switch=TRUE){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  # Fit on means
  y = unlist(block2['y_mean'])
  x = unlist(block2['doses'])
  
  
  # q_ic = p_ic
  # Equation of the fit
  # f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f} 
  
  if(viability_switch==T){
    m0 <- drm(y ~ x,  fct = LL.4()) # 4 parameter logistic fit
    xy_fit = ic_50_sigmoidfit(x, m0, p_ic)
    x_ic = xy_fit[[1]]
    y_ic = xy_fit[[2]]  
    function_ajusted = m0$curve[[1]]
  }
  if(viability_switch==F){
    m0 <- drm(-y ~ x,  fct = LL.4()) # 4 parameter logistic fit
    xy_fit = ic_50_sigmoidfit(x, m0, p_ic)
    x_ic = xy_fit[[1]]
    y_ic = -xy_fit[[2]]  
    function_ajusted = function(x) -m0$curve[[1]](x)
  }
  
  samples = block2[2:(m+1)]
  mse = sample_meansquarederror(function_ajusted(x), samples) # Mean squared error for samples is
  
  
  if(dose_dependent_auc==TRUE)
    auc = integrate(function_ajusted, lower = min(x), max(x))$value
  if(dose_dependent_auc==FALSE){
    if(viability_switch==TRUE){
      x_ = 1:length(x)
      m1 = drm(y ~ x_,  fct = LL.4())
      auc = integrate(m1$curve[[1]], lower = min(x), max(x))$value  
    }
    if(viability_switch==FALSE){
      x_ = 1:length(x)
      m1 = drm(-y ~ x_,  fct = LL.4())
      auc = integrate(-m1$curve[[1]], lower = min(x), max(x))$value
    }
    
  }
  
  list_stats = list( ic50 = x_ic[[1]], 
                     y_ic = y_ic,
                     mse = mse,
                     auc = auc,
                     coefficients = m0$coefficients,
                     logistic_curve = function_ajusted
  )
  
  return(list_stats)  
}



plot_sigmodiFit = function(block2, dose_dependent_auc=TRUE, p_ic = 50, title = '', viability_switch=TRUE){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1

  list_pl = sigmoid_fit(block2, dose_dependent_auc, p_ic, viability_switch)
  x_ic = list_pl$ic50
  y_ic = list_pl$y_ic
  mse = list_pl$mse
  auc = list_pl$auc
  logistic_curve = list_pl$logistic_curve
  
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
  
  if(title=='') title = 'Parametric Logistic'
  
  p = plot_initialize(block2)
  p = plot_point_samples(p, block2)
  
  colls <<- c(colls, "pL"="red", "IC"="red")
  linetypes <<- c(linetypes, "solid", "dotted")
  shapes <<- c(shapes, NA, NA)
  
  if(viability_switch==TRUE){
    p <- p +  
      ggtitle(title) +
      geom_function(fun = logistic_curve, aes(colour='pL')) + 
      geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
      geom_vline(  aes(xintercept =  x_ic, colour="IC"),  linetype="dotted", show.legend = F) + 
      annotate(geom = 'text', y= y_lim_right, x =max(block2$doses), 
               hjust=1,
               vjust=1,
               label = text, parse =T, size = 7, 
               color='red') +
      scale_colour_manual(name="Labels",values=colls,
                          guide = guide_legend(
                            override.aes =
                              list(
                                linetype = linetypes,
                                shape = shapes
                              )
                          )
      )
  }
  if(viability_switch==FALSE){
    p <- p +  
      ggtitle(title) +
      geom_function(fun = logistic_curve, aes(colour='pL')) + 
      geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
      geom_vline(  aes(xintercept =  x_ic, colour="IC"),  linetype="dotted", show.legend = F) + 
      annotate(geom = 'text', y= y_lim_right, x = min(block2$doses), 
               hjust=0,
               vjust=1,
               label = text, parse =T, size = 7, 
               color='red') +
      scale_colour_manual(name="Labels",values=colls,
                          guide = guide_legend(
                            override.aes =
                              list(
                                linetype = linetypes,
                                shape = shapes
                              )
                          )
      )
  }
  
  # We can add scalecolormanual since there is no other layer  to add on top
  
  
  return(p)
}

###########################
# Multiple input functions
plot_sigmodiFit_mult = function( block2, dose_dependent_auc=TRUE, p_ic=50, title = '', viability_switch=TRUE){
  # if(title=='') title = 'pL'
  n = length(block2)-1
  drugs = block2[[n+1]]
  plots = list()
  for(i in 1:n){
    plots[[i]] =  plot_sigmodiFit(block2[[i]], dose_dependent_auc=TRUE, p_ic=50, title = drugs[i], viability_switch)
  }
  # Create row of plots with given title 
  wid  = 6*4
  hei = 4*4 
  p = gridExtra::grid.arrange(grobs=plots, ncol=n, nrow=1, widths = rep(wid, n), heights=hei, top=textGrob(title, gp=gpar(fontsize=15,font=8)) )
  return(p)
}

