
# Sigmoid Fit (Parametric Adjustment)
# By blocks adjust linear piecewise function and save value, return plot with means and points
# Sigmoid fit is modified to be over the means of the data, return also mean square error. 

library(dplyr)
library(ggplot2)
library(nplr)
library(drc)

sigmoid_fit = function(block2, dose_dependent_auc=TRUE, p_ic = 50){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  # Fit on means
  y = unlist(block2['y_mean'])
  x = unlist(block2['doses'])
  
  q_ic = 100-p_ic
  # q_ic = p_ic
  
  # 4 parameter logistic fit,
  m0 <- drm(y ~ x,  fct = LL.4())
  
  
  # Mean squared error for samples is
  samples = block2[2:(m+1)]
  mse = sample_meansquarederror(m0$curve[[1]](x), samples)
  
  # Equation of the fit
  # f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f} 
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
  
  if(dose_dependent_auc==TRUE)
    auc = integrate(m0$curve[[1]], lower = min(x), max(x))$value
  if(dose_dependent_auc==FALSE){
    x_ = 1:length(x)
    m1 = drm(y ~ x_,  fct = LL.4())
    auc = integrate(m1$curve[[1]], lower = min(x), max(x))$value
  }
  
  list_stats = list( ic50 = x_ic[[1]], 
                     y_ic = y_ic,
                     mse = mse,
                     auc = auc,
                     coefficients = m0$coefficients,
                     logistic_curve = m0$curve[[1]]
  )
  
  return(list_stats)  
}



plot_sigmodiFit = function(block2, dose_dependent_auc=TRUE, p_ic = 50, title = ''){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  # Fit on means
  y = unlist(block2['y_mean'])
  x = unlist(block2['doses'])
  
  
  list_pl = sigmoid_fit(block2, dose_dependent_auc=TRUE, p_ic = p_ic)
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
  
  
  # We can add scalecolormanual since there is no other layer  to add on top
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
  
  return(p)
}

# p = plot_sigmodiFit(block2, T,30)
# p + theme(legend.box=)

# Use package nplr
# https://cran.r-project.org/web/packages/nplr/vignettes/nplr.pdf
# Fit on averages
# x = block2$doses
# y = block2$y_mean
# 
# np1 = nplr::nplr( x, y, useLog = TRUE)
# # np1 = nplr::nplr( x, convertToProp(y), useLog = TRUE)
# plot(np1, showSDerr = TRUE, lwd = 4 , cex.main=1.25, main="Cell line MCF-7. Response to Irinotecan")
# 
# 
#  # Fit on all available data
# y = unlist(block2[, 1:4])
# x = unlist(rep(block2['doses'],4))
# np2 = nplr::nplr( x=x, y=y,useLog=TRUE)
# plot(np2, showSDerr = TRUE, lwd = 4 , cex.main=1.25, main="Cell line MCF-7. Response to Irinotecan")
# 
# # Use drc traditional package
# y = unlist(block2[, 1:4])
# x = unlist(rep(block2['doses'],4))
# 
# # 4 parameter logistic fit,
# m0 <- drm(y ~ x,  fct = LL.4())
# doses = block2$doses
# plot(doses, m0$curve[[1]](doses), log = 'x', type='l')


