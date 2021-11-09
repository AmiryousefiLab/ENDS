
# Sigmoid Fit (Parametric Adjustment)
# By blocks adjust linear piecewise function and save value, return plot with means and points
# Sigmoid fit is modified to be over the means of the data, return also mean square error. 

library(dplyr)
library(ggplot2)
library(nplr)
library(drc)


plot_sigmodiFit = function(block2, dose_dependent_auc=TRUE){
  
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
  if(dose_dependent_auc)
    auc = integrate(m0$curve[[1]], lower = min(x), max(x))$value
  if(dose_dependent_auc==FALSE){
    x_ = 1:length(x)
    m1 = drm(y ~ x_,  fct = LL.4())
    auc = integrate(m1$curve[[1]], lower = min(x), max(x))$value
  }
    
  
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
  
  colls <<- c(colls, "P Logistic"="red", "IC50"="red")
  linetypes <<- c(linetypes, "solid", "dotted")
  shapes <<- c(shapes, NA, NA)
  
  
  # We can add scalecolormanual since there is no other layer  to add on top
  p <- p +  
    ggtitle('Parametric Logistic') +
    geom_function(fun = m0$curve[[1]], aes(colour='P Logistic')) + 
    geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic, colour="IC50"),  linetype="dotted", show.legend = F) + 
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

# p = plot_sigmodiFit(block2)



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


