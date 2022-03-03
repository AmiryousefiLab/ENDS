# Change directory to file save location
# library(rstudioapi)
# currentpath = rstudioapi::getActiveDocumentContext()$path
# setwd(dirname(currentpath))

library(ggplot2)
library(dplyr)
library(ggthemes)
library(fdrtool)
library(ggrepel)

### ------------------------------------------------------------------------------------------
# Implementing modeling ideas that Ali thought about


# We can add each of parts of the functions to the plots, so we have the implementation of building on top of each graph
# Change ylim from min to max of the samples
# p = p + new_function(p)



plot_initialize = function(block2, title=''){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  if(title=='') title = 'Nonparametric Spline'
  
  # Initialize global variables that will be filled in in legend
  colls <<- c()
  linetypes <<- c()
  shapes <<- c()
  alphas <<- c()
  
  p = block2 %>% 
    ggplot(aes(doses, y_mean) )+ 
    geom_point(size=2.5, stroke=0, shape=16, color='darkblue') +
    scale_x_log10(n.breaks=12, limits = c(min(block2[,1]), max(block2[,1]))) + 
    ylim(min(block2[,2:(m+1)], na.rm=T), y_lim_right) +
    ylab('Drug response') + 
    ggtitle( title ) +
    xlab(expression( paste(italic(Dose),phantom(x) ,mu, M )) ) +
    theme(plot.title = element_text(face="bold", size = 25, hjust=0.5),
          axis.title = element_text(face="italic", 
                                    color = "#267A43", 
                                    vjust=-0.35,  size=20),
          panel.grid = element_line(colour = "lightgray"),
          panel.border = element_rect(fill = NA, 
                                      colour = "black",
                                      size = rel(1)),
          axis.text.x = element_text(face="plain", 
                                     color="black", 
                                     size= 12, 
                                     angle = 40,
                                     vjust=0.8),
          axis.text.y = element_text(face="plain", 
                                     color="black", 
                                     size= 12),
          panel.background = element_rect(fill = "white",
                                          colour = "lightgray",
                                          size = 0.5, linetype = "solid"),
          legend.position = 'right',
          legend.text = element_text(size=15)
    ) 
  
  return(p)
}


nonparaametric_fit = function(block2, dose_dependent_auc=T, p_ic = 50){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  y_fit = block2$y_mean
  x_fit = block2$doses
  q_ic = 100-p_ic
  
  y_ic = min(y_fit)+ (q_ic/100)*(max(y_fit)-min(y_fit)) # y_ic value
  if(q_ic==0){
    x_ic = min(x_fit[y_fit==min(y_fit)])
  }
  if(q_ic==100){
    x_ic = min(x_fit[y_fit==max(y_fit)])
  }
  if(0<q_ic &  q_ic<100){
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
    xy_fit = ic_50_monotonefit(x_fit, y_fit, p_ic)
    x_ic_mono = xy_fit[[1]]
    y_ic_mono = xy_fit[[2]]
    
    if(is.nan(x_ic_mono)){ # Try again with flipped values, if monotonic fit failed
      mono1 = fdrtool::monoreg(x = log10(block2$doses), y = -block2$y_mean, type = 'antitonic')
      m = dim(block2)[2]-2
      if(m==0) m <- m+1
      y_fit = mono1$yf
      x_fit = (block2$doses)
      xy_fit = ic_50_monotonefit(x_fit, y_fit, p_ic)
      x_ic_mono = xy_fit[[1]]
      y_ic_mono = -xy_fit[[2]]
    }
    
    # ic50 is defined as the closest one to ic_50_mono
    x_ic = x_ic_multiple[which.min(abs(x_ic_multiple-x_ic_mono))]
    
  }
  
  
  if(dose_dependent_auc==TRUE) 
    auc = line_integral(block2$doses,block2$y_mean)
  
  if(dose_dependent_auc==FALSE) 
    auc = line_integral(1:length(block2$doses), block2$y_mean)
  
  
  samples = block2[2:(m+1)]
  y = block2$y_mean
  mse = sample_meansquarederror(y, samples)
  
  # MinMax Bands
  if(m==1){
    block2_y_max =  block2$y_mean
    block2_y_min =  block2$y_mean
  }
  if(m>1){
    block2_y_max = apply( block2[, 2:(dim(block2)[2]) ],1, function(x) max(x, na.rm = T) )
    block2_y_min = apply( block2[, 2:(dim(block2)[2]) ],1, function(x) min(x, na.rm = T) )  
  }
  
  # Absolute Doses
  
  # Drug Span Gradient
  x = log10(block2$doses)
  # divid eover 100 to standarize
  y = (block2$y_mean)/100
  
  model = lm(y~x)
  angle = atan(model$coefficients[2])
  angle_pi = round(angle/pi,2)
  angle_degrees = angle*180/pi
  
  # Absolute Doses
  angles = atan(diff(y)/diff((x)))
  angles_degrees = angles*180/pi
  
  list_stats = list( ic50 = x_ic, 
                     y_ic = y_ic,
                     mse = mse,
                     auc = auc,
                     drug_span_grad_angle = angle_degrees,
                     spline_angles = angles_degrees
  )
  
  return(list_stats)  
}


plot_NPDS = function(p, block2, dose_dependent_auc=TRUE, p_ic=50){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  list_nps = nonparaametric_fit(block2, dose_dependent_auc=T, p_ic )
  x_ic = list_nps$ic50
  y_ic = list_nps$y_ic
  mse = list_nps$mse
  auc = list_nps$auc
  angle_degrees = list_nps$drug_span_grad_angle
  angles_degrees = list_nps$spline_angles
  
  
  if(max(block2[,2:(m+1)], na.rm=T)==100){
    y_lim_right = 110
  }else{
    y_lim_right = max(block2[,2:(m+1)], na.rm=T)
  }
  
  # Delete last geom_text_repel from plot, so text is not overlaid
  if(length(p$layers)>0){
    flag=0
    for(i in 1:length(p$layers)){
      if(class(p$layer[[i]]$geom)[1] == 'GeomTextRepel'){
        j = i
        flag = 1
      }
    }
    if(flag==1){
      p$layers[[j]] = NULL  
    }
  }
  
  text0 = p_ic
  text1 = max(round(x_ic,2), signif(x_ic,3) )
  text2 = round(mse,2)
  text3 = round(auc)
  text = sprintf("atop(atop(IC[%s] == %s, AUC == %s),atop( MSE == %s, \t) )",text0,text1, text3, text2)
  
  colls <<- c(colls, "npS"="darkblue", "IC"="red")
  linetypes <<- c(linetypes, "solid", "dotted")
  shapes <<- c(shapes, NA, NA)
  alphas <<- c(alphas, 1, 1)
  p <-  p + 
    geom_line(aes(colour='npS') , size = 0.8) +
    geom_hline( yintercept =  y_ic, color='red',  linetype="dotted") +
    geom_vline(  aes(xintercept =  x_ic, colour="IC"),  linetype="dotted", show.legend = F) +
    ggrepel::geom_text_repel(aes(label=round(y_mean) ), size=3.0, force_pull = 2, seed=42) +
    annotate(geom = 'text', y= y_lim_right, x =max(block2$doses), 
             hjust=1,
             vjust=1,
             label = text, parse =T, size = 7, color='blue')
  
  
  return(p)
}


plot_point_samples =  function(p, block2){
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  colls <<- c(colls, "Samples"="blue")
  linetypes <<- c(linetypes, "blank")
  shapes <<- c(shapes, 16)
  alphas <<- c(alphas, 0.2)
  
  for(i in 2:(m + 1)){
    if(i ==2){
      # p = p + geom_point(aes(doses, !!sym(names(block2)[i])), size= 1.2, stroke=2, shape='\u2014', color='blue', alpha=0.2) 
      p = p + geom_point(aes(doses, !!sym(names(block2)[i]), colour = "Samples"), 
                         size= 1.5, 
                         alpha=0.3, 
                         shape=16)   
    }
    if(i>2){
      # p = p + geom_point(aes(doses, !!sym(names(block2)[i])), size= 1.2, stroke=2, shape='\u2014', color='blue', alpha=0.2) 
      p = p + geom_point(aes(doses, !!sym(names(block2)[i])), 
                         color='blue',
                         size= 1.5, 
                         alpha=0.3, 
                         shape=16)   
    }
    
  }
  return(p)
}


# Empirical variability band
plot_minmaxBands = function(p, block2){
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  if(m==1){
    block2_y_max =  block2$y_mean
    block2_y_min =  block2$y_mean
  }
  if(m>1){
    block2_y_max = apply( block2[, 2:(dim(block2)[2]) ],1, function(x) max(x, na.rm = T) )
    block2_y_min = apply( block2[, 2:(dim(block2)[2]) ],1, function(x) min(x, na.rm = T) )  
  }
  
  # Pretty ggplot visualization
  colls <<- c(colls, "MMB"="dimgrey")
  linetypes <<- c(linetypes, "solid")
  shapes <<- c(shapes, NA)
  alphas <<- c(alphas, 1)
  p = p + 
    geom_line(aes(doses, block2_y_min), color='dimgrey', linetype = 'dashed') +
    geom_line(aes(doses, block2_y_max, colour='MMB'), linetype = 'dashed')  
  
  return(p)
}

# Adjacent Mean Variability Band I dont understand the idea behind the variability bands for the adjacent mean lines
# I understand adjacent mean lines are the step function generated by the means of the drug response data


plot_empiricalVariabilityBand = function(p, block2){
  
  block2_doses_end  = lead(block2$doses)
  block2_adjMeanLine = (block2$y_mean + lead(block2$y_mean))/2
  # Pretty ggplot visualization
  colls <<- c(colls, "EVB"="darkgreen")
  linetypes <<- c(linetypes, "solid")
  shapes <<- c(shapes, NA)
  alphas <<- c(alphas, 1)
  p = p +
    geom_segment(aes(x=doses, y=block2_adjMeanLine, xend = block2_doses_end, yend=block2_adjMeanLine, color='darkgreen') ) +
    geom_line(aes(x=doses, y=block2_adjMeanLine) ,  color='darkgreen') +
    geom_line(aes(x = block2_doses_end, y=block2_adjMeanLine,  colour='EVB') )
  
  return(p)  
}

# Drug Span-Gradient

# we compute linear regression fit of means and then output slope/degree_angle of the fit. Same for max and min, generate plot and output the value
plot_drugSpanGradient = function(p, block2){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  x = log10(block2$doses)
  # divid eover 100 to standarize
  y = (block2$y_mean)/100
  
  model = lm(y~x)
  angle = atan(model$coefficients[2])
  angle_pi = round(angle/pi,2)
  angle_degrees = angle*180/pi
  # angle_degrees_plot = paste0(round(angle_degrees,1), intToUtf8(176))
  angle_degrees_plot = round(angle_degrees,1) 
  
  colls <<- c(colls, "DSG"="darkred")
  linetypes <<- c(linetypes, "solid")
  shapes <<- c(shapes, NA)
  alphas <<- c(alphas, 1)
  p = p + 
    geom_smooth(aes(colour = 'DSG'), method = "lm", se = F,  linetype = 'dashed', size=0.8, formula =  y ~ x) +
    # annotate('text',x = max(block2$doses)*0.8, y = max(block2[,2:(m+1)], na.rm=T)*0.8, size = 5, label =  bquote( theta~'='~ .(angle_degrees_plot) ) ) 
    annotate('text',x = min(block2$doses)*1.2, y = min(block2[,2:(m+1)], na.rm=T)*1.2,
             size = 5, 
             label =  bquote( theta~'='~ .(angle_degrees_plot)^o ),
             # label =  bquote( 'DSG='~ .(angle_degrees_plot)^o ),
             hjust=0) 
  return(p)  
}


# Most effective dose, we calculate the angle of each dose and return the one with highest negative slope
# Compare slopes with drug span-grad and label in positive or negative relative dose
# positive absolute doses, compared with line with slope 0
# concave ratio = #PRD/#NRD
# alive dose and dead dose are names for initial and final dose


# I propose constructing concave-convex ratio
# as the integral of function - mean line (how could we normalize this quantity)
# If positive then more concave than convex, also other way around
# Also the variance around the mean is the integral of the squared difference


plot_relativedoses = function(p, block2, relative = TRUE){
  
  m = dim(block2)[2]-2
  if(m==0) m <- m+1
  
  
  # For this funciton preferably call drug span gradient before
  # and NPDS as well
  
  x = log10(block2$doses)
  # divid eover 100 to standarize
  y = (block2$y_mean)/100
  
  model = lm(y~x)
  angle = atan(model$coefficients[2])
  angle_pi = round(angle/pi,2)
  angle_degrees = angle*180/pi
  
  angles = atan(diff(y)/diff((x)))
  angles_degrees = angles*180/pi
  
  
  
  # absolute doses
  abs_doses = 1*(angles_degrees<0)
  # relative doses, effective/positive if it is less than mean
  rel_doses = 1*(angles_degrees<angle_degrees)
  
  concave_ratio = round(sum(rel_doses)/(sum(rel_doses==0)),1)
  
  
  
  linetypes <<- c(linetypes, "blank")
  shapes <<- c(shapes, 16)
  alphas <<- c(alphas, 1)
  
  
  # Delete last geom_text_repel from plot, so text is not overlayed
  if(length(p$layers)>0){
    flag=0
    for(i in 1:length(p$layers)){
      if(class(p$layer[[i]]$geom)[1] == 'GeomTextRepel'){
        j = i
        flag = 1
      }
    }
    if(flag==1){
      p$layers[[j]] = NULL  
    }
  }
  
  if(relative == TRUE){
    colls <<- c(colls, "REED"="chocolate4")
    angles_degrees_plot = c(paste0(round(angles_degrees-angle_degrees,1), intToUtf8(176)) ,'')
    p = p +  
      ggrepel::geom_text_repel(aes(label = angles_degrees_plot, colour='RED'), show.legend = F, alpha=1) 
    #+  annotate('text',x = min(block2$doses)*1.2, y = min(block2[,2:(m+1)])*1.1, size = 5,label =  paste0('CR=', concave_ratio), hjust=0) 
  }
  if(relative == FALSE){
    colls <<- c(colls, "AED"="purple")
    angles_degrees_plot = c(paste0(round(angles_degrees,1), intToUtf8(176)) ,'')
    p = p +  
      ggrepel::geom_text_repel(aes(label = angles_degrees_plot, colour='AED'), show.legend = F, alpha=1) 
    #+  annotate('text',x = min(block2$doses)*1.2, y = min(block2[,2:(m+1)])*1.1, size = 5,label =  paste0('CR=', concave_ratio), hjust=0) 
  }
  
  
  
  
  return(p)  
}
