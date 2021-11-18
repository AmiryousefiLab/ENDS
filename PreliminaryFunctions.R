

# Preliminary Functions shared between scripts
# loaded before everything else


library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggthemes)


# ENDS (Epistemic nonparametric drug scoring)
# By blocks adjust linear piecewise function and save value, return plot with means and points

clean_doses = function(doses){
  
  idx_um = sapply(doses, function(x) grepl('uM',x))
  idx_nm = sapply(doses, function(x) grepl('nM',x))
  
  doses_numeric =  as.numeric(sapply(doses, function(x)( gsub("[^0-9.-]", "", x))))
  # why 58.6.2? typo or on purpose
  doses_numeric[is.na(doses_numeric)] = 58.6
  # Convert all to uM 
  doses_numeric[idx_nm] = doses_numeric[idx_nm]*1e-3
  return(doses_numeric)
}


piecewise_linear = function(x, block2){
  doses_numeric = block2['doses']
  rms = block2['y_mean']
  
  if(x<min(doses_numeric) | x>max(doses_numeric)) return(NA)
  # Minimum value such that x is greater than 
  x_i1 = min(doses_numeric[doses_numeric>=x])
  x_i0 = max(doses_numeric[doses_numeric<=x])
  
  if(x_i0==x_i1 ) return(x_i0)
  if(x_i0<x_i1){
    y_hat_i0 = as.numeric(rms[doses_numeric==x_i0])
    y_hat_i1 = as.numeric(rms[doses_numeric==x_i1])
    m = (y_hat_i1-y_hat_i0)/(x_i1-x_i0)
    value = m * (x-x_i0) +  y_hat_i0
    return(value)
  }
}



read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
# Devise function to extract data from the spredsheet


extract_dose_block = function(df_list, drug, patient, treatment, sample){
  sheets = names(df_list)
  # in sheets look for patient and drug
  # drug = '5FU'
  # patient = 'P1'
  # treatment = 'T1'
  # sample = 5
  block2=NULL
  try({
    idx = which(grepl(drug, sheets) & grepl(patient, sheets))
    df = df_list[[idx]]
    treatments = names(df)
    doses = df[,1]
    qry2 = paste(treatment, as.character(sample), sep='.')
    idt = which(grepl(qry2, treatments))
    # Always in blocks of 4
    block = df[, idt:(idt+3)]
    block2 = cbind(doses = clean_doses(doses), block)}
    ,T
  )
  return(block2)
}

preprocess_data = function(block, mean_samples = TRUE, keep_outliers = TRUE, over_viability = TRUE, drop_values = TRUE){
  m =  dim(block)[2]
  names(block)[1] = 'doses'
  block <- block %>% 
    arrange(doses)
  
  block2 = block
  
  # Further on can mark outliers on plot in a different color
  if(keep_outliers==FALSE & m>2 & drop_values==T){
    for(j in 2:m){
      y_toclean = block[,j]
      mu = mean(y_toclean, na.rm = T)
      stdev = sd(y_toclean, na.rm = T)
      y_clean = y_toclean
      y_clean[(y_toclean-mu)>2*stdev] <- NA
      y_clean[(y_toclean-mu)< -2*stdev] <- NA
      block2[,j] = y_clean
    }
  }
  
  if(keep_outliers==FALSE & m>2 & drop_values==F){
    for(j in 2:m){
      y_toclean = block[,j]
      mu = mean(y_toclean, na.rm = T)
      stdev = sd(y_toclean, na.rm = T)
      y_clean = y_toclean
      y_clean[(y_toclean-mu)>2*stdev] <- 2*stdev
      y_clean[(y_toclean-mu)< -2*stdev] <- -2*stdev
      block2[,j] = y_clean
    }
  }
  
  # For clipping viability
  # should we clip each of the samples or only average?
  if(over_viability == FALSE){
    block2[2:m][block2[2:m]>100] = 100
  }
  
  if(m==2){
    names(block2) = c('doses','y_mean')
  }
  if(m>2 & mean_samples == TRUE){
    block2['y_mean'] = apply(block2[,2:m], 1, mean, na.rm=T)
  }
  if(m>2 & mean_samples == FALSE){
    block2['y_mean'] = apply(block2[,2:m], 1, median, na.rm=T)
  }
  return(block2)
}


ic_50_monotonefit <- function(x_fit, y_fit, p_ic=50){
  
  # Correct from top to bottom
  q_ic = 100-p_ic
  y_ic = min(y_fit)+ (q_ic/100)*(max(y_fit)-min(y_fit))
  
  if(q_ic==0){
    x_ic = min(x_fit[y_fit==min(y_fit)])
  }
  if(q_ic==100){
    x_ic = min(x_fit[y_fit==max(y_fit)])
  }
  if(0<p_ic &  p_ic<100){
    # now the value on the xaxis for which we get that value on the yaxis
    idx1 = max(which(y_fit>=y_ic))
    idx2 = min(which(y_fit<=y_ic))
    
    if(idx1!=idx2){
      y1 = y_fit[idx1]
      y2 = y_fit[idx2]
      x1 = x_fit[idx1]
      x2 = x_fit[idx2]
      x_ic = (y_ic-y1)*((x2-x1)/(y2-y1)) + x1
    }else{
      print('IC_50 is not unique')
      x_ic = x_fit
    }
  }
  
  return(list(x_ic, y_ic))
}

# Build AUC function
line_integral <- function(x, y) {
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx *my)
} 

npds_auc = function(x, y){
  # x is doses
  # y is average drug resonse
  auc = line_integral(x,y)
  # auc = integrate(approxfun(x,y), range(x)[1], range(x)[2])
  return( auc )
}

# Mean squared error
sample_meansquarederror = function(y, samples){
  mean(as.matrix((samples-y)^2) , na.rm=T)
}

PlotOverlay = function(block2, check_boxes, dose_dependent_auc=TRUE, p_ic=50, title = ''){
  # controler for generated plot depending on checkboxes
  p = plot_initialize(block2, title )
  if(is.null(check_boxes) ) return(p)
  
  if("Point Samples" %in% check_boxes)
    p <- plot_point_samples(p, block2)
  if("Spline" %in% check_boxes)
    p <- plot_NPDS(p, block2, dose_dependent_auc, p_ic)
  if("Min-Max Bands" %in% check_boxes)
    p <- plot_minmaxBands(p, block2)
  if("Empirical Viability Bands" %in% check_boxes)
    p <- plot_empiricalVariabilityBand(p, block2)
  if("Drug Span Gradient" %in% check_boxes)
    p <- plot_drugSpanGradient(p, block2)
  if("Absolute Doses" %in% check_boxes)
    p <- plot_relativedoses(p, block2, relative=FALSE)
  if("Relative Doses" %in% check_boxes)
    p <- plot_relativedoses(p, block2)
  p <- p + 
    scale_colour_manual(name="Labels",values=colls,
                        guide = guide_legend(
                          override.aes =
                            list(
                              linetype = linetypes,
                              shape = shapes,
                              alpha = alphas
                            ))) 
  return(p)
}
