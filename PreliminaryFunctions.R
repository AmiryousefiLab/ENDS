

# Preliminary Functions shared between scripts
# loaded before everything else


library(ggplot2)
library(openxlsx)
library(dplyr)
library(ggthemes)
library(grid)
library(gridExtra)
library(plyr)


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
    p <- plot_relativedoses(p, block2, relative = TRUE)
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


model_Statistics = function(block2, SplinePlot, SigmoidPlot, MonotonePlot, NPBPlot, dosedep_auc, p_ic, mean_switch, outlier_switch, onehunda_switch){
  
    l1=l2=l3=l4=NULL
    if(SplinePlot){
      list_nps = nonparaametric_fit(block2, dosedep_auc, p_ic)
      l1 = unlist(list_nps)
    }
    if(SigmoidPlot){
      list_pl = sigmoid_fit(block2, dosedep_auc, p_ic)
      list_pl$logistic_curve = NULL
      l2 = unlist(list_pl)
      
    }
    if(MonotonePlot){
      list_npm = monotone_fit(block2, dosedep_auc, p_ic)
      l3 = unlist(list_npm)
    }
    if(NPBPlot){
      block2 = preprocess_data_mult(block, mean_samples = mean_switch, keep_outliers = outlier_switch, over_viability = onehunda_switch, drop_values=F)
      list_npb = npb_fit(block2, dosedep_auc, p_ic)
      n = length(list_npb)
      list_npb[[n]] = NULL
      names(list_npb$param_est) = c('C_est', 'sigma2_est', 'a_est')
      l4 = unlist(list_npb)
    }
    # DataFrame to save results
    m = sum(!is.null(l1),!is.null(l2),!is.null(l3),!is.null(l4))
    n = max(length(l1), length(l2), length(l3), length(l4))
    if(n>0){
      
      M = matrix(NA, nrow=n, ncol=2*4)
      ic_name = paste('ic', p_ic, sep='')
      if(!is.null(l1)){
        M[1:length(l1),1] = c(ic_name, 'y_ic',"mse", "auc", "drug_span_grad_angle", "spline_angles1",
                              "spline_angles2", "spline_angles3", "spline_angles4", "spline_angles5",
                              "spline_angles6", "spline_angles7", "spline_angles8", "spline_angles9",
                              "spline_angles10", "spline_angles11", "spline_angles12", "spline_angles13",
                              "spline_angles14", "spline_angles15", "spline_angles16", "spline_angles17",
                              "spline_angles18", "spline_angles19", "spline_angles20")
        M[1:length(l1),2] = l1
      }
      if(!is.null(l2)){
        M[1:length(l2),3] = c(ic_name, 'y_ic', "mse", "auc", "coefficients.b:(Intercept)", "coefficients.c:(Intercept)",
                              "coefficients.d:(Intercept)", "coefficients.e:(Intercept)")
        M[1:length(l2),4] = l2
      }
      if(!is.null(l3)){
        M[1:length(l3),5] = c(ic_name, 'y_ic', "mse", "auc", "y_fit1", "y_fit2", "y_fit3", "y_fit4",
                              "y_fit5", "y_fit6", "y_fit7", "y_fit8", "y_fit9", "y_fit10",
                              "y_fit11", "y_fit12", "y_fit13", "y_fit14", "y_fit15", "y_fit16",
                              "y_fit17", "y_fit18", "y_fit19", "y_fit20", "y_fit21")
        M[1:length(l3),6] = l3
      }
      if(!is.null(l4)){
        M[1:length(l4),7] = c(ic_name, 'y_ic', "mse", "auc", "lambda", "C_est", "sigma2_est",
                              "a_est1", "a_est2", "a_est3", "a_est4",
                              "a_est5", "a_est6", "a_est7", "a_est8",
                              "a_est9", "a_est10", "a_est11",
                              "a_est12", "a_est13", "a_est14",
                              "a_est15", "a_est16", "a_est17",
                              "a_est18", "a_est19", "a_est20",
                              "a_est21")
        M[1:length(l4),8] = l4
      }
      df_stats = as.data.frame(M)
      colnames(df_stats)=c('npS names','npS values','pL names','pL values','npM names','npM values','npB names','npB values')
      return(df_stats)
    }
}


##############################################
# Functions for many one or more drugs-response bundles

delete_na_columns = function(df){
  df <- df[,colSums(is.na(df) | df=='' )<nrow(df)]
  return(df)
}

create_blocks <- function(tbl){
  drugs = tbl[,1] 
  drugs_unique = unique(drugs)
  drugs_n = length(drugs_unique)
  blocks = list()
  for(i in 1:drugs_n){
    t_d = tbl[drugs==drugs_unique[i],]
    t_d2 = delete_na_columns(t_d)
    blocks[[i]] = t_d2
  }
  blocks[[drugs_n+1]] = drugs_unique
  return(blocks)
}

preprocess_data_mult = function( block, mean_samples = TRUE, keep_outliers = TRUE, over_viability = TRUE, drop_values = TRUE){
  n = length(block)-1
  block2 = list()
  for(i in 1:n){
    block2[[i]] =  preprocess_data(block[[i]][,-1], mean_samples, keep_outliers , over_viability , drop_values )
  }
  block2[[n+1]] = block[[n+1]]
  return(block2)
}


PlotOverlay_mult = function(block2, check_boxes, dose_dependent_auc=TRUE, p_ic=50, title = ''){
  # if(title=='') title = 'npS'
  
  n = length(block2)-1
  drugs = block2[[n+1]]
  plots = list()
  for(i in 1:n){
    plots[[i]] =  PlotOverlay(block2[[i]], check_boxes, dose_dependent_auc=TRUE, p_ic=50, title = drugs[i])
  }
  # Create row of plots with given title 
  wid  = 6*4
  hei = 4*4 
  p = gridExtra::grid.arrange(grobs=plots, ncol=n, nrow=1, widths = rep(wid, n), heights=hei, top=textGrob(title, gp=gpar(fontsize=15,font=8)) )
  return(p)
}


blocks_to_csv = function(block2){
  n = length(block2)-1
  drugs = block2[[n+1]]
  df_csv = data.frame()
  for(i in 1:n){
    temp = cbind(drugs[i], block2[[n-i+1]])
    df_csv = rbind.fill(df_csv, temp)
  }
  names(df_csv)[1] = 'drugs'
  return(df_csv)
}


model_Statistics_mult = function(block2, SplinePlot, SigmoidPlot, MonotonePlot, NPBPlot, dosedep_auc, p_ic, mean_switch, outlier_switch, onehunda_switch){
  n = length(block2)-1
  drugs = block2[[n+1]]
  df_stats = data.frame()
  for(i in 1:n){
    df_temp = model_Statistics(block2[[i]], SplinePlot, SigmoidPlot, MonotonePlot, NPBPlot, dosedep_auc, p_ic, mean_switch, outlier_switch, onehunda_switch)
    df_temp = cbind(drugs[i], df_temp) # add column at beginning
    df_stats = rbind.fill(df_stats, df_temp)
  }
 return(df_stats) 
}
