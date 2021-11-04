# Get requests from drugcomb to build into webpage

library(httr)
library(jsonlite)
library(dplyr)

# GET /response/{block_id}
# GET /response?from={offset}&to={limit}

get_df_url <- function(url){
  r <- GET(url)
  # Each blockid corresponds to a different combination
  # in order to find a single drug we could search for a specific drug and characteristics in data to pinpoint it
  raise = content(r, as='text', encoding = 'UTF-8')
  new = fromJSON(raise)
  return(new)
}


drug_comb_info_retreival <-  function(block_id, drug_row=T){
  # I can trace back the drug that have this specific combination65432453454332
  url = paste0("https://api.drugcomb.org/response/", block_id)
  df_comb = get_df_url(url)
  
  # From summary get drug id and cell line id
  url = paste0("https://api.drugcomb.org/summary/", block_id)
  df_sum = get_df_url(url)
  
  if(drug_row){
    df_single = df_comb %>% 
      filter(conc_c==0)
    # with drug id search for names of drugs
    url = paste0("https://api.drugcomb.org/drugs/", df_sum$drug_row_id)
    df_drug = get_df_url(url)  
  } else{
    df_single = df_comb %>% 
      filter(conc_c==0)
    url = paste0("https://api.drugcomb.org/drugs/", df_sum$drug_col_id)
    df_drug = get_df_url(url)  
  }
  
  # with cell line id search for cell line name
  url = paste0("https://api.drugcomb.org/cell_lines/", df_sum$cell_line_id)
  df_cell_line = get_df_url(url)
  
  df_single$drug = df_drug$dname
  df_single$cell_line = df_cell_line$name
  
  df_return = df_single %>% 
    select(drug, cell_line,conc_r, conc_c, inhibition )
  
  return(df_return)
}

# drug_comb_info_retreival(block_id = 1)
# drug_comb_info_retreival(block_id = 2)
# drug_comb_info_retreival(block_id = 80)

# This last 2 are the same, but differ by minimal margin from previous
# drug_comb_info_retreival(block_id = 90949)
# drug_comb_info_retreival(block_id = 90952)

# We conclude that we cannot generate different samples per cell_line from drug_comb database
# We would like to have different samples per cell_line, check if drug_comb offers this (as with prev database)

  # Get all blocks with specific cell line and drug
  
# library(data.table)
# df_sum = fread('data/summary_v_1_5.csv')  
# df_look = df_sum %>% 
#   filter(drug_row =="5-Fluorouracil" & cell_line_name == "A2058")

# Are these 2 the same?
 # Information is not duplicated, not the opposite drug combination exists
# df_look1 = df_sum %>% 
#   filter(drug_row =="5-Fluorouracil" & drug_col == "Veliparib")
# df_look2 = df_sum %>% 
#   filter(drug_row =="Veliparib" & drug_col == "5-Fluorouracil")


