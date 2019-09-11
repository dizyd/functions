


# Function: Check if packages are installed, if not install and load them
#'   Input: @packages character vector specifying the packages to install and/or load
#   Output: Require output indicated if packages were loaded succesfully (TRUE) or not (FALSE)

install_n_load_packages         <- function(packages=c("tidyverse")){
  
  ## Make list of wanted packages 
  
  wanted.packages <- packages
  
  ## compare installed and wanted packages
  
  new.packages <- wanted.packages[!(wanted.packages %in% installed.packages()[,"Package"])]
  
  # installed the not yet installed but required packages and load them
  
  if(length(new.packages)) install.packages(new.packages,dependencies = TRUE)
  sapply(wanted.packages, require, character.only = TRUE)
  
}


# Function: Save Results with current date
#'   Input: @results results object to save
#'          @path    path were to save the object
#'          @name    file name
#   Output: Require output indicated if packages were loaded succesfully (TRUE) or not (FALSE)

save_res <- function(results,path,name){
  
  date <- format(Sys.time(), "%d%m%Y")
  save(results,file=paste0(path,name,"_",date,".Rdata"))
  
}