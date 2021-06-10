#devtools::source_url("https://github.com/dizyd/functions/blob/master/helper_functions.R?raw=TRUE")


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

# Function: Transform sigma to precision and reverse
#'   Input: @sigma sigma
#'          @tau   tau
#   Output: either sigma or tau, depending on which input was provided

sigma_tau <- function(tau = NA,sigma = NA){
  
  if(is.na(tau)){
    tau= 1/(sigma^2)
    return(tau)
  }
  
  if(is.na(sigma)){
    sigma = 1/sqrt(tau)
    return(sigma)
  }
  
}

# Function: flatten a (nested) list to a data.frame
#'   Input: @x (nested) list
#   Output: Data.frame

flattenlist     <- function(x){  
  more_lists <- sapply(x, function(first_x) class(first_x)[1]=="list")
  out        <- c(x[!more_lists], unlist(x[more_lists], recursive=FALSE))
  if(sum(more_lists)){ 
    Recall(out)
  }else{
    return(out)
  }
}


# Function: prepare a Readme overview for a tidy data.frame
#'   Input: @df   tidy data.frame
#'          @desc variabel descriptions
#'          @info additional info you might want to add before the data.frame, prints also the dimensions of the data.frame
#'          @file name of the output file
#   Output: saves .txt file in working directory

make_df_readme     <- function(df,desc,info = NULL,file = "readme.txt"){  
  
  temp0 <- data.frame("Variable"    = names(df),
                      "Description" = desc)
  
  # Start writing to the file
  sink(file)
  if(length(info)){
    cat(info,"\n","\n","data.frame (",nrow(df),",",ncol(df),")","\n","\n",sep = "")
  }
  
  temp1 <- temp0 %>% knitr::kable(format = "markdown") 
  paste0(temp1,"\n") %>% cat(sep = "")
  
  # Stop writing to the file
  sink()
  
}