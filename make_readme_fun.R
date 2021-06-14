# Load with:
# devtools::source_url("https://github.com/dizyd/functions/blob/master/make_readme_fun.R?raw=TRUE")

# Function: prepare a Readme overview for a tidy data.frame
#'   Input: @df   data.frame
#'          @desc variable descriptions
#'          @info additional info you might want to add before the data.frame, prints also the dimensions of the data.frame
#'          @file name of the output file
#   Output: saves .txt file in working directory

make_df_readme     <- function(df,desc,info = NULL,file = "readme.txt"){  
  
  
  temp0 <- data.frame("Variable"    = names(df),
                      "Type"        = sapply(df, class),
                      "Description" = desc)
  
  row.names(temp0) <- NULL
  
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