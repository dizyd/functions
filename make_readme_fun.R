# Load with:
# devtools::source_url("https://github.com/dizyd/functions/blob/master/make_readme_fun.R?raw=TRUE")

# Function: prepare a Readme overview for a tidy data.frame
#'   Input: @df   data.frame
#'          @desc variable descriptions
#'          @info additional info you might want to add before the data.frame, prints also the dimensions of the data.frame
#'          @file name of the output file
#'          @add_examples if examples of data in the variables should be added to the readme
#   Output: saves .txt file in working directory

make_df_readme     <- function(df,desc,info = NULL,file = "readme.txt",add_examples=TRUE){  
  
  
  temp0 <- data.frame("Variable"    = names(df),
                      "Type"        = sapply(df, class),
                      "Description" = desc)
  
  row.names(temp0) <- NULL
  
  if(add_examples){
    
    temp_info <- df[sample(1:nrow(df),2),] %>%
      t() %>%
      as.data.frame() %>%
      apply(., 1, paste, collapse=",") %>% 
      unlist()
    
    names(temp_info) <- NULL
    temp0 <- temp0 %>% add_column("example" = temp_info,.before = "Description")
  }
  
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




# Example:

# df = data.frame("ID" = 1:50,
#                 "cond" = sample(c("A","B"),50,replace=T),
#                 "rt" = rnorm(50,500,40))
# 
# desc = c("unique numeric participant ID",
#          "condition [A,B]",
#          "average reaction time in ms")
# 
# make_df_readme(df = df,
#                desc = desc,
#                info = "Data of the experiment Hello World",
#                add_examples=TRUE)
