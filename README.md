# functions

<!-- badges: start -->
<!-- badges: end -->

This repository contains a collection of functions and code snippets which I regularly use in my work. Instead over copying & pasting the functions every time again, I decided to collect these helper and convinience functions in this repository. 


Until now, this repository contains the following files:


## `helper_functions.R`

This is the main file which contains several useful functions, for instance to calculate the RMSE or the mode, flatten a nested list as outputted when using `foreach()`, or for reporting mean and sd in one function call. 

## `init_new_project.R` 

This file contains the code to create a new project, the necessary folder structure, the github repository and the readme file. To do this, it uses the `usethis` package.

## `make_readme_fun.R`

This file contains the `R`-code of the `make_df_readme()` function. This function returns a description of the variables contained in a data.frame in a `.txt` function. This file is then a usefull addition when you share your data files with your colleagues or put them on OSF (or just for yourself). The table returned is in Markdown format and thus can be easily copied into a github readme file or to the wiki on your OSF project.
