# This code creates a new R-Project, with readme, adds it to github, and creates folders
library(usethis)

# define project name
project <- "myproject"
path    <- "D:/OneDrive/Promotion/1 Projects/"

# Creates and switches you to project
create_project(paste0(path, project))


# Code for after the Project was created 
library(usethis)

# define project name
project <- "myproject"

# Add documentation
use_readme_md()

# Set up git
use_git()
use_github

# Create directories
fs::dir_create(c("1 Models", "2 Data", "3 Scripts", "4 Plots", "5 Results", "6 Experiment"))

# Add data directory and results to .gitignore
cat("2 Data",   file = ".gitignore", append = TRUE, sep = "\n")
cat("5 Results",file = ".gitignore", append = TRUE, sep = "\n")

# Create and open up a file
file.edit(paste0(project,"_README.md"))