#! /usr/local/bin/Rscript --vanilla --default-packages=utils

getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

# Set the script path:
script_path<- getScriptPath()
setwd(script_path)

#Load packages
package_list<-c("ape")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

system('python --version')

system('conda init bash')

system('conda activate dlcpar_py27')

system('python --version')

system('conda deactivate')

system('python --version')