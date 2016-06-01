setwd("/home/tim/git/AdultCoverage/AdultCoverage")

x <- read.csv("Data/Mozambique.csv", stringsAsFactors = FALSE)

source("R/DDM/R/cdmltw.R")
source("R/DDM/R/ggb.R")
source("R/DDM/R/bh.R")
source("R/DDM/R/utils.R")

ggb(x, exact.ages = seq(15,50,by=5))
ggb(x, exact.ages = seq(30,75,by=5))

bh1(x,exact.ages=seq(15,55,by=5))



#		
#library(devtools)
#install_github("timriffe/AdultCoverage/AdultCoverage/R/DDM")
#


ggb(bla, bla, bla, by = list("Sex", "Period"))


