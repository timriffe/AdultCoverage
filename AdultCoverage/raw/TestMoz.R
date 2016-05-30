setwd("/home/tim/git/AdultCoverage/AdultCoverage")
x <- read.csv("Data/Mozambique.csv", stringsAsFactors = FALSE)

source("R/DDM/R/ggb.R")
source("R/DDM/R/bh.R")
source("R/DDM/R/utils.R")

ggb(x)
bh1(x,exact.ages=seq(15,55,by=5))
plot.ggb(tab1[[2]])$coverage


