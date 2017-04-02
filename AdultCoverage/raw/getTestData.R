# Author: tim
###############################################################################

Brasilfemales <- read.table("/home/tim/git/AdultCoverage/AdultCoverage/Data/UFdata_females_period2.txt",
		header = TRUE, stringsAsFactors = FALSE)
dim(Brasilfemales)
save(Brasilfemales,file="/home/tim/git/AdultCoverage/AdultCoverage/R/DDM/data/BrasilFemales.rda")

BrasilMales <- read.table("/home/tim/git/AdultCoverage/AdultCoverage/Data/UFdata_males_period1.txt",
		header = TRUE, stringsAsFactors = FALSE)
head(BrasilMales)
save(BrasilMales,file="/home/tim/git/AdultCoverage/AdultCoverage/R/DDM/data/BrasilMales.rda")

# documentation located in utils.R

