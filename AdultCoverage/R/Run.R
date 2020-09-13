
#ITA <- read.table("Data/italy.test.txt", 
#  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#
##Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)
#
## new data for Brasil to run. By regions.
#BR1 <- read.table(file.path("Data","data_Brazil_p1.txt"), 
#		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#BR2 <- read.table(file.path("Data","data_Brazil_p2.txt"), 
#		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#BR3 <- read.table(file.path("Data","data_Brazil_p3.txt"), 
#		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#
library(devtools)

setwd("AdultCoverage/R/DDM")

load_all()
document()
check()

#--------------------
library(dplyr)
library(magrittr)

# change data
ggb(Moz)
seg(Moz)
ggbseg(Moz)
res <- ddm(Moz)

codi <- Moz %>% ggbMakeColumns() 
ggbcoverageFromYear(codi,lm.method = "ols")
ZA <- readr::read_csv("/home/tim/workspace/AdultCoverage/AdultCoverage/Data/AM_SEG_South Africa_males_3_1.csv")
ZA <- 
ZA %>% 
  mutate(date1 = mdy(date1),
         date2 = mdy(date2))
ggb(ZA, exact.ages=seq(5,80,by=5),lm.method="deming", deaths.summed = TRUE, mig.summed=TRUE)
Moz$id<-1
ggb(Moz, exact.ages=seq(5,75,by=5),lm.method="deming", deaths.summed = FALSE, mig.summed=FALSE)

seg(ZA,maxA=80,minA=15,exact.ages=seq(5,80,by=5), deaths.summed = TRUE, mig.summed=TRUE, eOpen=4.35,exact.ages.ggb=seq(5,80,y=5),delta=TRUE)

seg(Moz, exact.ages=seq(5,75,by=5), deaths.summed = FALSE, mig.summed=FALSE)
seg(ZA)
ggbseg(ZA)
library(DependenciesGraphs)
# Prepare data
dep <- funDependencies("package:DDM","ggbgetRMS")

# visualization
plot(dep)

dep <- envirDependencies("package:DDM")

# visualization
plot(dep)
