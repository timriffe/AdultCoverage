
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
ddm(Moz,delta=TRUE)
ddm(BrasilMales)
ggb(ZA,
      deaths.summed=TRUE,
      mig.summed =TRUE,
      exact.ages=seq(5,80,by=5),
      lm.method = "oldschool")

load_all()
ggbMakeColumns(ZA) %>% View()
  
seg(ZA,
    deaths.summed=TRUE,
    mig.summed =TRUE,
    exact.ages=seq(25,60,by=5),
    exact.ages.ggb=seq(5,80,by=5),
    delta=TRUE,
    lm.method="oldschool",
    eOpen=4.34676243530911)


segMakeColumns(ZA,
               deaths.summed=TRUE,
               mig.summed =TRUE,
               exact.ages.ggb=seq(5,80,by=5),
               delta=TRUE,
               lm.method="ols",
               eOpen=4.35) %>% 
  View()

segCoverageFromYear(ZA,
               deaths.summed=TRUE,
               mig.summed =TRUE,
               exact.ages.ggb=seq(5,80,by=5),
               delta=TRUE,
               lm.method="tls",
               eOpen=4.35)$codi %>% View()

seg(ZA,
    deaths.summed=TRUE,
    mig.summed =TRUE,
    exact.ages.ggb=seq(5,80,by=5),
    delta=TRUE,
    lm.method="tls",
    eOpen=4.35)

#--------------------
library(dplyr)
library(magrittr)

# change data
ggb(Moz)
seg(Moz,delta=F)
ggbseg(Moz)
res <- ddm(Moz,delta=TRUE)

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









minA = 15 
maxA = 75 
minAges = 8 
exact.ages.ggb = NULL
exact.ages.seg = NULL 
eOpen = NULL
deaths.summed = FALSE
mig.summed = deaths.summed
lm.method = "oldschool"
nx.method = 2
delta = FALSE
