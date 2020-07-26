
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
getwd()
load_all("AdultCoverage/R/DDM")
document("AdultCoverage/R/DDM")
res <- seg(Moz)
res
# The Brasil data
BM <- seg(BrasilMales)
BF <- seg(BrasilFemales)
ggbseg(Moz)
ggb(Moz)
seg(Moz)

ggb_res    <- ggb(BrasilMales)
seg_res    <- seg(BrasilMales)
ggbseg_res <- ggbseg(BrasilMales)

ggb_res$method <- "GGB"
seg_res$method <- "SEG"
ggbseg_res$method <- "GGBSEG"
colnames(ggbseg_res) <- c("cod","coverage", "lower.ggb","upper.ggb","lower","upper","method")
res <- rbind(ggb_res[,c("cod","method","lower","upper")],
             seg_res[,c("cod","method","lower","upper")],
             ggbseg_res[,c("cod","method","lower","upper")]
             )
library(tidyverse)
library(ggplot2)
res %>% 
  as.data.frame() %>% 
  ggplot(aes(x=as.factor(cod),ymin=lower, ymax = upper,color = method)) + 
           geom_errorbar(position = "dodge")

segplot(subset(BrasilMales, cod == 51))
ggbsegplot
ggbChooseAges(subset(BrasilMales, cod == 51))
