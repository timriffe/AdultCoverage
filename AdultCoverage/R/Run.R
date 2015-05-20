# TODO: Add comment
# 
# Author: triffe
###############################################################################


ITA <- read.table("/hdir/0/triffe/workspace/AdultCoverage/Data/italy.test.txt", 
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(ITA)
#Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)
x <- with(ITA, data.frame(
    cod = Year, 
    age = Age, 
    pop1 = Pop.fem, 
    year1 = Year, 
    pop2 = Pop.fem.1, 
    year2 = Year.1, 
    death = (Death.fem + Death.fem.1)/2) )
x$age[x$age == 4] <- 1
# coerce to required format:

