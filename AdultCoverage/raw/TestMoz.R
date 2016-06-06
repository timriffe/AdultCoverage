setwd("/home/tim/git/AdultCoverage/AdultCoverage")


library(DDM)
x <- read.csv("Data/Mozambique.csv", stringsAsFactors = FALSE)
x <- addcod(x)
x <- x[,c("cod", "pop1", "pop2", "death", "Age", "Sex", "year1", "year2")]
head(x)
x$year1 <- NULL
x$year2 <- NULL
x$date1 <- as.Date("1997-08-01")
x$date2 <- as.Date("2007-08-01")
x <- x[,c("cod", "pop1", "pop2", "death", "Age", "Sex", "date1", "date2")]
head(x)

ggb(x, exact.ages = seq(15,50,by=5))
ggb(x, exact.ages = seq(30,75,by=5))

bh1(x,exact.ages=seq(15,55,by=5), eOpen = 5)
bh2(x,exact.ages=seq(15,55,by=5), eOpen = 5)


#		
#library(devtools)
#install_github("timriffe/AdultCoverage/AdultCoverage/R/DDM")
#
library(devtools)
load_all("/home/tim/git/AdultCoverage/AdultCoverage/R/DDM")

install_github("timriffe/AdultCoverage/AdultCoverage/R/DDM")
library(DDM)

x <- read.table("/home/tim/Dropbox/paper Lima Riffe and Queiroz/results/UF_estimates/UFdata_females_period1.txt",
		header=TRUE,sep="\t")

group01 <- function(X){
	X[1,c("pop1","pop2","deaths")] <- colSums(X[1:2, c("pop1","pop2","deaths")])
	X <- X[-2, ]
	X
}

x <- do.call(rbind,lapply(split(x,x$cod),group01))
x$sex <- "f"

# age-determination
ggbChooseAges(x[x$cod==53,])
allddm <- ddm(x)
ddmplot(allddm)
#or
ddm(x)

#or
ggb(x)
bh1(x)
bh2(x)


x$cod
#ggb(x, exact.ages =seq(20,60,by=5))
#bh1(x, exact.ages =seq(20,60,by=5))
#bh2(x, exact.ages =seq(20,60,by=5))

data.frame(
		ggb = ggb(x)$coverage, 
		bh1 = bh1(x)$coverage, 
		bh2 = bh2(x)$coverage)
ddm$x <- 1:nrow(ddm)

# some mean functions
h.mean <- function(x){
	1/mean(1/x)
}
g.mean <- function(x){
	n <- length(x)
	prod(x)^(1/n)
}

# now plot results

png("/home/tim/git/AdultCoverage/AdultCoverage/Figures/CompareMethodsTest.png")
plot(ddm$x,ddm$ggb, pch = 19,col="#FFA155", ylim=c(0,5), cex=.6, 
		xlab = "data row", ylab = "coverage estimate",main = "UFdata females period 1",
		panel.first = list(
				segments(ddm$x,apply(ddm[,1:3],1,min),ddm$x,apply(ddm[,1:3],1,max), lwd=.5,col=gray(.5))))
points(ddm$x,ddm$bh1,pch=19,col = "royalblue", cex=.6)
points(ddm$x,ddm$bh2,pch=19,col = "forestgreen", cex=.6)
points(ddm$x, apply(ddm[,1:3],1,g.mean), col = "gray", cex = .6, pch = 19)
points(ddm$x, apply(ddm[,1:3],1,h.mean), col = "magenta", cex = .6, pch = 19)
legend("topright", col = c("#FFA155","royalblue","forestgreen","magenta","gray"), 
		pch = 19, cex=.6, legend = c("GGB","BH1","BH2","Hmean","Gmean"),bty="n")
dev.off()


x <- read.csv("Data/testData53.csv", stringsAsFactors = FALSE)
X <- x[x$cod == 15, ]
X$sex <- "f"
bh1(x,exact.ages = seq(25,55,by=5),eOmega = 5)
bh1(x,exact.ages = ggbgetAgesFit(tab1[[2]]),eOpen = 2.37173661361496)
bh2(x,exact.ages = ggbgetAgesFit(tab1[[2]]),eOpen = 2.37173661361496)

bh1MakeColumns(tab1[[1]],eOpen = 2.37173661361496)$deathLT
bh2MakeColumns(tab1[[2]],agesFit = ggbgetAgesFit(tab1[[2]]),eOpen = 2.37173661361496)$deathLT


x <- read.csv("Data/testData53.csv", stringsAsFactors = FALSE)

library(devtools)
install_github("timriffe/AdultCoverage/AdultCoverage/R/DDM")
document()
library(DDM)
