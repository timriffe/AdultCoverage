
# Author: tim
###############################################################################

#'
#' @title run all three deaths registration coverage estimation methods
#' 
#' @description Estimate the generalized growth balance method, and the two Bennett-Horiuchi methods of estimating death registration coverage. This requires two censuses and an estimate of the deaths in each 5-year age group between censuses. This might be the arithmetic average of deaths in each age class, or simply the average of deaths around the time of the two censuses. All methods use some stable population assumptions. 

#' @details All methods require some specification about which age range to base results on. If not given, an optimal age range will be estimated automatically, and this information is returned to the user. To identify an age-range in the visually, see \code{plot.ggb()}, when working with a single year/sex/region of data. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared residual between the fitted line and right term.If you want coverage estimates for a variety of partitions (intercensal periods/regions/by sex), then stack them, and use a variable called \code{$cod} with unique values for each data partition. If data is partitioned using the variable \code{$cod}, then the age range automatically determined might not be the same for each partition. If user-specified, (using a vector of \code{exact.ages}) the age ranges will be the same for all partitions. If you want to specify particular age ranges for each data partition, then you'll need to loop it somehow. 
#' 
#' All three methods require time points of the two censuses. Census dates can be given in a variety of ways: 1) (preferred) using \code{Date} classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing respective integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month.  Different values of \code{$cod} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths for each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. 
#' 
#' The Bennett-Horiuchi methods require an estimate of remaining life expectancy in the open age group of the data provided. This is produced using a standard reference to the Coale-Demeny West model life tables. If 

#' @references Need to cite stuff here.

ddm <- function(X, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL, eOpen = NULL){
	ggbres <- ggb(X = X, 
					minA = minA, 
					maxA = maxA, 
					minAges = minAges, 
					exact.ages = exact.ages)
	bh1res <- bh1(X = X, 
					minA = minA, 
					maxA = maxA, 
					minAges = minAges, 
					exact.ages = exact.ages, 
					eOpen = eOpen)
	bh2res <- bh2(X = X, 
					minA = minA, 
					maxA = maxA, 
					minAges = minAges, 
					exact.ages = exact.ages, 
					eOpen = eOpen)
	# return all results
	results <- data.frame(	cod = ggbres$cod,
				ggb = ggbres$coverage,
				bh1 = bh1res$coverage,
				bh2 = bh2res$coverage,
				lower = ggbres$lower,
				upper = ggbres$upper
			)
	class(results) <- "ddm"
	results
}

plot.ddm <- function(X, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL, eOpen = NULL){
	if (class(X) == "data.frame"){
		X <- ddm(X, minA = minA, 
				maxA = maxA, 
				minAges = minAges, 
				exact.ages = exact.ages, 
				eOpen = eOpen)
	} 
	
	X$x <- 1:nrow(X)
	
	# these can move to utils.R if we decide it's worth keeping as a summary method
	h.mean <- function(x){
		1/mean(1/x)
	}
	g.mean <- function(x){
		n <- length(x)
		prod(x)^(1/n)
	}
	plot(X$x,X$ggb, pch = 19,col="#FFA155", ylim=range(X[,c("ggb","bh1","bh2")]), cex=.6, 
			xlab = "data row", ylab = "coverage estimate",main = "UFdata females period 1",
			panel.first = list(
					segments(X$x,apply(X[,1:3],1,min),X$x,apply(X[,1:3],1,max), lwd=.5,col=gray(.5))))
	points(X$x,X$bh1,pch=19,col = "royalblue", cex=.6)
	points(X$x,X$bh2,pch=19,col = "forestgreen", cex=.6)
	#points(X$x, apply(X[,1:3],1,g.mean), col = "gray", cex = .6, pch = 19)
	points(X$x, apply(X[,1:3],1,h.mean), col = "magenta", cex = .6, pch = 19)
	legend("topright", col = c("#FFA155","royalblue","forestgreen","magenta"), 
			pch = 19, cex=.6, legend = c("GGB","BH1","BH2","Hmean"),bty="n")
	
	
}

