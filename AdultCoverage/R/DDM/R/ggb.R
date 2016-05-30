
# Author: tim
###############################################################################
# contains functions related to growth balance method

ggbgetRMS <- function(agesi, codi){
	slope       <- with(codi, 
			sd(lefterm[age %in% agesi]) /  sd(rightterm[age %in% agesi])
	)
	intercept   <-  with(codi, 
			(mean(lefterm[age %in% agesi]) * slope - mean(rightterm[age %in% agesi]))
	) 
	codi$fitted <- codi$rightterm * slope + intercept
	# get root mean square of residuals
	sqrt(sum(((codi$lefterm[codi$age %in% agesi] - codi$fitted[codi$age %in% agesi])^2))/length(agesi))
	
}

ggbcoverageFromYear <- function(codi, exact.ages., minA., maxA., minAges.){
	
	# if exact.ages is given, we override other age-parameters
	if (!is.null(exact.ages.) & length(exact.ages.) >= 3){
		if (min(exact.ages.) < minA.){
			minA. <- min(exact.ages.)
		} 
		if (max(exact.ages.) > maxA.){
			maxA. <- max(exact.ages.)
		}
		if (minAges. < length(exact.ages.)){
			minAges. <- length(exact.ages.)
		}
	}
	
	codi    <- ggbMakeColumns(codi, minA., maxA.)
	
	if (!is.null(exact.ages.) & length(exact.ages.) >= 3){
		agesfit <- exact.ages.
	} else {
		agesfit <- ggbgetAgesFit(codi, minAges.)
	}
		
	# this is the basic formula
	coverage <- ggbcoverageFromAges(codi, agesfit)
	data.frame(cod = unique(codi$cod), coverage = coverage, lower = min(agesfit), upper = max(agesfit))
}

ggbColumnsFromAges <- function(codi, agesfit){
	# assumes ggbMakeColumns() has been run.
	slope       <- with(codi, 
			sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
	)
	intercept   <-  with(codi, 
			(mean(lefterm[age %in% agesfit]) * slope - mean(rightterm[age %in% agesfit]))
	) 
	codi$fitted <- codi$rightterm * slope + intercept
	codi
}

ggbcoverageFromAges <- function(codi, agesfit){
	codi <- ggbColumnsFromAges(codi, agesfit)	
	# this is the coverage estimate
	1 / with(codi, sd(lefterm[age %in% agesfit]) / sd(rightterm[age %in% agesfit]))
}

ggbMakeColumns <- function(codi, minA., maxA.){
	AgeInt                 <- detectAgeInterval(codi, MinAge =  minA., MaxAge = maxA., ageColumn = "age")
	ages                   <- codi$age
	# TR make this accept exact dates.
	dif.                   <- yint2(codi)
	
	codi$pop1cum           <- rev(cumsum(rev(codi$pop1)))  # like Tx
	codi$pop2cum           <- rev(cumsum(rev(codi$pop2)))  # like Tx
	codi$deathcum          <- rev(cumsum(rev(codi$deaths))) # like lx
	
	# define new column for birthdays between pop estimates
	# TR: loop removed
	codi$birthdays <- c(0, sqrt(codi$pop1[ -nrow(codi) ] * codi$pop2[ -1 ])) / AgeInt
	
	# create 'stationary' Tx as geometric avg of within-cohort consecutive ages
	codi$Tx               <- sqrt(codi$pop1cum * codi$pop2cum)
	
	# growth rate per annum
	codi$cumgrowth        <- log(codi$pop2cum / codi$pop1cum) / dif.
	
	# eqns from formula in Hill/Horiuchi
	codi$rightterm        <- codi$deathcum / codi$Tx
	# eqns from formula in Hill/Horiuchi
	codi$lefterm          <- (codi$birthdays / codi$Tx) - codi$cumgrowth
	# certain columns can be safely ignored in future operations
	codi$exclude          <-  codi$Tx != 0 & codi$birthdays != 0 & codi$age >= minA. & codi$age <= 75
	codi
}

ggbgetAgesFit <- function(codi, minAges.){
		
	maxAges   <- sum(codi$exclude)
	agesUniv  <- codi$age[codi$exclude]
	
	FirstAges <- agesUniv[agesUniv < 30]
	
	ind       <- 0
	agesL     <- list()
	# determine ages to test
	for (Nr in maxAges:minAges.){ #
		its <- length(agesUniv) - Nr + 1
		for (set in 1:its){ # 
			ind <- ind + 1
			agesL[[ind]] <- agesUniv[set:(set+Nr-1)]
		}
	}
	
	# these are the ages that give the best r2 or RMS
	
	agesfit     <- agesL[[which.min(unlist(lapply(agesL, ggbgetRMS, codi = codi)))]]
	agesfit
}


#'
#' @title estimate death registration coverage using the GGB method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method is based on finding a best-fitting linear relationship between two modeled parameters (right term and left term), but the fit, and resulting coverage estimate, depend on exactly which age range is taken. This function either finds a nice age range for you automatically, or you can specify an exact vector of ages. 
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month.
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, and \code{$age}
#' 


ggb <- function(X, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL){         ##  Data
	##  Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)
	tab         <- data.frame(X)           
	colnames(tab) <- tolower(colnames(tab))
	
	# in case there is no splitting var, this way we split anyway.
	tab         <- addcod(tab)
	
	# guess which column is the deaths column, rename it deaths
	tab         <- guessDeathsColumn(tab)
	# TR: account for decimal intervals
	tab$pop1    <- as.double(tab$pop1)
	tab$pop2    <- as.double(tab$pop2)
	tab$deaths  <- as.double(tab$deaths)
	
	# this splitting variable, invented if necessary
	tab1        <- split(tab, factor(tab$cod))

	# iterate over whatever it happens to be: regions, years
	coverages   <- as.data.frame(do.call(rbind,lapply(
					    tab1, 
					    ggbcoverageFromYear, 
					    exact.ages. = exact.ages,
						minA. = minA, 
						maxA. = maxA,
						minAges. = minAges
					)))
	
	# this has cod as a column, but no year, sex. 
	return(coverages)
}

############################################################
# Note: GGB is used to determine ages for all methods here #
# ergo, these plotting functions are only implemented for  #
# GGB method                                               #
############################################################

# functions related to ggb plotting.

#'
#' @title does a given pairlist of x and y coordinates fall within the plot region?
#' 
#' @description Check to see if a point clicked falls in the plot or outside it. This function is used by \code{plot.ggb()}.
#' 
#' @param USR, as given by \code{par("usr")}
#' @param click a pairlist with elements \code{$x} and \code{$y}, as returned by \code{locator(1)}
#' 
#' @return logical. \code{TRUE} if in the plot region.
#' 
#' @export
inUSR <- function(USR, click){
	x <- click$x
	y <- click$y
	if (x < USR[1] | x > USR[2] | y < USR[3] | y > USR[4]){
		return(FALSE)
	}
	TRUE
}

#'
#' @title which age is closest to the point clicked?
#' @description a utility function called by \code{plot.ggb()}.
#' 
#' @param xvec
#' @param yvec
#' @param click
#' @param

guessage <- function(xvec,yvec,click,age){
	age[which.min(sqrt((xvec - click$x) ^ 2 + (yvec - click$y) ^ 2))]
}

adjustages <- function(a, age, agesfit){
	maxa <- max(agesfit)
	mina <- min(agesfit)
	
	isitamin <- abs(a - mina) < abs(a - maxa)
	if (isitamin){
		mina <- a
	} else {
		maxa <- a
	}
	age[age >= mina & age <= maxa]
}


slopeint <- function(codi, agesfit){
	slope       <- with(codi, 
			sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
	)
	intercept   <-  with(codi, 
			(mean(lefterm[age %in% agesfit]) * slope - mean(rightterm[age %in% agesfit]))
	) 
	list(b=slope, a=intercept)
}

plot.ggb <- function(codi, exact.ages = NULL,minA = 15, maxA = 75, minAges = 8, maxit=15){
	# this is the automatic age selection.

	if (!is.null(exact.ages) & length(exact.ages) >= 3){
		if (min(exact.ages) < minA){
			minA <- min(exact.ages)
		} 
		if (max(exact.ages) > maxA){
			maxA <- max(exact.ages)
		}
		if (minAges < length(exact.ages)){
			minAges <- length(exact.ages)
		}
	}
	
	# 
	codi    <- ggbMakeColumns(codi, minA, maxA)
	
	if (!is.null(exact.ages) & length(exact.ages) >= 3){
		agesfit <- exact.ages
	} else {
		agesfit <- ggbgetAgesFit(codi, minAges)
	}
	
	
	si       <- slopeint(codi, agesfit)
	codi     <- ggbColumnsFromAges(codi, agesfit)
	# this is the basic formula
	coverage <- ggbcoverageFromAges(codi, agesfit)
	
	# some objects used throughout
	age     <- codi$age
	leftt   <- codi$lefterm
	rightt  <- codi$rightterm
	
	# age ranges used for fitting
	amin    <- min(agesfit); amax <- max(agesfit)
	plot(rightt, leftt, asp = 1, pch = 19, col = "#00000050", cex = 1.6,
			xlab = "right term",
			ylab = "left term",
			main = paste0("Age range [",amin,
					",",amax,"], est. coverage = ",round(coverage*100,1)),
			sub = "(optimized age range)")
	# automatically fit line (RMS of ggb)
	abline(a = si$a, b = si$b, col = "blue")
	# shows points used to fit line
	points(rightt[age %in% agesfit], 
			leftt[age %in% agesfit], col = "#FFFF00", pch = 19, cex = 1.6)
	text(rightt, leftt, age, cex = .6)
	legend("bottomright",lty=1,col="blue",legend="fitted line",bty="n")
	# message to user
	mtext("*Adjusting ages*\nClick any age to change bounds for fitting\nClick in margin to quit",
			side = 3, line = -3, outer = FALSE)
	
	
	USR    <- par()$usr
	agesfit. <- agesfit # keep old one...
	for (i in 1:maxit){
		# now enter into interactive loop
		clicki <- locator(1)
		# was the click inside the plot or outside the plot?
		IN     <- inUSR(USR, clicki)
		if (IN){
			
			# if it's inside the plot, then which is the closest age?
			a        <- guessage(rightt,leftt,clicki,age)
			# guess how to adjust ages based on which age was clicked
			agesfit. <- adjustages(a, age, agesfit.)
			# new range of ages used for fitting
			amin     <- min(agesfit.)
			amax     <- max(agesfit.)
			# regenerate columns of codi using new ages
			codi     <- ggbColumnsFromAges(codi, agesfit.)
			# an estimate of the resulting coverage
			coverage <- ggbcoverageFromAges(codi, agesfit.)
			# get params for abline..
			si       <- slopeint(codi, agesfit.)
			
			# redraw plot
			plot(   rightt, 
					leftt, 
					asp = 1, 
					pch = 19, 
					col = "#00000050", 
					cex = 1.6,
					xlab = "right term",
					ylab = "left term",
					main = paste0("Age range [", amin,
							",", amax, "], est. coverage = ",round(coverage * 100, 1)),
					sub = "(optimized age range)")
			# new fitted slope, intercept
			abline(a = si$a, b = si$b, col = "blue")
			# indicate which points used with color
			points(rightt[age %in% agesfit.], 
					leftt[age %in% agesfit.], col = "#FFFF00", pch = 19, cex = 1.6)
			text(rightt, leftt, age, cex = .6)
			legend("bottomright", lty = 1, col = "blue", legend = "fitted line", bty = "n")
		} else {
			break
		}
	}
	
	# click outside the margin to save results,
	# in this case return potentially useful stuff 
	list(agesfit = agesfit., codi = codi, coverage = coverage)
}

#
#makeggb <- function(x){
#	
#}
#codi <- DM3[DM3$cod == 53, ]
#class(codi) <- "ggb"
#plot(codi, maxit=25)




