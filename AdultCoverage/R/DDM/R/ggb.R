
# Author: tim
###############################################################################
# contains functions related to growth balance method


#' @title calcuate the root means square of the error to help find optimal age range
#' 
#' @description Called by \code{ggbgetAgesFit()} whenever the user doesn't want to manually determine the age range used to determine registration coverage. Probably no need to be called by top-level users. If a user would rather determine the optimal age range some other way, then look to \code{ggbcoverageFromYear()} where \code{ggbgetRMS} is called and add another condition or make it call something else.
#' 
#' @param agesi the set of ages used for this iteration
#' @param codi \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, and \code{$age}. 
#' 
#' @return the RMSE
#' 
#' @export 
#' 

ggbgetRMS <- function(agesi, codi){
	codi <- ggbFittedFromAges(codi, agesfit = agesi)
	# get root mean square of residuals
	# I'd do this with magrittr, but just to stay dependency-free...
	sqrt(
		sum(
			  (
				 (codi$leftterm[codi$age %in% agesi] - codi$fitted[codi$age %in% agesi]) ^ 2
			   )
	        ) / length(agesi)
         )
	# no checks here for NAs...
}


#'
#' @title estimate death registration coverage for a single year/sex/region using the GGB method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method is based on finding a best-fitting linear relationship between two modeled parameters (right term and left term), but the fit, and resulting coverage estimate, depend on exactly which age range is taken. This function either finds a nice age range for you automatically, or you can specify an exact vector of ages. Called by \code{ggb()}. Users probably don't need to call this directly. Just use \code{ggb()} instead.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. 
#' 
#' @param codi \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, and \code{$age}.
#' @param exact.ages. optional. use an exact set of ages to estimate coverage.
#' @param minA. the minimum of the age range searched. Default 15
#' @param maxA. the maximum of the age range searched. Default 75
#' @param minAges. the minimum number of adjacent ages needed as points for fitting. Default 8
#' 
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. 
#' 
#' @export

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


#'
#' @title make the growth-adjusted quasi lifetable columns required by GGB method
#' 
#' @description Called by \code{plot.ggb()} and \code{ggbcoverageFromYear()}. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. If columns produced by \code{ggbMakeColumns()} are not present, then we call it here just to keep things from breaking.
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{ggb()}
#' @param minA. the minimum of the age range searched. Default 15
#' @param maxA. the maximum of the age range searched. Default 75
#' 
#' @return codi, with many columns added, most importantly \code{$rightterm}, \code{$leftterm}, and \code{$exclude}.
#' 
#' @export

ggbFittedFromAges <- function(codi, agesfit){
	
	if (! "leftterm" %in% colnames(codi)){
		codi <- ggbMakeColumns(codi)
	}
	# assumes ggbMakeColumns() has been run.
	slope       <- with(codi, 
			sd(leftterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
	)
	intercept   <-  with(codi, 
			(mean(leftterm[age %in% agesfit]) * slope - mean(rightterm[age %in% agesfit]))
	) 
	codi$fitted <- codi$rightterm * slope + intercept
	codi
}

#' @title given a set of ages, what is the implied death registration coverage?
#' 
#' @description For a single year/sex/region of data (formatted as required by \code{ggb()}), what is the registration coverage implied by a given age range? Called by \code{ggbcoverageFromYear()} and \code{plot.ggb()}.

ggbcoverageFromAges <- function(codi, agesfit){
	if (! "leftterm" %in% colnames(codi)){
		codi <- ggbMakeColumns(codi, minA. = min(agesfit), maxA. = max(agesfit))
	}
	# this is the coverage estimate
	1 / with(codi, sd(leftterm[age %in% agesfit]) / sd(rightterm[age %in% agesfit]))
}


#'
#' @title make the growth-adjusted quasi lifetable columns required by GGB method
#' 
#' @description Called by \code{plot.ggb()} and \code{ggbcoverageFromYear()}. This simply modulates some cod ethat would otherwise be repeated. Users probably don't need to call this function directly. 
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{ggb()}
#' @param minA. the minimum of the age range searched. Default 15
#' @param maxA. the maximum of the age range searched. Default 75
#' 
#' @return codi, with many columns added, most importantly \code{$rightterm}, \code{$leftterm}, and \code{$exclude}.
#' 
#' @export
ggbMakeColumns <- function(codi, minA. = 15, maxA. = 75){
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
	codi$leftterm          <- (codi$birthdays / codi$Tx) - codi$cumgrowth
	# certain columns can be safely ignored in future operations
	codi$exclude          <-  codi$Tx != 0 & codi$birthdays != 0 & codi$age >= minA. & codi$age <= 75
	codi
}


#'
#' @title determine the age range that minimizes the mean squared error
#' @description Called by \code{ggbcoverageFromYear()} whenever \code{exact.ages} are not given. This automates what one typically does visually.
#' @seealso plot.ggb
#' 
#' @param

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
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimate for a variety of intercensal periods, then stack them, and use a variable called \code{$cod} with unique values for each data chunk. Different values of \code{$cod} could indicate sexes, regions, intercensal periods, etc.
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, and \code{$age}.
#' @param minA
#' @param maxA
#' @param minAges
#' @param exact.ages
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. Rows indicate data partitions, as indicated by the optional \code{$cod} variable.
#' 
#' @export

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
			sd(leftterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
	)
	intercept   <-  with(codi, 
			(mean(leftterm[age %in% agesfit]) * slope - mean(rightterm[age %in% agesfit]))
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
	codi     <- ggbFittedFromAges(codi, agesfit)
	# this is the basic formula
	coverage <- ggbcoverageFromAges(codi, agesfit)
	
	# some objects used throughout
	age     <- codi$age
	leftt   <- codi$leftterm
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
			# regenerate $fitted column of codi using new ages
			codi     <- ggbFittedFromAges(codi, agesfit.)
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




