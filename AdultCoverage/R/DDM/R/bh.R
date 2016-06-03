
# Author: tim
###############################################################################
# contains functions related to Bennet-Horiuchi methods

#'
#' @title estimate remaining life expectancy in the open age group
#' @description This calculation is based on an indirect method to reference the Coale-Demeny West model lifetable. First one makes a psuedo lifetable deaths column using some stable pop assumptions (different in BH1 vs BH2). Then take the ratio of the sum of ages 10-39 to 40-59. These ratios have been worked out for each model lifetable level, so we can pick the level based on the ratio we produce from the data. From there, we just pick out the remaining life expectancy that corresponds to the top age in our data, which for now hopefully isn't higher than 95. The model lifetables don't go higher than 95 for now, but that's well beyond the range for this method family. If your data go beyong 85 or so, then just group down to 85, say, and estimate using that instead of keeping a high open age. Called by \code{bh1MakeColumns()} and \code{bh2MakeColumns()}, and not intended for direct user interface, because you need to produce the \code{$deathsLT} column. You can skip calling this function by specifying \code{eOpen} in the top call to \code{bh1()} or \code{bh2()}.
#' 
#' @param codiaugmented the standard codi object being passed through the chain, but having been pre-processed in the course of \code{bh1MakeColumns()} or \code{bh2MakeColumns()}
#' 
#' @return numeric an estimate of remaining life expectancy in the open age group
#' 
#' @export

eOpenCD <- function(codiaugmented){
	# this needs to have the column deathsLT, as passed in 
	# by either bh1MakeColumns() or bh2MakeColumns(). These
	# are slightly different, but this piece of code is identical and therefore 
	# modularized.
	# this throws an error if sex isn't coded as expected
	sex        <- detectSex(Dat = codiaugmented, 
			sexColumn = "sex")
	ages       <- codiaugmented$age
	
	ratio      <- with(codiaugmented,
			sum(deathLT[ages %in% c(10:39)]) / 
					sum(deathLT[ages %in% c(40:59)]))
	
	if (sex == "f"){
		standardratios <- c(1.376,	1.3,	1.233,	1.171,
				1.115,	1.062,	1.012,	0.964,
				0.918,	0.872,	0.827,	0.787,
				0.729,	0.673,	0.617,	0.56,
				0.501,	0.438,	0.365,	0.298,
				0.235,	0.175,	0.117)
		ex <- cdmltw("F")$ex
	}
	if (sex == "m"){
		# need to change this:
		standardratios <- c(1.161,	1.094,	1.034,	0.98,	
				0.93,	0.885,	0.842,	0.802,	
				0.763,	0.725,	0.689,	0.648,
				0.609,	0.57,	0.53,	0.49,	
				0.447,	0.401,	0.352,	0.305,	
				0.255,	0.202,	0.147)
		
		ex <- cdmltw("M")$ex
	}
	
	# this spits back a C-D level that is a decimal estimate
	AllLevels <- 3:25
	CDlevel   <- splinefun(AllLevels~standardratios)(ratio)
	
	# which is the open age?
	OA        <- max(ages)
	availages <- as.integer(colnames(ex))
	
	stopifnot(OA %in% as.integer(colnames(ex)))
	
	eOpen     <- splinefun(ex[,as.character(OA)]~1:25)(CDlevel)
}


#' @title given a set of ages, what is the implied death registration coverage?
#' 
#' @description For a single year/sex/region of data (formatted as required by \code{bh1()}, \code{bh2()}), what is the registration coverage implied by a given age range? Called by \code{bh1CoverageFromYear()} and  \code{bh2CoverageFromYear()}. Here, the function simply takes the arithmetic mean of a given age range of \code{$Cx}, as returned by \code{bh1MakeColumns()} or \code{bh2MakeColumns()}. Not intended for top-level use.
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{ggb()}
#' @param agesFit an integer vector of ages, either returned from \code{ggbgetAgesFit} or user-supplied.
#' @param minA the minimum of the age range searched. Default 15
#' @param maxA the maximum of the age range searched. Default 75
#' @param eOpen optional. A value for remaining life expectancy in the open age group.
#'
#' @return numeric. the estimated level of coverage.
#' @export

bhCoverageFromAges <- function(codi, agesFit){
	stopifnot("Cx" %in% colnames(codi))
	inds    <- codi$age %in% agesFit
	sum(codi$Cx[inds]) / length(agesFit)
}

#'
#' @title make the Bennett-Horiuchi quasi lifetable columns required by the estimation method
#' 
#' @description Called by \code{bh1CoverageFromYear()}. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. 
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{bh1()}
#' @param minA the minimum of the age range searched. Default 15
#' @param maxA the maximum of the age range searched. Default 75
#' @param eOpen optional. A value for remaining life expectancy in the open age group.
#' 
#' @return codi, with many columns added, most importantly \code{$Cx}.
#' 
#' @export

bh1MakeColumns <- function(codi, minA = 15, maxA = 75, eOpen = NULL){
	

	# attempt to detect AgeInterval, should be obvious. And really, we only consider 5-yr intervals.
	# if this were done w single ages minAges would need to increase to at least 35 I guess.
	AgeInt                 <- detectAgeInterval(Dat = codi, MinAge =  minA, MaxAge = maxA, ageColumn = "age")
	ages                   <- codi$age
	# ages better be unique! I expect this will work unless $cod is mis-specified. 
	# This is a good catch for that
    stopifnot(length(ages) == length(unique(ages)))
	
	# codi can use date columns, or year, month, day columns...
	dif                    <- yint2(X = codi)
	N                      <- nrow(codi)
	
	codi <- within(codi,{
				# birthdays, as in GGB
				birthdays <- c(0, sqrt(pop1[ -N  ] * pop2[ -1 ])) / AgeInt
				# age-specific growth
	        	growth    <- log(pop2 / pop1) / dif
				growth[is.infinite(growth) | is.nan(growth)] <- 0
				# cumulative growth
				cumgrowth <- AgeInt * c(0,cumsum(growth[ -N ])) + AgeInt / 2 * growth
				deathLT   <- deaths * exp(cumgrowth)
			})
	
	if (is.null(eOpen)){
		eOpen <- eOpenCD(codiaugmented = codi)
	} else {
		eOpen <- eOpen
	}
	
	eON   <- eOpen * codi$growth[N]

	codi <- within(codi, {
				pop_a <- 0
				pop_a[N] <- deaths[N] * exp(eON) - (eON ^ (2/6))
				for(j in N:2){
					pop_a[j - 1] <- pop_a[j] * exp(AgeInt * growth[j - 1]) + 
							deaths[j - 1] * exp(AgeInt / 2 * growth[j - 1])
				}
				rm(j)
				Cx <- pop_a / birthdays
			})
	
	################
	codi
}

#'
#' @title estimate death registration coverage for a single year/sex/region using the Bennet-Horiuchi method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model lifetable, although the user may also supply a value. Called by \code{bh1()}. Users probably do not need to use this function directly.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. 
#' 
#' @param codi \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$sex}, \code{$age}, and \code{$cod} (to indicate regions, periods, sexes).
#' @param exact.ages optional. use an exact set of ages to estimate coverage.
#' @param minA the minimum of the age range searched. Default 15
#' @param maxA the maximum of the age range searched. Default 75
#' @param minAges the minimum number of adjacent ages needed as points for fitting. Default 8
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. 
#' 
#' @export


bh1CoverageFromYear <-  function(codi, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL, eOpen = NULL){        ##  Data
	# if exact.ages is given, we override other age-parameters
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
	
	# outsource automatic age selection to GGB method.
	if (is.null(exact.ages)){
		agesFit <- ggbgetAgesFit(codi = codi, minA = minA, maxA = maxA, minAges = minAges)
	} else {
		agesFit <- exact.ages
	}
	
	codi     <- bh1MakeColumns(	codi = codi, 
								minA = minA, 
								maxA = maxA,
								eOpen = eOpen )
	
	coverage <- bhCoverageFromAges(codi = codi, agesFit = agesFit)

	data.frame(cod = unique(codi$cod), coverage = coverage, lower = min(agesFit), upper = max(agesFit))
}

#'
#' @title estimate death registration coverage using the Bennet-Horiuchi method 
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model lifetable, although the user may also supply a value.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} columns are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$cod} with a unique values for each data chunk. Different values of \code{$cod} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths in each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range in the traditional visual way, see \code{plot.ggb()}, when working with a single year/sex/region of data. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared residual between the fitted line and right term. Finally, only specify \code{eOpen} when working with a single region/sex/period of data, otherwise the same value will be passed in irrespective of mortality and sex.
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$age}, \code{$sex}, and \code{$cod} (if there are more than 1 region/sex/intercensal period).
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param exact.ages optional. A user-specified vector of exact ages to use for coverage estimation
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. Rows indicate data partitions, as indicated by the optional \code{$cod} variable.
#' 
#' @export
bh1 <- function(X, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL, eOpen = NULL){
	
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
	
	tab1        <- split(tab, tab$cod)
	
	coverages <- as.data.frame(
					do.call(
						rbind,
						lapply(
							tab1, 
							bh1CoverageFromYear, 
							minA = minA, 
							maxA = maxA,
							minAges = minAges,  
							exact.ages = exact.ages,
							eOpen = eOpen
            )))
	#return(data.frame(Coverage = coverages,correctionFactor = 1/coverages))
	
	coverages
}

###################################################################################
# End BH1 functions
###################################################################################
# Start BH2 functions
###################################################################################



#'
#' @title make the Bennett-Horiuchi quasi lifetable columns required by the second estimation method
#' 
#' @description Called by \code{bh2CoverageFromYear()}. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. 
#' @details \code{agesFit} is a vector passed in from \code{bh2MakeColumns()}, and it was either estimated using the GGB automatic method there, or simply came from the argument \code{exact.ages} specified in \code{bh2()}. By default we just automatically estimate these. \code{eOpen} can be either user-specified, or it will be estimated automatically using \code{eOpenCD()}.
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{bh2()}
#' @param minA the minimum of the age range searched. Default 15
#' @param maxA the maximum of the age range searched. Default 75
#' @param eOpen optional. A value for remaining life expectancy in the open age group.
#' @param agesFit vector of ages as passed in by \code{bh2coverageFromYear)} 
#' @return codi, with many columns added, most importantly \code{$Cx}.
#' 
#' @export

bh2MakeColumns <- function(
		codi, 
		minA = 15, 
		maxA = 75,  
		agesFit, 
		eOpen = NULL){
	# this throws an error if sex isn't coded as expected
	sex                    <- detectSex(Dat = codi, sexColumn = "sex")
	
	AgeInt       <- detectAgeInterval(
							Dat = codi, 
							MinAge =  minA, 
							MaxAge = maxA, 
							ageColumn = "age")
	ages         <- codi$age
	dif          <- yint2(X = codi)
	
	# just get left term / right term
	# slopeint() is in the ggb() family
	ab           <- slopeint(codi = codi, agesfit = agesFit)
	
	# the only difference between this method and bh1 is that
	# in the next couple lines we use pop1adj instead of pop1...
	relcomp      <- exp(ab$a * dif)

	codi <- within(codi,{
				pop1adj   <- pop1 / relcomp
				# birthdays, as in GGB
				birthdays <- c(0, sqrt(pop1adj[ -N  ] * pop2[ -1 ])) / AgeInt
				# age-specific growth
				growth    <- log(pop2 / pop1adj) / dif
				growth[is.infinite(growth) | is.nan(growth)] <- 0
				# cumulative growth
				cumgrowth <- AgeInt * c(0,cumsum(growth[ -N ])) + AgeInt / 2 * growth
				deathLT   <- deaths * exp(cumgrowth)
			})
	
	if (is.null(eOpen)){
		eOpen <- eOpenCD(codiaugmented = codi)
	} else {
		eOpen <- eOpen
	}
	
	eON   <- eOpen * codi$growth[N]
	
	codi <- within(codi, {
				pop_a    <- 0
				pop_a[N] <- deaths[N] * exp(eON) - (eON ^ (2/6))
				for(j in N:2){
					pop_a[j - 1] <- pop_a[j] * exp(AgeInt * growth[j - 1]) + 
							deaths[j - 1] * exp(AgeInt / 2 * growth[j - 1])
				}
				rm(j)
				Cx <- pop_a / birthdays
			})
	
	codi
}


#'
#' @title estimate death registration coverage for a single year/sex/region using the modified Bennet-Horiuchi method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model lifetable, although the user may also supply a value. The difference between this method and \code{bh1()} is that here we adjust census 1 part way through processing, based on some calcs similar to GGB. Called by \code{bh2()}. Users probably do not need to use this function directly.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. 
#' 
#' @param codi \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$sex}, \code{$age}, and \code{$cod} (to indicate regions, periods, sexes).
#' @param exact.ages optional. use an exact set of ages to estimate coverage.
#' @param minA the minimum of the age range searched. Default 15
#' @param maxA the maximum of the age range searched. Default 75
#' @param minAges the minimum number of adjacent ages needed as points for fitting. Default 8
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. 
#' 
#' @export


bh2coverageFromYear <- function(codi, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL){
	codiggb      <- ggbMakeColumns(codi = codi, minA = minA, maxA = maxA)
	# Get age range using the GGB auto fitting
	if (is.null(exact.ages)){
		agesFit <- ggbgetAgesFit(codi = codiggb, 
								minA = minA, 
								maxA = maxA, 
								minAges = minAges)
	} else {
		agesFit <- exact.ages
	}
	
	codi         <- bh2MakeColumns(
								codi = codiggb, 
								minA = minA, 
								maxA = maxA, 
								sex = sex, 
								agesFit = agesFit)
	coverage <- bhCoverageFromAges(codi = codi, agesFit = agesFit )
	data.frame(cod = unique(codi$cod), coverage = coverage, lower = min(agesFit), upper = max(agesFit))
}


#'
#' @title estimate death registration coverage using the adjusted Bennet-Horiuchi method 
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model lifetable, although the user may also supply a value. The difference between this method and \code{bh1()} is that here we adjust census 1 part way through processing, based on some calcs similar to GGB.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} columns are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$cod} with a unique values for each data chunk. Different values of \code{$cod} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths in each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range in the traditional visual way, see \code{plot.ggb()}, when working with a single year/sex/region of data. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared residual between the fitted line and right term. Finally, only specify \code{eOpen} when working with a single region/sex/period of data, otherwise the same value will be passed in irrespective of mortality and sex.
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$age}, \code{$sex}, and \code{$cod} (if there are more than 1 region/sex/intercensal period).
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param exact.ages optional. A user-specified vector of exact ages to use for coverage estimation
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. Rows indicate data partitions, as indicated by the optional \code{$cod} variable.
#' 
#' @export

bh2 <- function(x, minA = 15, maxA = 75, minAges = 8, sex = "f", exact.ages = NULL, eOpen = NULL){

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
	
	tab1        <- split(tab, tab$cod)
	
	
	coverages <- as.data.frame(
			do.call(
					rbind,
					lapply(
						tab1, 
						bh2coverageFromYear, 
						minA = minA,  
						maxA = maxA,
						minAges = minAges,  	
						exact.ages = exact.ages,
						eOpen = eOpen
                  )))
	coverages
}
