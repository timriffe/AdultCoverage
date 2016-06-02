
# Author: tim
###############################################################################
# contains functions related to Bennet-Horiuchi methods

bh1CoverageFromAges <- function(codi, agesFit, minA = 15, maxA = 75){
	
	if (!"Cx" %in% colnames(codi)){
		codi <- bh1MakeColumns(codi = codi, minA = minA, maxA = maxA)
	}
	inds    <- codi$age %in% agesFit
	sum(codi$Cx[inds]) / length(agesFit)
}

bh1MakeColumns <- function(codi, minA = 15, maxA = 75){
	
	# this throws an error if sex isn't coded as expected
	sex                    <- detectSex(Dat = codi, sexColumn = "sex")
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
	
	###################################################################
	# TR: begin precarious chunk
	# optional specify eOpen?
	ratio                  <- with(codi,
			                   sum(deathLT[ages %in% c(10:39)]) / 
					           sum(deathLT[ages %in% c(40:59)]))
	
	if (sex == "f"){
		# TODO: expand ex in-ine out to actual open ages.
		# model lifetable
		# based on Bennett & Horiuchi (1984) 
		# "Mortality Estimateion fro Registered Deaths in Less Developed Countries", Demography
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
	###################################################################
    # TR: end precarious chunk
	###################################################################

	
	###############
	# change to within() statement? This is rather convoluted
#    minus          <- ((eOpen * codi$growth[N])^2) / 6


	eON   <- eOpen * codi$growth[N]

	codi <- within(codi, {
				pop_a <- 0
				pop_a[N] <- deaths[N] * exp(eON - (eON ^ (1/3)))
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

bh1CoverageFromYear <-  function(codi, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL){        ##  Data
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
								maxA = maxA )
	
	coverage <- bh1CoverageFromAges(codi = codi, agesFit = agesFit)

	data.frame(cod = unique(codi$cod), coverage = coverage, lower = min(agesFit), upper = max(agesFit))
}

# TODO: detect sex rather than specify as argument X <-x
bh1 <- function(X, minA = 15, maxA = 75, minAges = 8, exact.ages = NULL){
	
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
	
	tab1        <- split(tab, X$cod)
	
	coverages <- as.data.frame(
					do.call(
						rbind,
						lapply(
							tab1, 
							bh1CoverageFromYear, 
							minA = minA, 
							maxA = maxA,
							minAges = minAges,  
							exact.ages = exact.ages
            )))
	#return(data.frame(Coverage = coverages,correctionFactor = 1/coverages))
	
	coverages
}

###################################################################################
# End BH1 functions
###################################################################################
# Start BH2 functions
###################################################################################

bh2CoverageFromAges <- function(codi, agesFit){
	inds    <- codi$age %in% agesFit
	sum(codi$Cx[inds]) / length(agesFit)
}

# change name to GB
bh2MakeColumns <- function(codi, minA = 15, maxA = 75,  minAges = 8, sex, agesFit){
	
	AgeInt      <- detectAgeInterval(Dat = codi, MinAge =  minA, MaxAge = maxA, ageColumn = "age")
	ages        <- codi$age
	dif         <- yint2(X = codi)
	agesi       <- codi$age %in% agesFit
	
	# just get left term / right term
	slope       <- with(codi, 
			sd(lefterm[age %in% agesFit]) / sd(rightterm[age %in% agesFit])
	)
	intercept   <- with(codi, 
			(mean(lefterm[age %in% agesFit]) - 
						slope * mean(rightterm[age %in% agesFit]))
	) 
	
	relcomp     <- exp(intercept * dif)
	# relcomp <- .96
	# adjust the first population count
	codi$pop1adj <- codi$pop1 / relcomp
	
	codi$birthdays            <- 0
    # iterate over age groups >= 10
	
	for (j in seq_along(ages)[ages >= minA]) {
		# take geometric average of p1 pop vs p2 pop within same cohort
		codi$birthdays[j]       <- 
				1 / AgeInt * sqrt(codi$pop1adj[j - 1] * codi$pop2[j])
	} # end age loop
	
    # age-specific growth
	
	codi[["growth"]]	   <-  log(codi$pop2 / codi$pop1adj) / dif
	codi$growth[is.infinite(codi$growth)] <- 0
	
	codi$cumgrowth         <-  0
	codi$cumgrowth[1]      <- 2.5 * codi$growth[1]
	
	for (j in 2:length(ages)){                    # TR: why 5 times?
		codi$cumgrowth[j]  <- 2.5 * codi$growth[j] + 5 * sum(codi$growth[(j - 1):1])
	}
	
	codi$death_tab         <- codi$death * exp(codi$cumgrowth)
	
	ratio                  <- sum(codi$death_tab[ages %in% c(10:39)]) / sum(codi$death_tab[ages%in%c(40:59)])
	
	if (sex == "f"){
		# TODO: expand ex in-ine out to actual open ages.
		# model lifetable
		# based on Bennett & Horiuchi (1984) 
		# "Mortality Estimation from Registered Deaths in Less Developed Countries", Demography
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
	
	# TODO: modify definition of AllLevels to be more robust.
	AllLevels             <- 3:25
	CDlevel               <- splinefun(AllLevels~standardratios)(ratio)
	# open at age 110 ....
	OA                    <- max(ages)
	availages             <- as.integer(colnames(ex))
	
	stopifnot(OA %in% availages)
	
	eOpen                 <- splinefun(ex[,as.character(OA)]~1:25)(CDlevel)
	
	N                     <- nrow(codi)
	codi$growth[is.nan(codi$growth)] <- 0
	
	minus                 <- ((eOpen * codi$growth[N])^2) / 6
	codi$pop_a            <- codi$death[N] * (exp(eOpen * codi$growth[N]) - minus)
	
	# TODO: remove loop?
	for(j in N:1){
		codi$pop_a[j - 1] <- codi$pop_a[j] * exp(AgeInt * codi$growth[j - 1]) + 
				codi$death[j - 1] * exp(AgeInt / 2 * codi$growth[j - 1])
	}
	
	
	codi$Cx               <-  codi$pop_a / codi$birthdays
	codi
}

bh2coverageFromYear <- function(codi, minA = 15, maxA = 75, minAges = 8, sex = "f", exact.ages = NULL){
	codiggb      <- ggbMakeColumns(codi = codi, minA = minA, maxA = maxA)
	# this is a test
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
								minAges = minAges, 
								sex = sex, 
								agesFit = agesFit)
	bh2CoverageFromAges(codi = codi, agesFit = agesFit )
}



#######################################
# TODO: death column should be mean intercensal deaths within an age. Could be simple avg,
# or the avg of the death around the time of the first and second censuses. Would be nice
# to calculate it from a separate Deaths object. Rather than have everything in 'x',
# the functions should have a Pop data.frame and a Deaths data.frame, fully specified
# in terms of AgeIntervals and YearIntervals, etc.
#######################################

#'
#' @title estimate the coverage coefficients
#' 
#' @description x
#' 
#' @param x
#' 
#' @return x
#' 
#' @export 
#' 
#' 

bh2 <- function(x, minA = 10, maxA = 75, minAges = 8, sex = "f", exact.ages = NULL){
	tab        <- data.frame(x)           ##  Dat in data frame : cod, age, pop1, year1, pop2, year2, death
	tab$pop1   <- as.double(tab$pop1)
	tab$pop2   <- as.double(tab$pop2)
	tab$death  <- as.double(tab$death)
	tab1       <- split(tab,x$cod)
	
	ages       <- sort(unique(tab$age))

	coverages  <- unlist(lapply(
					tab1, 
					bh2coverageFromYear, 
					minA = minA,  
					maxA = maxA,
					minAges = minAges,  
					sex = sex,
					exact.ages = exact.ages))
	coverages
}
