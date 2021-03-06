
# Author: tim
###############################################################################
# contains functions related to Bennett-Horiuchi methods


#' @title estimate remaining life expectancy in the open age group
#' @description This calculation is based on an indirect method to reference the Coale-Demeny West model life table. First one makes a pseudo life table deaths column using some stable pop assumptions (different in SEG vs GGB-SEG). Then take the ratio of the sum of ages 10-39 to 40-59. These ratios have been worked out for each model life table level, so we can pick the level based on the ratio we produce from the data. From there, we just pick out the remaining life expectancy that corresponds to the top age in our data, which for now hopefully is not higher than 95. The model life tables do not go higher than 95 for now, but that's well beyond the range for this method family. If your data go beyond 85 or so, then just group down to 85, say, and estimate using that instead of keeping a high open age. Called by \code{segMakeColumns()} and \code{ggbsegMakeColumns()}, and not intended for direct user interface, because you need to produce the \code{$deathsLT} column. You can skip calling this function by specifying \code{eOpen} in the top call to \code{seg()} or \code{ggbseg()}.
#' 
#' @param codiaugmented the standard codi object being passed through the chain, but having been preprocessed in the course of \code{segMakeColumns()} or \code{ggbsegMakeColumns()}
#' 
#' @return numeric an estimate of remaining life expectancy in the open age group
#' 
#' @importFrom stats splinefun
#' @export

eOpenCD <- function(codiaugmented){
	# this needs to have the column deathsLT, as passed in 
	# by either segMakeColumns() or ggbsegMakeColumns(). These
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
	eOpen
}

#' @title select age range for SEG method automatically.
#' @description We try all possible age ranges given the constraints provided. Whichever age range has the lowest RMSE with respect to the corresponding coverage estimate is returned.
#' @inheritParams ggbgetAgesFit
#' @export
#' @return integer vector of ages that minimize the RMSE.

seggetAgesFit <- function (codi, 
                           minA = 15, 
                           maxA = 75, 
                           minAges = 8) {

  
  keep      <- codi$age >= minA & codi$age <= maxA
  maxAges   <- sum(keep)
  agesUniv  <- codi$age[keep]
  FirstAges <- agesUniv[agesUniv < 30]
  ind       <- 0
  agesL     <- list()
  for (Nr in maxAges:minAges) {
    its <- length(agesUniv) - Nr + 1
    for (set in 1:its) {
      ind <- ind + 1
      agesL[[ind]] <- agesUniv[set:(set + Nr - 1)]
    }
  }
  
  RMSE    <- lapply(agesL, seggetRMS, 
                    codi = codi) %>% 
             unlist()
  RMSEmin <- which.min(RMSE)
  agesfit <- agesL[[RMSEmin]]
  list(agesfit = agesfit, RMSE = RMSE[RMSEmin])
}

#' calculates RMSE for SEG and a specfied age vector.
#' @description return the RMSE of `Cx` with respect to the coverage estimate for a given range of ages.
#' @inheritParams ggbgetAgesFit
#' @param agesFit ineger vector of age trim for SEG method.
#' @export
#' @return numeric. RMSE

seggetRMS <- function(codi, agesFit){
  coverage      <- segCoverageFromAges(codi, agesFit = agesFit)
  Cx            <- codi$Cx
  keep          <- codi$age %in% agesFit  
  sqrt(mean((Cx[keep] - coverage["mean"]) ^ 2, na.rm = TRUE))
}

#' @title given a set of ages, what is the implied death registration coverage?
#' 
#' @description For a single year/sex/region of data (formatted as required by \code{seg()}, \code{ggbseg()}), what is the registration coverage implied by a given age range? Called by \code{segCoverageFromYear()} and  \code{ggbsegCoverageFromYear()}. Here, the function simply takes the arithmetic mean of a given age range of \code{$Cx}, as returned by \code{segMakeColumns()} or \code{ggbsegMakeColumns()}. Not intended for top-level use.
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{ggb()}
#' @param agesFit an integer vector of ages, either returned from \code{seggetAgesFit} or user-supplied.
#'
#' @return numeric. the estimated level of coverage.
#' @export

segCoverageFromAges <- function(codi, agesFit){
	stopifnot("Cx" %in% colnames(codi))
	inds    <- codi$age %in% agesFit
	
	Cx <- codi$Cx[inds]
	mn <- mean(Cx)
	qt <- quantile(Cx, c(.25,.5,.75),na.rm=TRUE)
	c(mean = mn,
	  median = qt[2],
	  l25 = qt[1],
	  u25 = qt[3],
	  weighted = sum(qt*c(.25,.5,.25)))
}


#' @title make the Bennett-Horiuchi quasi life table columns required by the estimation method
#' 
#' @description Called by \code{segCoverageFromYear()}. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. 
#' 
#' @inheritParams segCoverageFromYear
#' @param minAges.ggb the minimum number of adjacent ages to be used in estimating GGB coverage.
#'
#' @return codi, with many columns added, most importantly \code{$Cx}.
#' 
#' @export
#' @import dplyr
#' @import magrittr
#' @importFrom lubridate decimal_date
#' @importFrom DemoTools lt_id_L_T
#' @importFrom DemoTools age2int

segMakeColumns <- function(codi, 
                           minA = 15, 
                           maxA = 75, 
                           eOpen = NULL, 
                           nx.method = 2,
                           deaths.summed = FALSE, #
                           mig.summed = deaths.summed,
                           # New args to support delta method.
                           delta = FALSE,
                           exact.ages.ggb = NULL,
                           lm.method = "oldschool",
                           opt.method  = "r2",
                           minAges.ggb = 8){
  
  calc_bdaymean <- function(birthdays, AgeInt){
    mn <- rep(NA,length(birthdays))
    N <- length(birthdays)
    for (i in 1:(N-1)){
      mn[i] <- (birthdays[i] + birthdays[i+1]) *AgeInt[i] / 2
    }
    mn
  }
  
  # TR: this little thing is ugly for my taste, hoping it gets replaced.
  calc_pop_a <- function(deathsAvg,eON,AgeInt,growth){
    N        <- length(deathsAvg)
    pop_aa   <- rep(0,N)          
    # PJ fix, modified
    pop_aa[N] <-  deathsAvg[N] * ((exp(eON) - (eON ^ 2) / 6)) 
    #pop_aa[N] <- deathsAvg[N] * exp(eON) - ((eON ^ 2) ^ (1 / 6))
    for(j in N:2){
      pop_aa[j - 1] <- pop_aa[j] * exp( AgeInt[j-1] * growth[j - 1]) + 
        deathsAvg[j - 1] * exp( AgeInt[j-1] / 2 * growth[j - 1])
    }
    pop_aa
  }
  
  
  if (delta){
    ggb.res <- ggb(X = codi,
                   minA = minA,
                   maxA = maxA,
                   nx.method = nx.method,
                   exact.ages = exact.ages.ggb,
                   lm.method = lm.method,
                   deaths.summed = deaths.summed,
                   mig.summed = mig.summed,
                   opt.method = opt.method,
                   minAges = minAges.ggb)
    del <- ggb.res$a
  } else {
    del <- 0
  }
  
  # group inf if necessary
  codi                   <- group01(codi)
  codi                   <- reduceOpen(codi, maxA = 95, group = TRUE)
  N                      <- nrow(codi)
  
  if ((!"mig" %in% colnames(codi))){
    codi <- codi %>% mutate(mig = 0)
  }
  
  # TR: add migration arithmetic identical to in GGB changes.
  codi <-
    codi %>% 
    mutate(	date1          = ifelse(is.numeric(.data$date1), .data$date1, decimal_date(.data$date1)+1),
            date2          = ifelse(is.numeric(.data$date2), .data$date2, decimal_date(.data$date2)+1),
            AgeInt         = age2int(.data$age),
            dif            = .data$date2 - .data$date1,
            avg            = deaths.summed,
            deathsAvg      = ifelse(.data$avg, .data$deaths / .data$dif, .data$deaths),
            avg            = mig.summed,
            migAvg         = ifelse(.data$avg, .data$mig / .data$dif, .data$mig),
            pop1           = as.double(.data$pop1),
            pop2           = as.double(.data$pop2),
            mig            = as.double(.data$mig),
            #obs_rd         = sqrt(.data$pop1*.data$pop2) * .data$dif,
            # deathcum       = lt_id_L_T(.data$deathsAvg),
            # TR: follows RD column G calc
            growth         = (log(.data$pop2)- log(.data$pop1)) / .data$dif - 
                             .data$migAvg / sqrt(.data$pop2 * .data$pop1) + 
                              del,
            growth         = ifelse(is.infinite(.data$growth) | is.nan(.data$growth), 0, .data$growth),
            # cumulative growth, needed for estimating eOpen..
            cumgrowth      = .data$AgeInt * c(0,cumsum(.data$growth[ -N ])) + .data$AgeInt / 2 * .data$growth,
            deathLT        = .data$deathsAvg * exp(.data$cumgrowth))
  

	if (is.null(eOpen)){
		eOpen <- eOpenCD(codiaugmented = codi)
	} else {
		eOpen <- eOpen
	}
	
	eON   <- eOpen * codi$growth[N]

	codi <-
	  codi %>% 
	  mutate(
	    estNx = calc_pop_a(.data$deathsAvg,eON,.data$AgeInt,.data$growth) * .data$dif,
	    est5Nx = calc_bdaymean(.data$estNx, .data$AgeInt),
	    obs5Nx =  sqrt(.data$pop2 * .data$pop1) * .data$dif,
	    Cx = .data$est5Nx / .data$obs5Nx)


			
	################
	codi
}


#' @title estimate death registration coverage for a single year/sex/region using the Bennett-Horiuchi method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model life table, although the user may also supply a value. Called by \code{seg()}. Users probably do not need to use this function directly.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. 
#' 
#' @param codi \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$sex}, \code{$age}, and \code{$id} (to indicate regions, periods, sexes).
#' @inheritParams seg
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. 
#' 
#' @importFrom stats quantile
#' 
#' @export


segCoverageFromYear <-  function(
                 codi, 
								 minA = 15, 
								 maxA = 75, 
								 minAges = 8, 
								 exact.ages = NULL, 
								 eOpen = NULL, 
								 nx.method = 2,
								 deaths.summed = FALSE,
								 mig.summed = deaths.summed,
								 # New args to support delta method.
								 delta = FALSE,
								 exact.ages.ggb = NULL,
								 lm.method = "oldschool",
								 opt.method = "r2"){       
  
  
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

	codi     <- segMakeColumns(	
	              codi = codi, 
								minA = minA, 
								maxA = maxA,
								eOpen = eOpen,
								nx.method = nx.method,
								deaths.summed = deaths.summed,
								mig.summed = mig.summed,
								# New args to support delta method.
								delta = delta,
								exact.ages.ggb = exact.ages.ggb,
								lm.method = lm.method,
								opt.method = opt.method,
								minAges.ggb = minAges)
	
	# TR: 09.09.2020 switched order, for whatever reason agesFit was
	# previously determined before segMakeColumns().
	# TR: 19.06.2020 pick ages by minimizing RMS
	# assuming a horizontal coverage line on 
	# a sequential range of ages of Cx
	if (is.null(exact.ages)){
	  fit.res <- seggetAgesFit(codi = codi, 
	                           minA = minA, 
	                           maxA = maxA,
	                           minAges = minAges)
	  agesFit <- fit.res$agesfit
	} else {
	  agesFit <- exact.ages
	}
	
	coverage <- segCoverageFromAges(codi = codi, agesFit = agesFit)

	meta <- 	data.frame(
	       id = unique(codi$id), 
			   Mxcoverage = coverage[1], 
			   median = coverage[2],
			   l25 = coverage[3],
			   u25 = coverage[4],
			   weighted = coverage["weighted"],
			   lower = min(agesFit), 
			   upper = max(agesFit)
			   )
	if (is.null(exact.ages)){
	  meta$RMSE <- fit.res$RMSE
	}
	out <- list(coverages = meta, codi = codi)
	out
}


#' @title estimate death registration coverage using the synthetic extinct generation method 
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the ages to fit against are not specified, then these are optimized. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model life table, although the user may also supply a value.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} columns are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$id} with a unique values for each data chunk. Different values of \code{$id} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths in each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range in the traditional visual way, see \code{plot.ggb()}, when working with a single year/sex/region of data. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared residual between the fitted line and right term. Finally, only specify \code{eOpen} when working with a single region/sex/period of data, otherwise the same value will be passed in irrespective of mortality and sex.
#' 
#' If \code{exact.ages} is specified as \code{NULL}, coverage is estimated by minimizing the RMSE of the coverage estimate versus \code{$Cx}.
#' 
#' @inheritParams ggb
#' @param delta logical. Do we perform the so-called delta adjustment?
#' @param exact.ages.ggb optional vector of ages used to estimate GGB coverage (if \code{delta} is \code{TRUE})
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' @return a \code{data.frame} with columns for the coverage coefficient \code{$coverage}, and the minimum \code{$lower} and maximum \code{$upper} of the age range on which it is based. Rows indicate data partitions, as indicated by the optional \code{$id} variable. \code{$l25} (\code{$u25}) give the mean of the lower (upper) quartile of the distribution of age-specific coverage estimates.
#' 
#' @export
#' 
#' @references Bennett Neil G, Shiro Horiuchi. Estimating the completeness of death registration in a closed population. Population Index. 1981; 1:207-221.
#' 
#' Preston, S. H., Coale, A. J., Trussel, J. & Maxine, W. Estimating the completeness of reporting of adult deaths in populations that are approximately stable.  Population Studies, 1980; v.4: 179-202
#' 
#' @examples 
#' # The Mozambique data
#' res <- seg(Moz)
#' res
#' # The Brasil data
#' BM <- seg(BrasilMales)
#' BF <- seg(BrasilFemales)
#' head(BM)
#' head(BF)

seg <- function(X, 
				minA = 15, 
				maxA = 75, 
				minAges = 8, 
				exact.ages = NULL, 
				eOpen = NULL, 
				nx.method = 2,
				deaths.summed = FALSE,
				mig.summed = deaths.summed,
				delta = FALSE,
				exact.ages.ggb = NULL,
				lm.method = "oldschool",
				opt.method = "r2"){
	
	tab1        <- headerPrep(X)
	coverages <- as.data.frame(
					do.call(
						rbind,
						lapply(
							tab1,
							function(Y){
								segCoverageFromYear(
								    Y,
								    minA = minA, 
										maxA = maxA,
										minAges = minAges,  
										exact.ages = exact.ages,
										eOpen = eOpen,
										nx.method = nx.method,
										deaths.summed = deaths.summed,
										mig.summed = mig.summed,
										delta = delta,
										exact.ages.ggb = exact.ages.ggb,
										lm.method = lm.method,
										opt.method = 	opt.method
										)$coverages	
							}
            )))
	
	coverages
}

#' plot the age-pattern of coverage estimates
#' 
#' @description the SEG method works by averaging the coverage estimates over a range of ages. 
#' Users may wish to see the age pattern for diagnostic purposes.
#' 
#' @details All arguments are essentially the same as those given to \code{seg()}
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$age}, \code{$sex}, and \code{$id} (if there are more than 1 region/sex/intercensal period).
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param exact.ages optional. A user-specified vector of exact ages to use for coverage estimation
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' @param deaths.summed logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param log logical. should we log the y axis?
#' @param nx.method integer 2 or 4. 4 is smoother
#' 
#' @return Function called for its graphical side effects
#' 
#' @importFrom grDevices gray
#' @importFrom graphics abline legend locator mtext par plot points rect segments text
#' @export 
#' 
#' @examples 
#' \dontrun{
#' segplot(Moz)
#' }


segplot <- function(
		X, 
		minA = 15, 
		maxA = 75, 
		minAges = 8, 
		exact.ages = NULL, 
		eOpen = NULL, 
		deaths.summed = FALSE,
		log = FALSE,
		nx.method = 2){
	tab1        <- headerPrep(X)
	
	if (length(tab1) > 1){
		warning("codi was not unique, taking first element")
	}
	tab1 <- tab1[[1]]
	goods <- segCoverageFromYear(
	                        tab1,
		                      minA = minA, 
		                      maxA = maxA,
		                      minAges = minAges,  
		                      exact.ages = exact.ages,
		                      eOpen = eOpen,
		                      deaths.summed = deaths.summed,
		                      nx.method = nx.method)
	codi     <- goods$codi
	coverage <- goods$coverages
	keep 	   <- codi$age >= minA &  codi$age <= maxA
	ages 	   <- codi$age[keep]
	yvals 	 <- codi$Cx[keep]
	keep2 	 <- codi$age >= coverage$lower & codi$age <= coverage$upper
	ages2 	 <- codi$age[keep2]
	yvals2 	 <- codi$Cx[keep2]
	
	# TR: 29 May ylim changes
	if (!log){
		ymin <- min(c(0, min(yvals)))
		ymax <- max(1, max(yvals))
		ylim <- c(ymin, ymax)
	} else {
		ylim <- range(yvals)
	}

	
	plot(ages, yvals, pch = 16, ylim = ylim, col = gray(.5), ylab="Cx", log = ifelse(log,"y",""))
	rect(ages2[1]-2.5,coverage$l25,ages2[length(ages2)]+2.5,coverage$u25, 
			border = NA, col = "#00000040")
	points(ages2, yvals2, pch = 16, col = "black")
	segments(ages2[1]-2.5,coverage$coverage,ages2[length(ages2)]+2.5,coverage$coverage,lwd=2)
	
	# TODO let user interactively select points similar to ggbplot, for purposes of
	# refitting, etc. This would be a second age-range diagnostic.
}


###################################################################################
# End seg functions
###################################################################################
# Start ggbseg functions
###################################################################################


#' @title make the Bennett-Horiuchi quasi life table columns required by the second estimation method
#' 
#' @description Called by \code{ggbsegCoverageFromYear()}. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. 
#' @details \code{agesFit} is a vector passed in from \code{ggbsegMakeColumns()}, and it was either estimated using the GGB automatic method there, or simply came from the argument \code{exact.ages} specified in \code{ggbseg()}. By default we just automatically estimate these. \code{eOpen} can be either user-specified, or it will be estimated automatically using \code{eOpenCD()}.
#' 
#' @inheritParams ggbseg
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by \code{ggbseg()}
#' @param agesFit.ggb vector of ages for calculating coverage in the GGB method.
#' @return codi, with many columns added, most importantly \code{$Cx}.
#' 
#' @export
#' @import dplyr
#' @import magrittr

ggbsegMakeColumns <- function(
		codi, 
		minA = 15, 
		maxA = 75,  
		agesFit.ggb, 
		eOpen = NULL,
		deaths.summed = FALSE,
		mig.summed = deaths.summed,
		lm.method = "oldschool",
		nx.method = 2){
  
  codi                   <- group01(codi)
  codi                   <- reduceOpen(codi, maxA = 95, group = TRUE)

  codi1                  <- ggbMakeColumns(codi, 
                                           minA = minA, 
                                           maxA = maxA,	
                                           nx.method = nx.method,
                                           deaths.summed = deaths.summed,
                                           mig.summed = mig.summed)
                                           
  
	dif          <- codi1$dif[1]
	
	# just get left term / right term
	# slopeint() is in the ggb() family
	# TR: this function will likely change.
	ab           <- slopeint(codi = codi1, 
	                         agesfit = agesFit.ggb, 
	                         lm.method = lm.method)
	
	# the only difference between this method and seg is that
	# in the next couple lines we use pop1adj instead of pop1...
	relcomp      <- exp(ab$a * dif)
   	
  codi <-
    codi %>% 
    # just call it pop1 still, used to be pop1adj, but this lets us recycle code easier.
    mutate(
      pop1in = .data$pop1,
      pop1 = .data$pop1in / relcomp) %>% 
    segMakeColumns(minA = minA,
                   maxA = maxA, 
                   eOpen = eOpen, 
                   deaths.summed = deaths.summed,
                   mig.summed = mig.summed,
                   nx.method = nx.method,
                   delta = FALSE) %>% 
    # rename back to original intended names
    dplyr::rename(pop1adj = .data$pop1, pop1 = .data$pop1in)
  # TR: add in migraion cals to SEG then activate this.
  #,
  #                 mig = mig)
  #,
   #                delta = FALSE)
	
	codi
}


#' @title estimate death registration coverage for a single year/sex/region using the modified Bennett-Horiuchi method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model life table, although the user may also supply a value. The difference between this method and \code{seg()} is that here we adjust census 1 part way through processing, based on some calculations similar to GGB. Called by \code{ggbseg()}. Users probably do not need to use this function directly.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. 
#' 
#' @inheritParams ggbseg
#' @param codi a chunk of data (single sex, year, region, etc) with all columns created by \code{ggbsegMakeColumns()}
#' @return a \code{data.frame} with columns for the coverage coefficient, and the min and max of the age range on which it is based. 
#' 
#' @export

ggbsegCoverageFromYear <- function(codi, 
								minA = 15, 
								maxA = 75, 
								minAges = 8, 
								exact.ages.ggb = NULL, 
								exact.ages.seg = NULL,
								eOpen = NULL, 
								deaths.summed = FALSE,
								mig.summed = deaths.summed,
								lm.method = "oldschool",
								nx.method = 2,
								opt.method = "r2"){
	
	# Get age range using the GGB auto fitting
	if (is.null(exact.ages.ggb)){
	  codiggb      <- ggbMakeColumns(codi = codi, 
	                                 minA = minA, 
	                                 maxA = maxA, 
	                                 nx.method = nx.method,
	                                 deaths.summed = deaths.summed,
	                                 mig.summed = mig.summed)
		fit.res.ggb <- ggbgetAgesFit(
		                codi = codiggb, 
								    minA = minA, 
							    	maxA = maxA, 
								    minAges = minAges,
								    lm.method = lm.method,
								    opt.method =opt.method )
		agesFit.ggb <- fit.res.ggb$agesfit
	} else {
		agesFit.ggb <- 	exact.ages.ggb
	}
	# Get age range using the SEG auto fitting

	codi         <- ggbsegMakeColumns(
								    codi = codi, 
								    minA = minA, 
								    maxA = maxA, 	
								    agesFit.ggb = agesFit.ggb,
								    eOpen = eOpen,
								    deaths.summed = deaths.summed,
								    mig.summed = mig.summed,
								    nx.method = nx.method)

	# TR: this moved to after column creation.
	if (is.null(exact.ages.seg)){
	  fit.res.seg <- seggetAgesFit(codi = codi, 
	                               minA = minA, 
	                               maxA = maxA, 
	                               minAges = minAges)
	  
	  agesFit.ggbseg <-  fit.res.seg$agesfit
	} else {
	  agesFit.ggbseg <- exact.ages.seg
	}
	
	coverage     <- segCoverageFromAges(codi = codi, agesFit = agesFit.ggbseg)
	
	res <- data.frame( 
	           id = unique(codi$id), 
	           Mxcoverage = coverage[1], 
	           median = coverage[2],
	           l25 = coverage[3],
	           u25 = coverage[4],
	           weighted = coverage["weighted"], 
	           lower.ggb = min(agesFit.ggb), 
	           upper.ggb = max(agesFit.ggb),
	           lower.ggbseg = min(agesFit.ggbseg), 
	           upper.ggbseg = max(agesFit.ggbseg))
	
	if (is.null(exact.ages.ggb)){
	  res$RMSE.ggb    <- fit.res.ggb$RMSE
	}
	if (is.null(exact.ages.seg)){
	  res$RMSE.ggbseg <-  fit.res.seg$RMSE
	}
	
	res
}


#' @title estimate death registration coverage using the hybrid generalized growth balance and synthetic extinct generation
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model life table, although the user may also supply a value. The difference between this method and \code{seg()} is that here we adjust census 1 part way through processing, based on some calculations similar to GGB.
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} columns are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$id} with a unique values for each data chunk. Different values of \code{$id} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths in each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range for the GGB method in the traditional visual way, see \code{ggbChooseAges()}, when working with a single year/sex/region of data. A similar visual picker for SEG is not yet implemented. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared error (RMSE) between the fitted line and right term for the GGB (first stage), and the SEG stage minimizes RMSE for \code{$Cx} against the coverage estimate implied by the age range. Finally, only specify \code{eOpen} when working with a single region/sex/period of data, otherwise the same value will be passed in irrespective of mortality and sex.
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$age}, \code{$sex}, and \code{$id} (if there are more than 1 region/sex/intercensal period).
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param exact.ages.ggb optional. A user-specified vector of exact ages to use for coverage estimation in the GGB (first stage) part of the estimation.
#' @param exact.ages.seg optional. A user-specified vector of exact ages to use for coverage estimation in the SEG (second stage) part of the estimation.
#' @param eOpen optional. A user-specified value for remaining life-expectancy in the open age group.
#' @param deaths.summed logical. Is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param mig.summed logical. Is the (optional) net migration column \code{mig} given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param lm.method character, one of:\itemize{
#'   \item{\code{"oldschool"}} default sd ratio operation of still unknown origin
#'   \item{\code{"lm"} or \code{"ols"}} for a simple linear model
#'   \item{\code{"tls"}, \code{"orthogonal"}, or \code{"deming"}} for total least squares
#'   \item{\code{"tukey"}, \code{"resistant"}, or "\code{"median"}} for Tukey's resistant line method
#' }
#' @param nx.method integer. either 2 or 4. 4 is smoother.
#' @param opt.method What should we try to minimize when picking ggb age trims? Current options are `"RMSE"`, `"ORSS"`, `"MAE"`, `"MAPE"`, or `"r2"`. Default `"r2"`.
#' 
#' @return a \code{data.frame} with columns for the coverage coefficient \code{$coverage}, and the minimum \code{$lower} and maximum \code{$upper} of the age range on which it is based. Rows indicate data partitions, as indicated by the optional \code{$id} variable.
#' 
#' @export
#' @references Hill K. Methods for measuring adult mortality in developing countries: a comparative review. The global burden of disease 2000 in aging populations. Research paper; No. 01.13; 2001.
#' 
#' Hill K, You D, Choi Y.  Death distribution methods for estimating adult mortality: sensitivity analysis with simulated data errors. Demographic Research. 2009; 21:235-254.
#' 
#' Preston, S. H., Coale, A. J., Trussel, J. & Maxine, W. Estimating the completeness of reporting of adult deaths in populations that are approximately stable.  Population Studies, 1980; v.4: 179-202
#' 
#' 
#' @examples 
#' # The Mozambique data
#' res <- ggbseg(Moz)
#' res
#' # The Brasil data
#' BM <- ggbseg(BrasilMales)
#' BF <- ggbseg(BrasilFemales)
#' head(BM)
#' head(BF)
ggbseg <- function(X, 
				minA = 15, 
				maxA = 75, 
				minAges = 8, 
				exact.ages.ggb = NULL, 
				exact.ages.seg = NULL,
				eOpen = NULL, 
				deaths.summed = FALSE,
				mig.summed = deaths.summed,
				nx.method = 2,
				lm.method = "oldschool",
				opt.method = "r2"){
	# TR: modularized Apr 2, 2017
    tab1        <- headerPrep(X)
	
	coverages <- as.data.frame(
			do.call(
					rbind,
					lapply(
						tab1, 
						ggbsegCoverageFromYear, 
						minA = minA,  
						maxA = maxA,
						minAges = minAges,  	
						exact.ages.ggb = exact.ages.ggb,
						exact.ages.seg = exact.ages.seg,
						eOpen = eOpen,
						deaths.summed = deaths.summed,
						mig.summed = mig.summed,
						lm.method = lm.method,
						nx.method = nx.method,
						opt.method = opt.method
                  )))
	coverages
}
