
# Author: tim
###############################################################################
# contains functions related to the generalized growth balance method


#' @title calculate the root means square of the error to help find optimal age range
#' 
#' @description Called by \code{ggbgetAgesFit()} whenever the user does not want to manually determine the age range used to determine registration coverage. Probably no need to be called by top-level users. If a user would rather determine the optimal age range some other way, then look to \code{ggbcoverageFromYear()} where \code{ggbgetRMS} is called and add another condition or make it call something else.
#' 
#' @param agesi the set of ages used for this iteration
#' @param codi \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, and \code{$age}. 
#' @inheritParams slopeint
#' 
#' @return the RMSE
#' 
#' @export 
ggbgetRMS <- function(agesi, codi, lm.method = "oldschool", opt.method = "RMS"){

	codi <- ggbFittedFromAges(codi, agesfit = agesi, lm.method = lm.method) %>% 
	  filter(.data$age %in% agesi)

	if (opt.method == "RMS"){
	  out <- sqrt(mean((codi$leftterm - codi$fitted)^2) )
	}
	if (opt.method == "logRMS"){
	  out <- sqrt(mean((log(codi$leftterm) - log(codi$fitted))^2) )
	}
	if (opt.method == "ORSS"){
	  out <- prcomp(cbind(codi$fitted,codi$leftterm))$sdev[1]
	}
	if (opt.method == "ORSS"){
	  out <- prcomp(cbind(log(codi$fitted),log(codi$leftterm)))$sdev[1]
	}
  out
}




#' @title estimate death registration coverage for a single year/sex/region using the GGB method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method is based on finding a best-fitting linear relationship between two modeled parameters (right term and left term), but the fit, and resulting coverage estimate, depend on exactly which age range is taken. This function either finds a nice age range for you automatically, or you can specify an exact vector of ages. Called by \code{ggb()}. Users probably don't need to call this directly. Just use \code{ggb()} instead.
#' 
#' Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. 
#' 
#' @inheritParams ggb
#' @param codi a chunk of data from a single \code{id}
#' @export

ggbcoverageFromYear <- function(codi, 
								exact.ages = NULL, 
								minA = 15, 
								maxA = 75, 
								minAges = 8, 
								deaths.summed = FALSE,
								mig.summed = deaths.summed,
								lm.method = "oldschool",
								nx.method = 2
								){
	
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
	
	# TR: test add this step, just in case
	codi    <- codi[ order(codi$age), ]
	
  if (!"leftterm" %in% colnames(codi)){
	   codi <- ggbMakeColumns(codi = codi, 
	                          minA = minA, 
	                          maxA = maxA, 
	                          deaths.summed = deaths.summed,
	                          mig.summed = mig.summed,
	                          nx.method = nx.method)
  }
	
	if (!is.null(exact.ages) & length(exact.ages) >= 3){
		agesfit <- exact.ages
	} else {                          
		fit.res <- ggbgetAgesFit(codi = codi, 
								 minA = minA, 
								 maxA = maxA, 
								 minAges = minAges,
								 lm.method = lm.method)
		agesfit <- fit.res$agesfit
	}
	
	# TR: added 17 June, 2016. Get Lambda to adjust first census:
	# TR: 3-4-2017 no longer optional

	coefs    <- slopeint(codi, 
	                     agesfit = agesfit, 
	                     lm.method = lm.method)
	# (time between censuses)
	dif      <- codi$dif[1]
	# TR: 3-4-2017 this is k1/k2
	# Relative completeness of Census 1 to Census 2
	delta    <- exp(coefs$a * dif)
	# TR: 2-9-2020
	
  if (delta > 1){
    k1 <- 1
    k2 <- 1 / delta
  }
	if (delta == 1){
	  k1 <- 1
	  k2 <- 1
	}
	if (delta < 1){
	  k1 <- delta
	  k2 <- 1
	}
	# now get everything from the k parameters 
#	a        <- log(delta) / dif
	coverage <- sqrt(k1 * k2) / coefs$b

	result   <- data.frame(
	       id = unique(codi$id), 
			   Mxcoverage = coverage, 
			   lower = min(agesfit), 
			   upper = max(agesfit),
			   a = coefs$a, 
			   b = coefs$b, 
			   delta = delta, 
			   k1 = k1, 
			   k2 = k2,
			   k3 = 1 / coefs$b,
			   t1 = codi$date1[1],
			   t2 = codi$date2[2],
			   t = dif,
			   lm.method = lm.method,
			   nx.method = nx.method)
	if (is.null(exact.ages)){
	  result$RMSE = fit.res$RMSE
	}
	# TR: Add r2 and RMSE to output
   
	result
}



#' @title make the growth-adjusted quasi life table columns required by GGB method
#' 
#' @description Called by `ggbChooseAges()` and `ggbcoverageFromYear()`. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. If columns produced by `ggbMakeColumns()` are not present, then we call it here just to keep things from breaking.
#' 
#' @inheritParams slopeint
#' @return codi, with the column `fitted` added.
#' 
#' @export

ggbFittedFromAges <- function(codi, 
                              agesfit, 
                              lm.method = "oldschool"){
	
	coefs       <- slopeint(codi, 
	                        agesfit, 
	                        lm.method = lm.method)
	codi$fitted <- coefs$a + codi$rightterm * coefs$b 
	codi
}

#' @title make the growth-adjusted quasi life table columns required by GGB method
#' 
#' @description Called by `ggbChooseAges()` and `ggbcoverageFromYear()`. This simply modulates some code that would otherwise be repeated. Users probably don't need to call this function directly. 
#' 
#' @param codi a chunk of data (single sex, year, region, etc) with all columns required by `ggb()`
#' @inheritParams ggb
#' 
#' @return codi, with many columns added, most importantly `rightterm`, `leftterm`, and `keep`.
#' 
#' @export
#' @importFrom DemoTools age2int 
#' @importFrom DemoTools lt_id_L_T 
#' @importFrom lubridate decimal_date
#' @import dplyr
#' @import magrittr

ggbMakeColumns <- function(codi, 
                           minA = 15, 
                           maxA = 75, 
                           nx.method = 2,
                           deaths.summed = FALSE, #
                           mig.summed = deaths.summed
                           ){           # maybe add a flag for whether summed or avg
  # group inf if necessary
  codi                   <- group01(codi)
  N                      <- nrow(codi)
  
  
  if ((!"mig" %in% colnames(codi))){
    codi <- codi %>% mutate(mig = 0)
  }
  
  codi <-
    codi %>% 
    mutate(	date1          = ifelse(is.numeric(.data$date1), .data$date1, decimal_date(.data$date1)),
            date2          = ifelse(is.numeric(.data$date2), .data$date2, decimal_date(.data$date2)),
            AgeInt         = age2int(.data$age),
            dif            = .data$date2 - .data$date1,
            avg            = deaths.summed,
            deathsAvg      = ifelse(.data$avg, .data$deaths / .data$dif, .data$deaths),
            pop1           = as.double(.data$pop1),
            pop2           = as.double(.data$pop2),
            mig            = as.double(.data$mig),
            avg            = mig.summed,
            migAvg         = ifelse(.data$avg, .data$mig / .data$dif, .data$mig),
            # RD: take geometric mean of pop1 and pop2
            # refer to RD spreadsheet 
            pop1cum        = lt_id_L_T(.data$pop1),
            pop2cum        = lt_id_L_T(.data$pop2) ,
            deathcum       = lt_id_L_T(.data$deathsAvg),
            migcum         = lt_id_L_T(.data$migAvg),               # new
            migcum2        = rev(cumsum(rev(.data$migAvg * .data$dif))),
            # Nx
            birthdays      = est_birthdays(pop1 = .data$pop1, pop2 = .data$pop2, 
                                           AgeInt = .data$AgeInt, nx.method = nx.method) * .data$dif,
            Tx             = sqrt(.data$pop1cum * .data$pop2cum) * .data$dif, # PYL in RD spreadsheet
            # this follows the SEG way of incorporating mig info, guarantees
          
            #cumgrowth      = (.data$pop2cum - .data$pop1cum - .data$migcum) / .data$Tx,
            #rightterm      = .data$dif * .data$deathcum / .data$Tx,   # from RD spreadsheet (X)
           # leftterm       = .data$bxp - .data$cumgrowth, # from RD spreadsheet (Y)
            
            Dxp = rev(cumsum(rev(.data$deathsAvg * .data$dif))),
            Mxp = rev(cumsum(rev(.data$migAvg * .data$dif))),
            PYL = .data$dif * sqrt(.data$pop1cum * .data$pop2cum),
            Nx = c(0,sqrt(.data$pop1[-N] * .data$pop2[-1]))/.data$AgeInt * .data$dif,
            bxp = .data$Nx / .data$PYL,
            rxp_m_ixp = (.data$pop2cum - .data$pop1cum - .data$migcum2) / .data$PYL,
            rightterm = .data$Dxp / .data$PYL,
            leftterm = .data$bxp - .data$rxp_m_ixp,
            # growth         = log(.data$pop2 / .data$pop1) / .data$dif - .data$migAvg / .data$Tx / .data$dif,
            # growth         = ifelse(is.infinite(.data$growth) | is.nan(.data$growth), 0, .data$growth),
            #cumgrowth      = .data$AgeInt * c(0,cumsum(.data$growth[ -N ])) + .data$AgeInt / 2 * .data$growth,
            #cumgrowth_old  = log(.data$pop2cum / .data$pop1cum) / .data$dif,
            #growth_old     = c(0,diff(.data$cumgrowth_old)),
            #cumgrowth_bq   =  (.data$pop2cum - .data$pop1cum - .data$migcum) / .data$Tx,
           # same
            #rightterm      = .data$deathcum / .data$Tx,
            #rightterm_test = .data$rightterm / .data$dif,
            #leftterm       = (.data$birthdays / .data$Tx) - .data$cumgrowth,
            #leftterm_old   = (.data$birthdays / .data$Tx) - .data$cumgrowth_old,
            #leftterm_bq   = (.data$birthdays / .data$Tx) - .data$cumgrowth_bq,

            keep           = .data$Tx != 0 & .data$birthdays != 0 & .data$age >= minA & .data$age <= maxA)
  
  codi
}



#' @title determine the age range that minimizes the mean squared error
#' @description Called by \code{ggbcoverageFromYear()} whenever \code{exact.ages} are not given. This automates what one typically does visually.
#' @seealso code{\link{ggbChooseAges}}
#' 
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @inheritParams slopeint
#' 
#' @return a vector of ages that minimizes the RMSE
#' 
#' @export

ggbgetAgesFit <- function(codi, 
                          minA = 15, 
                          maxA = 75, 
                          minAges = 8,
                          lm.method = "oldschool"){
	
	maxAges   <- sum(codi$keep)
	agesUniv  <- codi$age[codi$keep]
	
	FirstAges <- agesUniv[agesUniv < 30]
	
	ind       <- 0
	agesL     <- list()
	# determine ages to test
	for (Nr in maxAges:minAges){ #
		its <- length(agesUniv) - Nr + 1
		for (set in 1:its){ # 
			ind <- ind + 1
			agesL[[ind]] <- agesUniv[set:(set+Nr-1)]
		}
	}
	
	# these are the ages that give the best r2 or RMS
	# this would be much better using pipes, but hey, one less dependency...
	
	RMSE        <- lapply(agesL, ggbgetRMS, codi = codi, lm.method = lm.method) %>% unlist()
	minRMSE     <- which.min(RMSE)
	agesfit     <- agesL[[minRMSE]]
	
	# PG: also return RMSE?
	list(agesfit = agesfit, RMSE = RMSE[minRMSE])
}



#' @title estimate death registration coverage using the GGB method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method is based on finding a best-fitting linear relationship between two modeled parameters (right term and left term), but the fit, and resulting coverage estimate, depend on exactly which age range is taken. This function either finds a nice age range for you automatically, or you can specify an exact vector of ages. 
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$id} with unique values for each data chunk. Different values of \code{$id} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths for each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range in the traditional visual way, see \code{ggbChooseAges()}, when working with a single year/sex/region of data. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared residual between the fitted line and right term.
#' 
#' @param X \code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$age}, and \code{$id} (if there are more than 1 region/sex/intercensal period).
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param exact.ages optional. A user-specified vector of exact ages to use for coverage estimation
#' @param nx.method either 2 or 4. 4 is smoother.
#' @param deaths.summed logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @param mig.summed logical. Is the (optional) net migration column \code{mig} given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.
#' @inheritParams slopeint

#' @return a \code{data.frame} with columns for: \itemize{
#'   \item{id} group id
#'   \item{Mxcoverage} coverage of the intercensal Mx values: \code{sqrt(k1*k2)/b}
#'   \item{lower} lower bound of ages used for fitting
#'   \item{upper} upper bound of ages used for fitting
#'   \item{a} intercept
#'   \item{b} slope
#'   \item{delta} empirical link: \code{exp(t*a) = k1/k2}
#'   \item{k1} completeness of census 1
#'   \item{k2} completeness of census 2
#'   \item{k3} completeness of deaths relative to census 2 (\code{1/b})
#'   \item{t1} decimal date of census 1
#'   \item{t2} decimal date of census 2
#'   \item{t} intercensal interval (\code{t2 - t1})
#'   \item{lm.method} line fitting method used
#'   \item{nx.method} birthday (Nx) approximation used (2 or 4 points)
#'   \item{RMSE} the root mean square error of the optimized ages used for fitting (only if ages were automatically selected)
#' }
#' @export
#' @references 
#' Hill K. Estimating census and death registration completeness. Asian and Pacific Population Forum. 1987; 1:1-13.
#' 
#' Brass, William, 1975.  Methods for Estimating Fertility and Mortality from Limited and Defective Data, Carolina Population Center, Laboratory for Population Studies, University of North Carolina, Chapel Hill.
#' 
#' @examples 
#' # The Mozambique data
#' res <- ggb(Moz)
#' res
#' # The Brasil data
#' BM <- ggb(BrasilMales)
#' BF <- ggb(BrasilFemales)
#' head(BM)
#' head(BF)

ggb <- function(
		X, 
		minA = 15, 
		maxA = 75, 
		minAges = 8, 
		exact.ages = NULL, 
		lm.method = "oldschool",
		nx.method = 2,
		deaths.summed = FALSE,
		mig.summed = deaths.summed){         
	
	# TR: modularized Apr 2, 2017
	tab1        <- headerPrep(X)
	# iterate over whatever it happens to be: regions, years
	coverages   <- as.data.frame(
			      do.call(
					rbind,
					lapply(
					    tab1, 
					    ggbcoverageFromYear, 
					    exact.ages = exact.ages,
						  minA = minA, 
						  maxA= maxA,
						  minAges = minAges,
						  deaths.summed = deaths.summed,
						  mig.summed = mig.summed,
						  lm.method = lm.method,
						  nx.method = nx.method
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


#' @title does a given pairlist of x and y coordinates fall within the plot region?
#' 
#' @description Check to see if a point clicked falls in the plot or outside it. This function is used by \code{ggbChooseAges()}.
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


#' @title which age is closest to the point clicked?
#' @description a utility function called by \code{ggbChooseAges()}.
#' 
#' @param xvec \code{$rightterm}, as given by \code{ggbMakeColumns()}
#' @param yvec \code{$lefttterm}, as given by \code{ggbMakeColumns()}
#' @param click a point given by \code{locator(1)}
#' @param age ages present in dataset
#' 
#' @return the age corresponding to the x,y pair of \code{$rightterm}, \code{$lefttterm} closest to the point clicked.
#' 
#' @export
guessage <- function(xvec,yvec,click,age){
	age[which.min(sqrt((xvec - click$x) ^ 2 + (yvec - click$y) ^ 2))]
}


#' @title adjust the range of ages used by \code{ggbChooseAges()}
#' @description a utility function called by \code{ggbChooseAges()}. After clicking a point, this function readjusts the age range
#' 
#' @param a an age specified by the user, as returned by \code{guessage()}
#' @param age ages present in dataset
#' @param agesfit the age range used for calculating the coverage coefficient
#' 
#' @return the adjusted set of ages used for calculating the coverage coefficient
#' 
#' @export
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

#' @title get the slope the slope and intercept implied by a set of ages
#' @description Called by \code{ggbFittedFromAges()} and \code{ggbChooseAges()}
#' @param codi \code{data.frame} as produced by \code{ggbMakeColumns()}
#' @param agesfit a set of continuous ages to estimate coverage from 
#' @param lm.method character, one of:\itemize{
#'   \item{\code{"oldschool"}} default sd ratio operation of still unknown origin
#'   \item{\code{"lm"} or \code{"ols"}} for a simple linear model
#'   \item{\code{"tls"}, \code{"orthogonal"}, or \code{"deming"}} for total least squares
#'   \item{\code{"tukey"}, \code{"resistant"}, or "\code{"median"}} for Tukey's resistant line method
#' }
#' @return a pairlist with elements \code{$a} for the intercept and \code{$b} for the slope
#' @importFrom tukeyedar eda_rline
#' @importFrom stats lm
#' @importFrom stats cov
#' @importFrom stats var
#' @importFrom stats sd
#' @export
slopeint <- function(codi, 
                     agesfit, 
                     lm.method = "oldschool"){
  
  lm.method <- match.arg(lm.method, choices = c("oldschool", "lm", "ols", "tls","orthogonal","deming","tukey","resistant","median"))
  stopifnot("leftterm" %in% colnames(codi))

  codi <- codi %>% 
    filter(.data$age %in% agesfit)
  
  if (lm.method == "oldschool"){
  	#age <- codi$age
  	# TODO: find eq numbers to cite here
  	slope       <- 	sd(codi$leftterm) /  sd(codi$rightterm)

  	# PJ: fix https://github.com/timriffe/AdultCoverage/issues/3
  	intercept <- mean(codi$leftterm) - mean(codi$rightterm) * slope
  	
  	#coefs <- with(codi,lm(leftterm~rightterm))$coef
  }
  
  
  if (lm.method %in% c("lm","ols")){
    
   ab        <- lm("leftterm~rightterm", data = codi)$coef
   slope     <- ab[2]
   intercept <- ab[1]
  }
  
  if (lm.method %in% c("tls","orthogonal","deming")){
    y <- codi$leftterm
    x <- codi$rightterm
    delta <- 1 # assumed
    # let's avoid a new package dependency
    # PJ's implementation via Wikipedia formulas
    slope <- (var(y)-delta*var(x)+sqrt((var(y)-delta*var(x))^2+4*delta*cov(x,y)^2)) /
      (2*cov(x,y))
    
    intercept <- mean(y) - slope * mean(x)
  }
  
  if (lm.method %in% c("tukey","resistant","median")){
    ab.etc <- eda_rline(x = codi$rightterm, y = codi$leftterm)
    slope     <- ab.etc$b
    intercept <- ab.etc$a
  }
  
  # PG: return r2 also
	list(a = intercept, b = slope)
}



#' @title interactively determine ages to use for estimating coverage
#' @description In a spreadsheet one would typically set up the GGB method to produce a plot that updates as the user changes the age range. This function implements that kind of work flow. This will be intuitive for spreadsheet users, but it does not scale well. Imagine you have 200 territorial units, then you would not want to repeat this task. \code{ggb()} does the same thing automatically. You can compare the age range you select manually with the one given back by \code{ggb()} as a diagnostic, for instance. To set up the plot device, just give a single year/region/sex of data. By default it will give the RMSE-optimized age range to start with, but you can specify a  vector of exact ages to use as well. All points are plotted, with a fitted line that has been set to a subset of the points, which is plotted in a different color. You can click any point to change the age range, and the plot updates accordingly, up to a maximum of 15 clicks so you don't waste your time. You can stop the plot by either clicking on the graphics device outside the plot area or clicking out the 15 tries (or more if you increase \code{maxit}).
#' @details If you want to send the results of this into \code{ggb()}, you can do so by setting \code{Exact.ages} to \code{seq(lower,upper,by=5)}, where \code{$lower}, and \code{$upper} are the results returned from \code{ggbChooseAges()} after you're done manually determining the age range.
#' 
#' @inheritParams ggb
#' @param maxit up to how many times do you want to let yourself fiddle with the age range?
#' @return \code{data.frame} containing elements \code{$coverage}, \code{$lower}, \code{$upper}, and \code{ages}.
#' 
#' @importFrom grDevices gray
#' @importFrom graphics abline legend locator mtext par plot points rect segments text
#' @export
#' 
#' @examples
#' \dontrun{
#' # for interactive sessions only
#' # *click points to adjus age range used (yellow)
#' # *click in margin to stop and return coverage results
#' ggbChooseAges(Moz)
#' }

ggbChooseAges <- function(
              X, 
		          minA = 15, 
						  maxA = 75, 
						  minAges = 8, 
						  exact.ages = NULL, 
						  maxit = 15, 
						  deaths.summed = FALSE,
						  mig.summed = deaths.summed,
						  lm.method = "oldschool",
						  nx.method = 2){
	# this is the automatic age selection.
	
	# only run if in anteractive r session...
    stopifnot(interactive())
	
	# reset ages if necessary
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
  
  codi_in <- X
	codi    <- data.frame(X)           


	# guess which column is the deaths column, rename it deaths
	# codi    <- guessDeathsColumn(codi)
	
	# start GGB stuff
	codi    <- ggbMakeColumns(codi, 
	                          minA = minA, 
	                          maxA = maxA, 
	                          nx.method = nx.method,
	                          deaths.summed = deaths.summed,
	                          mig.summed = mig.summed)
	
	# some potential starting ages. either auto or self-supplied
	if (!is.null(exact.ages) & length(exact.ages) >= 3){
		agesfit <- exact.ages
	} else {
		fit.res <- ggbgetAgesFit(codi, 
		                         minA = minA, 
		                         maxA = maxA, 
		                         minAges = minAges,
		                         lm.method = lm.method)
		agesfit <- fit.res$agesfit
	}
	
	codi     <- ggbFittedFromAges(codi = codi, 
	                              agesfit = agesfit, 
	                              lm.method = lm.method)
	
	
	coverage <- ggbcoverageFromYear(codi = codi, 
	                                exact.ages = agesfit, 
	                                minA = minA, 
	                                maxA = maxA,
	                                minAges = minAges,
	                                deaths.summed = deaths.summed,
	                                mig.summed = mig.summed,
	                                lm.method = lm.method,
	                                nx.method = nx.method)

	
	age      <- codi$age
	# rm border ages re PJ issue #2
	ages.rm  <- age < minA | age > maxA 
	leftt    <- codi$leftterm
	rightt   <- codi$rightterm
	leftt[ages.rm]  <- NA
	rightt[ages.rm] <- NA
	
	
	# age ranges used for fitting
	amin    <- min(agesfit); amax <- max(agesfit)
	plot(rightt, leftt, asp = 1, pch = 19, col = "#00000050", cex = 1.6,
			xlab = "right term",
			ylab = "left term",
			main = paste0("Age range [",amin,
					",",amax,"], \nLine method = ",lm.method,"\nEst. coverage = ",round(coverage$Mxcoverage*100,1)),
			sub = "(optimized age range)")
	# automatically fit line (RMS of ggb)
	abline(a = coverage$a, b = coverage$b, col = "blue")
	# shows points used to fit line
	points(rightt[age %in% agesfit], 
			leftt[age %in% agesfit], col = "#FFFF00", pch = 19, cex = 1.6)
	text(rightt, leftt, age, cex = .6)
	legend("bottomright",lty=1,col="blue",legend="fitted line",bty="n")
	# message to user
	mtext("*Adjusting ages*\nClick any age to change bounds for fitting\nClick in margin to quit",
			side = 3, line = -3, outer = FALSE)
	
	
	USR    <- par()$usr
	agesfit <- agesfit # keep old one...
	for (i in 1:maxit){
		# now enter into interactive loop
		clicki <- locator(1)
		# was the click inside the plot or outside the plot?
		IN     <- inUSR(USR, clicki)
		if (IN){
			
			# if it's inside the plot, then which is the closest age?
			a        <- guessage(rightt,leftt,clicki,age)
			# guess how to adjust ages based on which age was clicked
			agesfit  <- adjustages(a, age, agesfit)
			# new range of ages used for fitting
			amin     <- min(agesfit)
			amax     <- max(agesfit)
			
			coverage <- ggbcoverageFromYear(codi = codi, 
			                                exact.ages = agesfit, 
			                                minA = minA, 
			                                maxA = maxA,
			                                minAges = minAges,
			                                deaths.summed = deaths.summed,
			                                mig.summed = mig.summed,
			                                lm.method = lm.method,
			                                nx.method = nx.method)
			
			# redraw plot
			plot(rightt, 
					leftt, 
					asp = 1, 
					pch = 19, 
					col = "#00000050", 
					cex = 1.6,
					xlab = "right term",
					ylab = "left term",
					main = paste0("Age range [", amin,
							",", amax, "], est. coverage = %",round(coverage$Mxcoverage * 100, 1)),
					sub = "(optimized age range)")
			# new fitted slope, intercept
		    abline(a=0,b=1,col=gray(.8)) # line of perfection
			#
			abline(a = coverage$a, b = coverage$b, col = "blue")
			# indicate which points used with color
			points(rightt[age %in% agesfit], 
					leftt[age %in% agesfit], col = "#FFFF00", pch = 19, cex = 1.6)
			text(rightt, leftt, age, cex = .6)
			legend("bottomright", lty = 1, col = "blue", legend = "fitted line", bty = "n")
		} else {
			break
		}
	}
	
	# click outside the margin to return results
    out <- ggb(X = codi_in,
               minA = minA, 
               maxA = maxA, 
               minAges = minAges, 
               exact.ages = agesfit, 
               lm.method = lm.method,
               nx.method =  nx.method,
               deaths.summed = deaths.summed,
               mig.summed = deaths.summed)
	out
}

# end




