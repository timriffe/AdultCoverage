
# Author: tim
###############################################################################
# contains functions related to the generalized growth balance method


#' @title calculate the root means square of the error to help find optimal age range
#' 
#' @description Called by `ggbgetAgesFit()` whenever the user does not want to manually determine the age range used to determine registration coverage. Probably no need to be called by top-level users. If a user would rather determine the optimal age range some other way, then look to` ggbcoverageFromYear() `where `ggbgetResidual` is called and add another condition or make it call something else.
#' 
#' @details Given a vector of ages and a line fitting method, we take one of several potential residuals. RMSE is the root of the mean squared error, MAE is the mean absolute error, MAPE is the mean absolute proportional error, ORSS is the standard deviation of the orthogonal residuals, and r2 is the r2 of an OLS fit.
#' 
#' @param agesi the vector of ages used for this iteration
#' @param codi `data.frame` as returned by `ggbMakeColumns()`
#' @param opt.method What should we try to minimize when picking age trims? Current options are `"RMSE"`, `"ORSS"`, `"MAE"`, `"MAPE"`, or `"r2"`. Default `"r2"`.
#' @param scale scale factor for the objective
#' @inheritParams slopeint
#' 
#' @return the RMSE
#' 
#' @importFrom stats prcomp
#' @export 
ggbgetResidual <- function(agesi, 
                      codi, 
                      lm.method = "oldschool", 
                      opt.method = "RMSE",
                      scale = 1){

  opt.method <-
    match.arg(opt.method,
              choices = c("RMSE","ORSS","r2","MAE","MAPE"))
	codi <- ggbFittedFromAges(codi, 
	                          agesfit = agesi, 
	                          lm.method = lm.method) %>% 
	  filter(.data$age %in% agesi)

	# if (opt.method == "logRMS"){
	#   #out <- sqrt(mean((log(codi$leftterm) - log(codi$fitted))^2) )
	#   out <- var(expit(codi$leftterm) - expit(codi$fitted)) / length(agesi)
	# }
	if (opt.method == "ORSS"){
	  out <- prcomp(cbind(codi$fitted,codi$leftterm))$sdev[1] * scale
	}
	# if (opt.method == "logORSS"){
	#   out <- prcomp(cbind(log(codi$fitted),log(codi$leftterm)))$sdev[1]
	# }
	if (opt.method == "r2"){
	  out <- lm(codi$fitted~codi$leftterm) %>% 
	    summary() %>% 
	    '[['("r.squared") %>% 
	    
	    '*'(-scale)
	}
	
	if (opt.method == "RMSE"){
	  out <- sqrt(mean((codi$leftterm - codi$fitted)^2) ) * scale
	}
	if (opt.method %in% c("RMSE","MAE","MAPE")){
	  resids <- codi$leftterm - codi$fitted
	  if (opt.method == "RMSE"){
	    out <- sqrt(mean(resids^2) ) * scale
	  }
	  if (opt.method == "MAE"){
	    out <- mean(abs(resids))
	  }
	  if (opt.method == "MAPE"){
	    out <- mean(abs(resids / codi$leftterm))
	  }
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
								opt.method = "RMSE",
								scale = 1,
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
								 lm.method = lm.method,
								 opt.method = opt.method,
								 scale = scale)
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

	rsq <- ggbgetResidual(agesi = agesfit,
	                 codi = codi,
	                 lm.method = lm.method,
	                 opt.method = "r2")
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
			   nx.method = nx.method,
			   opt.method = opt.method,
			   r2 = -rsq)

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
    mutate(	date1          = ifelse(is.numeric(.data$date1), 
                                    .data$date1, 
                                    decimal_date(.data$date1)),
            date2          = ifelse(is.numeric(.data$date2), 
                                    .data$date2, 
                                    decimal_date(.data$date2)),
            AgeInt         = age2int(.data$age),
            dif            = .data$date2 - .data$date1,
            avg            = deaths.summed,
            deathsAvg      = ifelse(.data$avg, 
                                    .data$deaths / .data$dif, 
                                    .data$deaths),
            pop1           = as.double(.data$pop1),
            pop2           = as.double(.data$pop2),
            mig            = as.double(.data$mig),
            avg            = mig.summed,
            migAvg         = ifelse(.data$avg, 
                                    .data$mig / .data$dif, 
                                    .data$mig),
            pop1cum        = lt_id_L_T(.data$pop1),
            pop2cum        = lt_id_L_T(.data$pop2),
            deathcum       = lt_id_L_T(.data$deathsAvg) * .data$dif,
            migcum         = lt_id_L_T(.data$migAvg) * .data$dif,
            Nx             = est_birthdays(pop1 = .data$pop1, 
                                          pop2 = .data$pop2, 
                                          AgeInt = .data$AgeInt, 
                                          nx.method = nx.method) * .data$dif,
            PYL = .data$dif * exp((log(.data$pop1cum) + log(.data$pop2cum)) / 2),
            bxp = .data$Nx / .data$PYL,
            rxp_m_ixp = (.data$pop2cum - .data$pop1cum - .data$migcum) / .data$PYL,
            rightterm = .data$deathcum / .data$PYL,
            leftterm = .data$bxp - .data$rxp_m_ixp,
            keep = .data$age >= minA & 
              .data$age <= maxA & 
              .data$Nx !=0 & 
              !is.na(.data$Nx))
  codi
}



#' @title determine the age range that minimizes the mean squared error
#' @description Called by `ggbcoverageFromYear()` whenever `exact.ages` are not given. This automates what one typically does visually.
#' @seealso code{\link{ggbChooseAges}}
#' 
#' @param codi `data.frame` as returbed by `ggbMakeColumns()`
#' @inheritParams ggb
#' 
#' @return a vector of ages that minimizes the RMSE
#' 
#' @export

ggbgetAgesFit <- function(codi, 
                          minA = 5, 
                          maxA = 75, 
                          minAges = 8,
                          lm.method = "tukey",
                          opt.method = "r2",
                          scale = 1){
	
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
	if (opt.method == "hybrid"){
	  hyb <- TRUE
	  opt.method <- "RMS"
	} else {
	  hyb <- FALSE
	}
	RMSE        <- lapply(agesL, 
	                      ggbgetResidual, 
	                      codi = codi, 
	                      lm.method = lm.method,
	                      opt.method = opt.method,
	                      scale = scale) %>% unlist()
	minRMSE     <- which.min(RMSE)
	agesfit     <- agesL[[minRMSE]]
	
	if (hyb){
	  opt.method <- "ORSS"
	  RMSE        <- lapply(agesL, 
	                        ggbgetResidual, 
	                        codi = codi, 
	                        lm.method = lm.method,
	                        opt.method = opt.method,
	                        scale = scale) %>% unlist()
	  minRMSE     <- which.min(RMSE)
	  agesfit2     <- agesL[[minRMSE]]
	  agesfit <- c(agesfit,agesfit2) %>% unique() %>% sort()
	}

	# PG: also return RMSE?
	list(agesfit = agesfit)
}



#' @title estimate death registration coverage using the GGB method
#' 
#' @description Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method is based on finding a best-fitting linear relationship between two modeled parameters (right term and left term), but the fit, and resulting coverage estimate, depend on exactly which age range is taken. This function either finds a nice age range for you automatically, or you can specify an exact vector of ages. 
#' 
#' @details Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$id} with unique values for each data chunk. Different values of \code{$id} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths for each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range in the traditional visual way, see \code{ggbChooseAges()}, when working with a single year/sex/region of data. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared residual between the fitted line and right term.
#' 
#' @param X `data.frame` with columns, `pop1`, `pop2`, `deaths`, `mig` (optional), `date1`, `date2`, `age`, and `id` (if there are more than 1 region/sex/intercensal period).
#' @param minA the lowest age to be included in search
#' @param maxA the highest age to be included in search (the lower bound thereof)
#' @param minAges the minimum number of adjacent ages to be used in estimating
#' @param exact.ages optional. A user-specified vector of exact ages to use for coverage estimation
#' @param nx.method either 2 or 4. 4 is smoother.
#' @param deaths.summed logical. is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume `FALSE`, i.e. that the average annual was given.
#' @param mig.summed logical. Is the (optional) net migration column `mig` given as the total per age in the intercensal period (`TRUE`). By default we assume `FALSE`, i.e. that the average annual was given.
#' @param opt.method what kind of residual do we minimize? choices `"RMS"`,`"logRMS"`, `"ORSS"`, `"logORSS"` (experimental)
#' @param scale multiplicative scale factor for the minimized residual
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
#'   \item{r2} the r2 of the optimized ages used for fitting (only if ages were automatically selected)
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
		lm.method = "tukey",
		opt.method = "r2",
		scale = 1,
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
						  nx.method = nx.method,
						  opt.method = opt.method,
						  scale = scale
					)))
	
	# this has cod as a column, but no year, sex. 
	return(coverages)
}






