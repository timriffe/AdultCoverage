
# 
# Author: riffe
###############################################################################


# deprecated, doesn't give good results
#ggbgetr2 <- function(agesi, codi){
#	slope       <- with(codi, 
#			sd(lefterm[age %in% agesi]) /  sd(rightterm[age %in% agesi])
#	)
#	intercept   <-  with(codi, 
#			(mean(lefterm[age %in% agesi]) * slope - mean(rightterm[age %in% agesi]))
#	) 
#	codi$fitted <- codi$rightterm * slope + intercept
#	
#	summary(lm(codi$lefterm[codi$age %in% agesi] ~ codi$fitted[codi$age %in% agesi]))$r.squared
#}

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

# codi <- tab1[[1]]
ggbcoverageFromYear <- function(codi, minA., AgeInt., minAges., ages.){
	
	codi    <- ggbMakeColumns(codi, minA., AgeInt., minAges., ages.)
	agesfit <- ggbgetAgesFit(codi, ages., minAges.)
	# this is the basic formula
    ggbcoverageFromAges(codi, agesfit)
	
}

ggbColumnsFromAges <- function(codi, agesfit){
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

ggbMakeColumns <- function(codi, minA., AgeInt., minAges., ages.){
	dif.                   <- codi$year2[1] - codi$year1[1] 
	codi$pop1cum           <- rev(cumsum(rev(codi$pop1))) # Tx
	codi$pop2cum           <- rev(cumsum(rev(codi$pop2))) # Tx
	codi$deathcum          <- rev(cumsum(rev(codi$death))) # lx
	
	# define new column for birthdays between pop estimates
	codi$birthdays         <- 0
	# iterate over age groups >= 10
	
	for (j in seq_along(ages.)[ages. >= minA.]) {
		# take geometric average of p1 pop vs p2 pop within same cohort
		codi$birthdays[j] <- 1 / AgeInt. * sqrt(codi$pop1[j - 1] * codi$pop2[j])
				
	} # end age loop
	
	# create stationary Lx as geometric avg of within-cohort consecutive ages
	codi$Lx               <- sqrt(codi$pop1cum * codi$pop2cum)
	
	# growth rate per annum
	codi$cumgrowth        <- log(codi$pop2cum / codi$pop1cum) / dif.
	
	# eqns from formula in Hill/Horiuchi
	codi$rightterm        <- codi$deathcum / codi$Lx
	# eqns from formula in Hill/Horiuchi
	codi$lefterm          <- (codi$birthdays / codi$Lx) - codi$cumgrowth
	# certain columns can be safely ignored in future operations
	codi$exclude          <-  codi$Lx != 0 & codi$birthdays != 0 & codi$age >= 15 & codi$age <= 75
	codi
}

ggbgetAgesFit <- function(codi, ages., minAges.){
	
#	codi <- ggbMakeColumns(codi, minA., AgeInt., minAges., ages.)
	# intercept
	
	maxAges   <- sum(codi$exclude)
	agesUniv  <- ages.[codi$exclude]
	
	FirstAges <- agesUniv[agesUniv < 30]
	
	ind       <- 0
	agesL     <- list()
	# determine ages to test
	for (Nr in maxAges:minAges.){ # Nr <- maxAges
		its <- length(agesUniv) - Nr + 1
		for (set in 1:its){ # set <- its[1]
			ind <- ind + 1
			agesL[[ind]] <- agesUniv[set:(set+Nr-1)]
		}
	}
	
	# these are the ages that give the best r2 or RMS

	agesfit     <- agesL[[which.min(unlist(lapply(agesL, ggbgetRMS, codi = codi)))]]
	agesfit
}




# codi <- tab1[[50]]
# names(tab1)
#minA. = 10; AgeInt. = 5; minAges. = 8; ages. = ages
ggb <- function(x, minA = 10, AgeInt = 5, minAges = 8){         ##  Data
	tab      <- data.frame(x)           ##  Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)
	# TR: account for decimal intervals
	#cod      <- factor(x$cod)
	tab$pop1   <- as.double(tab$pop1)
	tab$pop2   <- as.double(tab$pop2)
	tab$death  <- as.double(tab$death)
	tab1       <- split(tab, factor(x$cod))
	# length(tab1)
	# iterate over whatever it happens to be: regions, years
	ages <- sort(unique(tab$age))
	
	coverages <- unlist(lapply(
					tab1, 
					ggbcoverageFromYear, 
					minA. = minA, 
					AgeInt. = AgeInt, 
					minAges. = minAges, 
					ages. = ages))
	
	return(coverages)
}
#  head(x)
#coveragesr2 <- ggb(x, fit = "r2")
#coverageggb <- ggb(x, fit = "RMS")
#plot(as.integer(names(coverages)), coverages, type = 'l', log = 'y' )
#abline(h = 1)
#rect(1872,)
#plot(coverages,log='y')
#abline(h=1)





inUSR <- function(USR, click){
	x <- click$x
	y <- click$y
	if (x < USR[1] | x > USR[2] | y < USR[3] | y > USR[4]){
		return(FALSE)
	}
	TRUE
}


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

plot.ggb <- function(codi, minA. = 10, AgeInt. = 5, minAges. = 8,maxit=15){
	
	# generate the basic components used to get GGB estimates,
	# using the automatic process to find optimal ages
	# used for adjustment
	ages.    <- sort(unique(codi$age))
	codi     <- ggbMakeColumns(codi, minA., AgeInt., minAges., ages.)
	# this is the automatic age selection.
	agesfit  <- ggbgetAgesFit(codi, ages., minAges.)
	si       <- slopeint(codi, agesfit)
	codi     <- ggbColumnsFromAges(codi, agesfit)
	
	# this is the basic formula
	coverage <- ggbcoverageFromAges(codi, agesfit)
	
	# some objects used throughout
	age     <- codi$age
	leftt   <- codi$lefterm
	rightt  <- codi$rightterm
	fitd    <- codi$fitted
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