#####################################################
b2hgetRMSbh2getRMS <- function(agesi, codi){
	# get root mean square of residuals
	sqrt(sum(codi$RelDiff[codi$age %in% agesi]^2)/length(agesi))
}

bh2CoverageFromAges <- function(codi, agesFit){
	inds    <- codi$age %in% agesFit
	sum(codi$Cx[inds]) / length(agesFit)
}

# change name to GB
bh2getRMS <- function(agesi, codi){
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
#head(codi)
bh2coverageFromYear <- function(codi, minA., AgeInt., minAges., ages.,sex.){
	dif.                   <- codi$year2[1] - codi$year1[1] 
	codi$pop1cum           <- rev(cumsum(rev(codi$pop1))) # Tx
	codi$pop2cum           <- rev(cumsum(rev(codi$pop2))) # Tx
	codi$deathcum          <- rev(cumsum(rev(codi$death))) # lx
	
	# define new column for birthdays between pop estimates
	codi$birthdays            <- 0
	# iterate over age groups >= 10
	
	for (j in seq_along(ages)[ages >= minA.]) {
		# take geometric average of p1 pop vs p2 pop within same cohort
		codi$birthdays[j]       <- round(
				1 / AgeInt. * sqrt(codi$pop1[j - 1] * codi$pop2[j]), 
				digits = 2)
	} # end age loop
	
	# create stationary Lx as geometric avg of within-cohort consecutive ages
	codi$Lx                <- sqrt(codi$pop1cum * codi$pop2cum)
	
	# growth rate per annum
	codi$cumgrowth        <- log(codi$pop2cum/ codi$pop1cum) / dif.
	
	# eqns from formula in Hill/Horiuchi
	codi$rightterm         <- codi$deathcum / codi$Lx
	# eqns from formula in Hill/Horiuchi
	codi$lefterm           <- (codi$birthdays / codi$Lx) - codi$cumgrowth
	# certain columns can be safely ignored in future operations
	codi$exclude          <-  codi$Lx != 0 & codi$birthdays != 0 & codi$age >= 15 & codi$age <= 75
	
	# intercept
	
	maxAges   <- sum(codi$exclude)
	agesUniv  <- ages.[codi$exclude]
	
	FirstAges <- agesUniv[agesUniv < 30]
	
	ind       <- 0
	agesL     <- list()
	# determine ages to test
	for (Nr in maxAges:minAges){ # Nr <- maxAges
		its <- length(agesUniv) - Nr + 1
		for (set in 1:its){ # set <- its[1]
			ind <- ind + 1
			agesL[[ind]] <- agesUniv[set:(set+Nr-1)]
		}
	}
	
#	# these are the ages that give the best r2 or RMS
#	if (fit. == "r2"){
#		agesfit     <- agesL[[which.max(unlist(lapply(agesL, getr2, codi = codi)))]]
#	}

	agesfit     <- agesL[[which.min(unlist(lapply(agesL, ggbgetRMS, codi = codi)))]]
	#agesfit     <- agesL[[which.max(unlist(lapply(agesL, getr2, codi = codi)))]]
	# this is the basic formula
	#coverageFromAges(codi, agesfit)

## copied and pasted
#slope       <- with(codi, 
#		sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
#)
#intercept   <-  with(codi, 
#		(mean(lefterm[age %in% agesfit]) * slope - mean(rightterm[age %in% agesfit]))
#) 
#codi$fitted <- codi$rightterm * slope + intercept
#
## this is the coverage estimate
#1/with(codi, sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit]))
## end copy and paste


    agesi       <- codi$age %in% agesfit

	slope       <- with(codi, 
						sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
					)
					#codi$lefterm
					#codi$rightterm
					#with(codi,sd(lefterm[age %in% agesfit]))
                    #with(codi,sd(rightterm[age %in% agesfit]))
	intercept   <-  with(codi, 
						(mean(lefterm[age %in% agesfit]) - slope * mean(rightterm[age %in% agesfit]))
					) 

	# relative completeness				
    relcomp     <- exp(intercept * dif.)
	
	# relcomp <- 1
	# adjust the first population count
# codi$pop1adj <- codi$rightterm * slope + intercept
	codi$pop1adj <- codi$pop1 / relcomp
	
	codi$birthdays            <- 0
# iterate over age groups >= 10
	
	for (j in seq_along(ages)[ages >= minA.]) {
		# take geometric average of p1 pop vs p2 pop within same cohort
		codi$birthdays[j]       <- 
				1 / AgeInt. * sqrt(codi$pop1adj[j - 1] * codi$pop2[j])
	} # end age loop
	
# age-specific growth
	
	codi[["growth"]]	          <-  log(codi$pop2 / codi$pop1adj) / dif.
	codi$growth[is.infinite(codi$growth)] <- 0
	
	codi$cumgrowth       <-  0
	codi$cumgrowth[1]    <-  2.5*codi$growth[1]
	
	for (j in 2:length(ages)){
		codi$cumgrowth[j]  <-  2.5*codi$growth[j]+5*sum(codi$growth[(j-1):1])
	}
	
# stopped here
	
	codi$death_tab       <- codi$death * exp(codi$cumgrowth)
	
	ratio           <- sum(codi$death_tab[ages%in%c(10:39)])/sum(codi$death_tab[ages%in%c(40:59)])
	
	if (sex. == "f"){
		# TODO: expand ex in-line out to actual open ages..
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
	if (sex. == "m"){
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
	AllLevels <- 3:25
	CDlevel   <- splinefun(AllLevels~standardratios)(ratio)
	# open at age 110 ....
	eOpen     <- splinefun(ex[,ncol(ex)]~1:25)(CDlevel)
	
	N         <- nrow(codi)
	codi$growth[is.nan(codi$growth)] <- 0
	if (sign( codi$growth[N ]) == -1){
		minus <- Re(((eOpen * codi$growth[N ])+.0i)^(1/3)) * 2
	} else {
		minus <- (eOpen * codi$growth[N ])^(1/3)
	}
	codi$pop_a     <- codi$death[N ] * (exp(eOpen * codi$growth[N ]) - minus)
	
	
	for(j in N:1){
		codi$pop_a[j - 1] <- codi$pop_a[j] * exp(AgeInt. * codi$growth[j-1]) + 
				codi$death[j - 1] * exp(AgeInt. / 2 * codi$growth[j - 1])
	}
	# plot(codi$pop_a, type = 'l')
	
	codi$Cx      <-  codi$pop_a / codi$birthdays
	codi$RelDiff <-  (codi$pop_a - codi$birthdays) / ( codi$pop_a + codi$birthdays) / 2
	
	codi$exclude          <-  codi$birthdays != 0 & codi$age >= 15 & codi$age <= 75 & !is.nan(codi$Cx) & !is.nan(codi$RelDiff)
	
	# intercept
	
	maxAges   <- sum(codi$exclude)
	agesUniv  <- codi$age[codi$exclude]
	
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
	# TODO: these ages might be different from 'agesfit'. Do we care?
	agesFit     <- agesL[[which.min(unlist(lapply(agesL, bh2getRMS, codi = codi)))]]
	
	bh2CoverageFromAges(codi, agesFit)
	
}

# TODO: detect AgeInt rather than specify as argument
bh2 <- function(x, minA = 10, AgeInt = 5, minAges = 8, sex = "f"){
	tab     <- data.frame(x)           ##  Dat in data frame : cod, age, pop1, year1, pop2, year2, death
	tab$pop1   <- as.double(tab$pop1)
	tab$pop2   <- as.double(tab$pop2)
	tab$death  <- as.double(tab$death)
	tab1    <- split(tab,x$cod)

	ages    <- sort(unique(tab$age))
	#minA. = minA;AgeInt. = AgeInt;minAges. = minAges;ages. = ages;sex. = sex
	# codi <- tab1[[1]]
	coverages <- unlist(lapply(
					tab1, 
					bh2coverageFromYear, 
					minA. = minA, 
					AgeInt. = AgeInt, 
					minAges. = minAges,  
					ages. = ages,
					sex. = sex)) # sex <- "f"
	coverages
}
#
#coveragesggb <- ggb(x)
#coveragesbh1 <- bh1(x)
#coveragesbh2 <- bh2(x)
#years <- as.integer(names(coveragesggb))
#getwd()
#
#pdf("Figures/ITAfemalesCompareMethods.pdf", height = 5, width = 10)
#plot(years, coveragesggb, type = 'l', col = "green", ylim = c(.7,1.4), main = "Italy, females 1872-2008", sub = "Death registration coverage detected using 3 methods")
#lines(years, coveragesbh1, col = "blue")
#lines(years, coveragesbh2, col = "red")
#abline(v=c(1872, 1924, 1947, 1954, 1981),col="magenta")
#abline(h=1,col="magenta")
#legend("topright", col = c("green","blue","red"), lty=1,legend=c("GGB","BH unadj", "BH adj"))
#dev.off()