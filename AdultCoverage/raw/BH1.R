bhgetRMS <- function(agesi, codi){
	# get root mean square of residuals
	sqrt(sum(codi$RelDiff[codi$age %in% agesi]^2)/length(agesi))
}

bhCoverageFromAges <- function(codi, agesFit){
	inds    <- codi$age %in% agesFit
	sum(codi$Cx[inds]) / length(agesFit)
}

bhCoverageFromYear <-  function(codi, minA. = 10, AgeInt. = 5, minAges. = 8, ages., sex. = "f"){        ##  Data
	
	dif. <- codi$year2[1] - codi$year1[1]
	
	codi$birthdays            <- 0
# iterate over age groups >= 10
	
	for (j in seq_along(ages)[ages >= minA.]) {
		# take geometric average of p1 pop vs p2 pop within same cohort
		codi$birthdays[j]       <- 
				1 / AgeInt. * sqrt(codi$pop1[j - 1] * codi$pop2[j])
	} # end age loop
	
# age-specific growth
	
	codi[["growth"]]	          <-  log(codi$pop2 / codi$pop1) / dif.
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
		# TODO: expand ex in-ine out to actual open ages..
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
	
	N              <- nrow(codi)
	codi$growth[is.nan(codi$growth)] <- 0
	if (sign( codi$growth[N ]) == -1){
		minus <- Re(((eOpen * codi$growth[N ])+.0i)^(1/3)) * 2
	} else {
		minus <- (eOpen * codi$growth[N ])^(1/3)
	}
	codi$pop_a     <- codi$death[N ] * (exp(eOpen * codi$growth[N ]) -minus)
	
	
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
	
	agesFit     <- agesL[[which.min(unlist(lapply(agesL, bhgetRMS, codi = codi)))]]
	
	bhCoverageFromAges(codi, agesFit)
	
}

# TODO: detect AgeInt rather than specify as argument
bh <- function(x, minA = 10, AgeInt = 5, minAges = 8, sex = "f"){
	tab     <- data.frame(x)           ##  Dat in data frame : cod, age, pop1, year1, pop2, year2, death
	tab1    <- split(tab,x$cod)
	cods    <- unique(tab$cod)
	ages    <- sort(unique(tab$age))
	#minA. = minA;AgeInt. = AgeInt;minAges. = minAges;ages. = ages;sex. = sex
	# codi <- tab1[[1]]
	coverages <- unlist(parallel::mclapply(
					tab1, 
					bhCoverageFromYear, 
					minA. = minA, 
					AgeInt. = AgeInt, 
					minAges. = minAges,  
					ages. = ages,
					sex. = sex,
					mc.cores = 4))
	coverages
}
