#bh1getRMS <- function(agesi, codi){
#	# get root mean square of residuals
#	sqrt(sum(codi$RelDiff[codi$age %in% agesi]^2)/length(agesi))
#}

bh1CoverageFromAges <- function(codi, agesFit){
	inds    <- codi$age %in% agesFit
	sum(codi$Cx[inds]) / length(agesFit)
}
# ages. <- unique(codi$age); sex. <- "m"
bh1MakeColumns <- function(codi, minA. = 10, AgeInt. = 5, minAges. = 8, ages., sex. = "f"){
	dif. <- codi$year2[1] - codi$year1[1]
	
	codi$birthdays            <- 0
# iterate over age groups >= 10
	# dif. <- 10.00137
	for (j in seq_along(ages.)[ages. >= minA.]) {
		# take geometric average of p1 pop vs p2 pop within same cohort
		codi$birthdays[j]       <- 
				1 / AgeInt. * sqrt(codi$pop1[j - 1] * codi$pop2[j])
	} # end age loop
	
# age-specific growth

	codi[["growth"]]	          <-  log(codi$pop2 / codi$pop1) / dif.
	codi$growth[is.infinite(codi$growth)] <- 0
	
	codi$cumgrowth         <-  0
	codi$cumgrowth[1]      <-  AgeInt. / 2 * codi$growth[1]
	
	for (j in 2:length(ages.)){
		codi$cumgrowth[j]  <-  AgeInt. / 2 * codi$growth[j] + AgeInt. * sum(codi$growth[(j - 1):1])
	}
	
# stopped here
	
	codi$deathLT       <- codi$death * exp(codi$cumgrowth)
	
	ratio                <- sum(codi$deathLT[ages. %in% c(10:39)]) / sum(codi$deathLT[ages. %in% c(40:59)])
	
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
	# which is the open age?
	OA        <- max(ages.)
	availages <- as.integer(colnames(ex))
	
	stopifnot(OA %in% availages)
	
	eOpen     <- splinefun(ex[,as.character(OA)]~1:25)(CDlevel)
	# eOpen <- 4.9
	N         <- nrow(codi)
	codi$growth[is.nan(codi$growth)] <- 0
#	if (sign( codi$growth[N]) == -1){
#		minus <- Re(((eOpen * codi$growth[N]) + .0i)^(2 / 6)) * 2
#	} else {
#		minus <- (eOpen * codi$growth[N])^(2 / 6)
#	}
	minus          <- ((eOpen * codi$growth[N])^2) / 6
	codi$pop_a     <- codi$death[N] * (exp(eOpen * codi$growth[N]) - minus)
	
	
	for(j in N:2){
		codi$pop_a[j - 1] <- codi$pop_a[j] * exp(AgeInt. * codi$growth[j - 1]) + 
				codi$death[j - 1] * exp(AgeInt. / 2 * codi$growth[j - 1])
	}
	
	
	codi$Cx               <-  codi$pop_a / codi$birthdays
		
	codi
}

bh1CoverageFromYear <-  function(codi, minA. = 10, AgeInt. = 5, minAges. = 8, ages., sex. = "f", exact.ages. = NULL){        ##  Data

    # this is a test
	if (is.null(exact.ages.)){
		agesFit <- ggbgetAgesFit(ggbMakeColumns(codi, minA., AgeInt., minAges., ages.), 
				ages., minAges.)
	} else {
		agesFit <- exact.ages.
	}
	codi    <- bh1MakeColumns(codi, minA. = 10, AgeInt. = 5, minAges. = 8, ages., sex. = "f")
	bh1CoverageFromAges(codi, agesFit)
	# agesFit <- seq(15,65,by=5)
}

# names(tab1)
#x <- BR3[BR3$sex == "f", ]
#minA. = 10; AgeInt. = 5; minAges. = 8; ages. = unique(codi$age);.exact.ages=exact.ages
# exact.ages <- seq(30,65,by=5)
#codi <- tab1[[55]]
#names(tab1)[80]
# TODO: detect AgeInt rather than specify as argument
bh1 <- function(x, minA = 10, AgeInt = 5, minAges = 8, sex = "f", exact.ages = NULL){
	
	tab        <- data.frame(x)           ##  Dat in data frame : cod, age, pop1, year1, pop2, year2, death
	tab$pop1   <- as.double(tab$pop1)
	tab$pop2   <- as.double(tab$pop2)
	tab$death  <- as.double(tab$death)
	tab1       <- split(tab,x$cod)

	ages       <- sort(unique(tab$age))
	#minA. = minA;AgeInt. = AgeInt;minAges. = minAges;ages. = ages;sex. = sex
	# codi <- tab1[[1]]
	coverages <- unlist(lapply(
					tab1, 
					bh1CoverageFromYear, 
					minA. = minA, 
					AgeInt. = AgeInt, 
					minAges. = minAges,  
					ages. = ages,
					sex. = sex,
					exact.ages. = exact.ages))
	coverages
}
