#####################################################
getr2 <- function(agesi, codi){
  slope       <- with(codi, 
    sd(lefterm[age %in% agesi]) /  sd(rightterm[age %in% agesi])
  )
  intercept   <-  with(codi, 
    (mean(lefterm[age %in% agesi]) * slope - mean(rightterm[age %in% agesi]))
  ) 
  codi$fitted <- codi$rightterm * slope + intercept
  
  summary(lm(codi$lefterm[codi$age %in% agesi] ~ codi$fitted[codi$age %in% agesi]))$r.squared
}
getRMS <- function(agesi, codi){
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
coverageFromYear <- function(codi, minA., AgeInt., minAges., dif., ages., fit. = "RMS"){
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
  codi$Lx                <- round(
    sqrt(codi$pop1cum * codi$pop2cum),
    digits = 2)
  
  # growth rate per annum
  codi$cumgrowth        <- round(
    log(codi$pop2cum/ codi$pop1cum) / dif.,
    digits = 5 )
  
  # eqns from formula in Hill/Horiuchi
  codi$rightterm         <- round(
    codi$deathcum / codi$Lx, 
    digits = 5)
  # eqns from formula in Hill/Horiuchi
  codi$lefterm           <- round(
    (codi$birthdays / codi$Lx) - codi$cumgrowth,
    digits = 5)
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
  
  # these are the ages that give the best r2 or RMS
  if (fit. == "r2"){
    agesfit     <- agesL[[which.max(unlist(lapply(agesL, getr2, codi = codi)))]]
  }
  if (fit. == "RMS"){
    agesfit     <- agesL[[which.min(unlist(lapply(agesL, getRMS, codi = codi)))]]
  }
  #agesfit     <- agesL[[which.max(unlist(lapply(agesL, getr2, codi = codi)))]]
  # this is the basic formula
  coverageFromAges(codi, agesfit)
  
}

coverageFromAges <- function(codi, agesfit){
  slope       <- with(codi, 
    sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit])
  )
  intercept   <-  with(codi, 
    (mean(lefterm[age %in% agesfit]) * slope - mean(rightterm[age %in% agesfit]))
  ) 
  codi$fitted <- codi$rightterm * slope + intercept
  
  # this is the coverage estimate
  1/with(codi, sd(lefterm[age %in% agesfit]) /  sd(rightterm[age %in% agesfit]))
}

# codi <- tab1[["2003"]]
# names(tab1)
ggb <- function(x, minA = 10, AgeInt = 5, minAges = 8, fit = "RMS"){         ##  Data
  tab      <- data.frame(x)           ##  Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)
  dif      <- x$year2[1] - x$year1[1] # TR: account for decimal intervals
  #cod      <- factor(x$cod)
  x$pop1   <- as.double(x$pop1)
  x$pop2   <- as.double(x$pop2)
  x$death  <- as.double(x$death)
  tab1     <- split(tab, factor(x$cod))
  
  # iterate over whatever it happens to be: regions, years
  cods <- unique(tab$cod)
  ages <- sort(unique(tab$age))
  
  coverages <- unlist(parallel::mclapply(
      tab1, 
      coverageFromYear, 
      minA. = minA, 
      AgeInt. = AgeInt, 
      minAges. = minAges, 
      dif. = dif, 
      ages. = ages,
      fit. = fit,
      mc.cores = 4))

    return(coverages)
  }

#  head(x)
#coveragesr2 <- ggb(x, fit = "r2")
#coveragesRMS <- ggb(x, fit = "RMS")
#plot(as.integer(names(coverages)), coverages, type = 'l', log = 'y' )
#abline(h = 1)
#rect(1872,)
#
#
#plot(coveragesr2, coveragesRMS, asp = 1)
#
#plot(as.integer(names(coveragesr2)), coveragesr2, type = 'l')
#lines(as.integer(names(coveragesRMS)), coveragesRMS, col = "red")
#abline(h=1)
#
#var(coveragesRMS)


