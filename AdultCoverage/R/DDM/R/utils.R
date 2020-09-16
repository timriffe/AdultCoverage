
# Author: tim
###############################################################################
# contains utilities used throught.


#' @title Detect the sex for some demographic data
#' 
#' @description The column name can be \code{"sex"} or \code{"Sex"} and nothing else. If coded with integer, the number 1 is recognized as male and numbers, 0, 2, or 6 are assumed to be female. Any other integer will throw an error. If character, if the first letter is \code{"f"}, then we assume female, and if the first letter is \code{"m"} we assume male. Case does not matter. Anything else will throw an error. This function allows for just a little flexibility.
#' 
#' @param Dat a \code{data.frame} containing a column called \code{Sex}, or code{sex}.
#' @param sexColumn character string giving the name of the Sex column \code{"Sex"} assumed.

#' @return either \code{"f"} or \code{"m"}
#' @export 


detectSex <- function(Dat, sexColumn = "Sex"){
	
	tempnames <- tolower(colnames(Dat))
	stopifnot(tolower(sexColumn) %in% tempnames)
	
	# make column easier to call
	colnames(Dat)[grepl(tolower(sexColumn),tempnames)] <- "Sex"

	sex <- unique(Dat$Sex)
	
	# if it's integer, accept 0,2,6 for female, 1 for male
	if (is.integer(sex)){
		stopifnot(sex %in% c(0,1,2,6))
		sex <- ifelse(sex %in% c(0,2,6),"f","m")
	}
	if (is.character(sex)){
		sex <- tolower(sex)
		sex <- unlist(strsplit(sex,""))[1]
		if (sex %in% c("m","f")){
			return(sex)
		} else {
			stop("must indicate sex using 'm' and 'f'")
		}
	}
	stop("must specify sex using a character class column called sex, with values 'm' and/or 'f'")
}


#' @title append a \code{$cod} column if missing
#' @description Only handles the case of missing \code{$cod} splitting variable for data of a single year/region. This is not super robust. If you have many regions or whatever then do it yourself. This function was just written to make \code{ggb()} robust to the case of a user specifying data that don't have any territorial or other subgroups, aside from sex.
#' 
#' @param X a \code{data.frame}, possibly but not necessarily with column \code{$sex}.
#' @return X with a new column, \code{$cod} appended. 
#' 
#' @export

addcod <- function(X){
	stopifnot(is.data.frame(X))
	
	if (!"cod" %in% colnames(X)){
		X$cod <- 1
		if ("sex" %in% colnames(X)){
			sexes         <- unique(X$sex)
			recvec        <- 1:length(sexes)
			names(recvec) <- sexes
			X$cod         <- recvec[X$sex]
		}
	}
    X
}


#' @title Figure out which column is the deaths column
#' @description This function will pick up \code{"death"}, \code{"deaths"}, \code{"Death"}, or \code{"Deaths"} (and maybe some others?) and rename it \code{"deaths"} for easier internal usage.
#' 
#' @param X \code{data.frame} that we want to check for a deaths column
#' 
#' @return The same \code{data.frame}, returned, with the deaths column renamed as \code{"deaths"}
#' 
#' @export
guessDeathsColumn <- function(X){
	tempcols <- tolower(colnames(X))
	dcol     <- grepl("death", tempcols)
	stopifnot(any(dcol))
	colnames(X)[dcol] <- "deaths"
	X
}


#' @title chop down or group down the open age
#' 
#' @description These methods are not intended to be applied to ages greater than, say 90 or 95. Usually, we'd top out in the range 75 to 85. In any case, the Coale-Demeny life table implementation that we have only goes up to age 95, so there is a practical limitation to deriving a remaining life expectancy for the open age group. If a user tries to apply the Bennett-Horiuchi methods to data with higher open ages, stuff breaks for the time being. So this function chops the data off at \code{min(maxA,95)}, after having (optionally) grouped data down. This function needs to work with a single partition of data (intercensal period, sex, region, etc).
#' 
#' @param X data formatted per the requirements of \code{bh1()}, \code{bh2()}
#' @param group logical. If \code{TRUE} we sum down to \code{min(maxA,95)}. If \code{FALSE}, we just chop off data above that age.
#' @param maxA integer ignore ages above this age.
#' 
#' @return X, with the open age having been reduced either with or without aggregation.
#' @export 

reduceOpen <- function(X, maxA = 75, group = TRUE){
	ages   <- X$age
	topper <- min(maxA,95)
	if (max(ages) > min(maxA,95)){
		if (group){
			colsgroup <- c("pop1","pop2","deaths")
			if ("deathsAvg" %in% colnames(X)){
				colsgroup <- c(colsgroup, "deathsAvg")
			}
			X[X$age == topper, colsgroup] <- 
					colSums(X[X$age >= topper, colsgroup])
			
		}
		X <- X[X$age <= topper, ]
	}
	X
}

#' @title group down standard abridged data in child mortality group
#' @description We want 5-year age groups starting from 0. Standard abridged data has 0i,1,5. So we need to group together 0 and 1. Just for the sake of getting comparable results.
#' 
#' @param X standard input as required by \code{ddm()}, \code{ggb()}, \code{bh1()}, or \code{bh2()}
#' @return X, with child ages grouped as necessary (or not)
#' 
#' @export
group01 <- function(X){
	ages  <- X$age
	if (1 %in% ages){
		ind0 <- ages == 0
		ind1 <- ages == 1
		cnames <- c("pop1","pop2","deaths")
		if ("deathsAvg" %in% colnames(X)){
			cnames <- c(cnames,"deathsAvg")
		}
		# TR: look for age interval column and do that too?
		X[ind0,cnames] <- colSums(X[ ind0 | ind1, cnames])
		
		X <- X[!ind1, ]
	}
    X
}



#' @title a utility function to prep the header
#' @description This is an internal utility function, to save on redundant lines of code. Not so useful for hand-processing.
#' @param X this is any codi-style \code{data.frame}
#' @return a list of codi chunks (by intercensal period, region, etc), with standardized names, dates, etc.
#' @export

headerPrep <- function(X){
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
	
	tab1
}



#' helper function to calculate Nx (birthdays)
#' @description Approximate average birthdays per annum in a reference period and age group between two censuses. We do this by taking geometric means within ages, either from age x in census 1 and age x + 5 in census 2 (\code{nx.method = 2}), or from ages x and x+5 in both censuses. We also divide out the age group width, presumably 5.
#' 
#' @param pop1 numeric vector of population counts by age from census 1. 
#' @param pop2 numeric vector of population counts by age from census 2. 
#' @param AgeInt numeric vector of age interval widths
#' @param nx.method either 2 or 4. 4 is smoother.
#' @return numeric vector of estiamte birthdays per annum per age
#' @export

est_birthdays <- function(pop1, pop2, AgeInt = NULL, nx.method = 2){
  stopifnot(nx.method %in% c(2,4))
  
  pop1l <- log(pop1)
  pop2l <- log(pop2)
  N <- length(pop1)
  if (nx.method == 2){
    birthdays <- c(0, exp((pop1l[ -N  ] + pop2l[ -1 ])/2)) / AgeInt
  }
  if (nx.method == 4){
    birthdays <-   c(0,exp((pop1l[ -N  ] + pop1l[ -1 ] +  pop2l[ -N  ] + pop2l[ -1 ]) / 4)) / AgeInt
  }
  birthdays
}

