
# Author: tim
###############################################################################
# contains utilities used throught.

#' @title Detect the age interval for some demographic data
#' 
#' @description Since death distribution methods are primarily used in adult ages, it's OK to chop off the irregular infant and child age intervals (0,1], (1,5]. Further, if high ages are in different intervals this might also be a non-issue. In principal, the user should set \code{MinAge} and \code{MaxAge} to the same values used in the death distribution methods. Here we have some defaults that should almost always return the result \code{5} for standard abridged data, or \code{1} for single age data. Really there are not any other common age-specifications, but it is best to identify these and be explicit about them. We return a warning and \code{NA} if more than one age interval is used. It is assumed that ages refer to the lower bounds of age intervals, as is the standard in demography.
#' 
#' @param Dat a \code{data.frame} containing a column called \code{Age}, or code{age}. 
#' @param MinAge integer ignore ages below this age.
#' @param MaxAge integer ignore ages above this age.
#' @param ageColumn character string giving the name of the Age column \code{"Age"} assumed.
#' 
#' @return integer the age interval. \code{NA} if this is not unique.
#' @export

detectAgeInterval <- function(Dat, MinAge = 5, MaxAge = 70, ageColumn = "Age"){
	
	stopifnot(tolower(ageColumn) %in% tolower(colnames(Dat)))
	colnames(Dat)[grepl(tolower(ageColumn),tolower(colnames(Dat)))] <- "Age"

	Ages <- with(Dat, unique(Age[Age >= MinAge & Age <= MaxAge]))
	Interval <- unique(diff(sort(Ages)))
	if (length(Interval) > 1){
		warning("You have more than one interval here!")
		return(NA)
	}
	Interval
}


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

#' @title a cheap way to choose which column to assign a \code{NoteCode} to
#' 
#' @description One property of the LexisDB scripts that might be useful for downstream checks is the ability to trace which functions have modified a given data object. These can go into NoteCode slots. This function writes \code{code} to the first unoccupied \code{NoteCode} column. If all three \code{NoteCode} columns are occupied, it concatenates the end of the third column. This way we preserve a full history. Unfortunately it gets split between columns. Oh well. Good for eyeballing. This function written for the sake of modularity. Function copied from Human Mortality Database collection directly as-is.
#' 
#' @param X the HMD data object that presumably has three \code{NoteCode} columns
#' @param code character string to assign to the column, typically the name of the function operating on \code{X}.
#' 
#' @export


assignNoteCode <- function(X, code){
	
	Free <- colSums(is.na(as.matrix(X[,c("NoteCode1","NoteCode2","NoteCode3")]))) > 0
	if (any(Free)){
		NoteCol <- paste0("NoteCode",min(which(Free)))
		X[[NoteCol]] <- code
	} else {
		X$NoteCode3 <- paste(X$NoteCode3, code, sep = " ")
	}
	X
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



#' @title if necessary divide deaths by intercensal interval
#' @description ideally \code{deaths} is the average annual deaths in the intercensal period, but it is also common to give it as the sum. If this was the case, set \code{deaths.summed} to \code{TRUE} and we take care of it.
#' @param codi the standard object as described in e.g. \code{ggb()}.
#' @param deaths.summed logical. If \code{TRUE} then \code{deaths} was specified as the sum over the intercensal period. Otherwise it was the mean.
#' @return codi a new column, \code{deathsAvg} will be appended.
#' 
#' @export

avgDeaths <- function(codi, deaths.summed = FALSE){
	if (!"deathsAvg" %in% colnames(codi)){
		if (deaths.summed){
			dif                    <- yint2(codi)
			codi$deathsAvg         <- codi$deaths / dif
		} else {
			codi$deathsAvg <- codi$death
		}
	}
	codi
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

