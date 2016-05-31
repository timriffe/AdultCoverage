
# Author: tim
###############################################################################
# contains utilities used throught.

#'
#' @title Detect the age interval for some demographic data
#' 
#' @description Since death distribution methods are primarily used in adult ages, it's OK to chop off the irregular infant and child age intervals (0,1], (1,5]. Further, if high ages are in different intervals this might also be a non-issue. In principal, the user should set \code{MinAge} and \code{MaxAge} to the same values used in the death distribution methods. Here we have some defaults that should almost always return the result \code{5} for standard abridged data, or \code{1} for single age data. Really there aren't any other common age-specifications, but it is best to identify these and be explicit about them. We return a warning and \code{NA} if more than one age interval is used. It is assumed that ages refer to the lower bounds of age intervals, as is the standard in demography.
#' 
#' @param Dat a \code{data.frame} containing a column called \code{Age}, or code{age}. 
#' @param MinAge integer ignore ages below this age.
#' @param MinAge integer ignore ages above this age.
#' @return integer the age interval. \code{NA} if this is not unique.
#' @export
#' 

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


#' @title a cheap way to choose which column to assign a NoteCode to
#' 
#' @description One property of the LexisDB scripts that might be useful for downstream checks is the ability to trace which functions have modified a given data object. These can go into NoteCode slots. This function writes \code{code} to the first unoccupied \code{NoteCode} column. If all three \code{NoteCode} columns are occupied, it concatenates the end of the third column. This way we preserve a full history. Unfortunately it gets split between columns. Oh well. Good for eyeballing. This function written for the sake of modularity. Function copied from HMD collection directly as-is.
#' 
#' @param X the HMD data object that presumably has three \code{NoteCode} columns
#' @param code character string to assign to the column, typically the name of the function operating on \code{X}.
#' 
#' @export
#' 

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


#' Logical utility functions
#'
#' @aliases logic HMDutils
#' 
#' @description These logical functions are like the usual ones, but \code{NA} values are treated as \code{FALSE} by default. This is not an exhaustive list, but these are the ones that speed our coding, and reduce code clutter. Functions copied from HMD collection directly as-is.
#' 
#' @param x,y any two vector that can be logically compared.
#' @name HMDlogic
#' 
#' @examples
#' \dontrun{
#' c(1,2,NA,4,5) == c(1,NA,3,4,NA)
#' # compare
#' c(1,2,NA,4,5) %==% c(1,NA,3,4,NA)
#' }
NULL
#' 

#' @rdname HMDlogic
'%==%' <- function(x,y){
	x == y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%!=%' <- function(x,y){
	x != y & !is.na(x) & !is.na(y)
}

# note this is incompatible with magrittr!
#' @rdname HMDlogic
'%>%' <- function(x,y){
	x > y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%<%' <- function(x,y){
	x < y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%>=%' <- function(x,y){
	x >= y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%<=%' <- function(x,y){
	x <= y & !is.na(x) & !is.na(y)
}

#' @title append a \code{$cod} column if missing
#' @description Only handles the case of missing \code{$cod} splitting variable for data of a single year/region. This isn't super robust. If you have many regions or whatever then do it yourself. This function was just written to make \code{ggb()} robust to the case of a user specifying data that don't have any territorial or other subgroups, aside from sex.
#' 
#' @param X a \code{data.frame}, possibly but not necessarily with column \code{$sex}.
#' @return X with a new column, \code{$cod} appended. 
#' 
#' @export
#' 
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

#'
#' @title a function to determine whether a year is a leap year. 
#' 
#' @description In order to remove lubridate dependency, we self-detect leap years and adjust February accordingly.
#' 
#' @param Year integer of year to query
#' 
#' @return logical is the Year a leap year or not
#' 
#' @export
#' @author Carl Boe

isLeapYear <- function (Year){      # CB: mostly good algorithm from wikipedia
	ifelse(
			( (Year %% 4) == 0  &  (Year %% 100) != 0   ) | ( (Year %% 400) == 0 ),
			TRUE, FALSE )
}

#'
#' @title ypart function to determine the proportion of a year passed as of a particular date
#' 
#' @description The fraction returned by this is used e.g. for intercensal estimates. Function uses 'lubridate' package to handle dates elegantly.
#' 
#' @param Year 4-digit year (string or integer)
#' @param Month month digits (string or integer, 1 or 2 characters)
#' @param Day Day of month digits (string or integer, 1 or 2 characters)
#' @param reproduce.matlab logical. Default TRUE. Assume 365 days in a year.
#' @param detect.mid.year logical. if \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' 
#' @export
#' 

ypart <- function(Year, Month, Day, reproduce.matlab = TRUE, detect.mid.year = FALSE, detect.start.end = TRUE){
	M <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
	if (reproduce.matlab){
		Day   <- as.integer(Day)
		Month <- as.integer(Month)
		if (is.na(Day) & is.na(Month)){
			return(.5)
		}
		if (Day == 1 & Month == 1){
			return(0)
		}
		return((M[Month] + Day) / 365)
	}
	# this chunk written just to avoiding writing p_my()
	if (detect.mid.year){
		.d <- as.integer(Day)
		.m <- as.integer(Month)
		if ((.d == 30 & .m == 6) | (.d == 1 & .m == 7)){
			return(.5)
		}
	}
	
	if (detect.start.end){
		.d <- as.integer(Day)
		.m <- as.integer(Month)
		if (.d == 1 & .m == 1){
			return(0)
		}
		if(.d == 31 & .m == 12){
			return(1)
		}
	}
	
	monthdur    <- diff(c(M,365))
	monthdur[2] <- monthdur[2] + isLeapYear(Year)
	M           <- cumsum(monthdur) - 31
	return((M[Month] + Day) / sum(monthdur))
	
	# get into date class
#  Date  <- ymd(paste(Year, Month, Day, sep = "/"))
#  # what was jan 1st?
#  Jan1  <- floor_date(Date, "year") 
#  # how many days have passed?
#  Days  <- yday(Date) - 1
#  # if we want to account for possible leap years, we get next year's Jan 1st (our Dec 31st)
#  Dec31 <- floor_date(ymd(paste(Year + 1, Month, Day, sep = "/")), "year") 
#  # and the difference is the year length
#  Denom <- as.integer(Dec31 - Jan1)
#  # only gives same as matlab on Jan 1st.
#  Days / Denom
}


#'
#' @title yint get interval as fraction of full years
#' 
#' @description Either assume 365 days in the year, or get the precise duration.
#' 
#' @param Day1 Day of first date
#' @param Month1 Month of first date
#' @param Year1 Year of first date
#' @param Day1 Day of second date
#' @param Month1 Month of second date
#' @param Year1 Year of second date
#' @param reproduce.matlab logical. default \code{TRUE}. Assume 365 days in all years?
#' @param detect.mid.year logical. default \code{FALSE}. Should June 30 and July 1 be considered .5?
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' 
#' @return decimal value of year fraction (can be greater than 1)
#' 
#' @export
#' 
yint <- function(Day1, Month1, Year1, Day2, Month2, Year2, reproduce.matlab = TRUE, 
		detect.mid.year = FALSE, detect.start.end = TRUE){
	if (reproduce.matlab){
		return(abs(Year1 - Year2 + (Day1 - Day2) / 365 +  (Month1 - Month2) / 12))
	}
	
	# we can be more exacting, if desired:
#  abs(decimal_date(ymd(paste0(Year1, sprintf("%02d", Month1), sprintf("%02d", Day1)))) -
#      decimal_date(ymd(paste0(Year2, sprintf("%02d", Month2), sprintf("%02d", Day2)))))
#  
	Ypart1 <- ypart(Year = Year1, 
			Month = Month1, 
			Day = Day1, 
			reproduce.matlab = reproduce.matlab, 
			detect.mid.year = detect.mid.year, 
			detect.start.end = detect.start.end)
	Ypart2 <- ypart(Year = Year2, 
			Month = Month2, 
			Day = Day2, 
			reproduce.matlab = reproduce.matlab, 
			detect.mid.year = detect.mid.year, 
			detect.start.end = detect.start.end)
#  
	(1 - Ypart1) + abs(Year2 - Year1) + Ypart2
	
}


#'
#' @title assume Jan 1 if no month or day given
#' @description We still require two year columns, \code{year1} and \code{year2}, at a minimum. If this function is called, and if month and day columns are missing we add these columns, with values of 1. If date columns are given, then these must be either in an unambiguous character format ("YYYY-MM-DD", e.g. \code{"2016-05-30"} is unambiguous). Date columns will override the presence of other year, month, day columns.
#' 
#' @param X a \code{data.frame} with at least columns \code{year1} and \code{year2}.
#' 
#' @return X the same \code{data.frame}, possibly with columns for year, month, or day added.
#' 
#' @export
fakeDates <- function(X){
	tempnames <- tolower(colnames(X))
	
	if (all(c("date1","date2") %in% tempnames)){
		# make sure we have unique dates
		stopifnot(length(unique(X$date1)) == 1 & length(unique(X$date2)) == 1)
		
		if (class(X$date1) == "character"){
			# get date class from character. Must be
			# in unambiguous format! This will throw error 
			# otherwise.
			# "YYYY-MM-DD", e.g. "2016-05-30" is unambiguous.
			X$date1 <- as.Date(X$date1)
			X$date2 <- as.Date(X$date2)
		}
		# one final stopcheck.
		stopifnot(class(X$date1) == "Date" & class(X$date2) == "Date")
		
		X$day1 		<- as.numeric(format(unique(X$date1),'%d'))
		X$month1 	<- as.numeric(format(unique(X$date1),'%m'))
		X$year1 	<- as.numeric(format(unique(X$date1),'%Y'))
		X$day2 		<- as.numeric(format(unique(X$date2),'%d'))
		X$month2 	<- as.numeric(format(unique(X$date2),'%m'))
		X$year2 	<- as.numeric(format(unique(X$date2),'%Y'))
	} else {	
		stopifnot(all(c("year1","year2") %in% tempnames))
		
		if (!"month1" %in% tempnames){
			X$month1 	<- 1
		}
		if (!"month2" %in% tempnames){
			X$month2 	<- 1
		}
		if (!"day1" %in% tempnames){
			X$day1 		<- 1
		}
		if (!"day2" %in% tempnames){
			X$day2 		<- 1
		}
	}
	X
}

#'
#' @title get the time interval without having to specify so many args
#' @description We accept dates, and fake them otherwise. Dates must be unique. Iterate over data if necessary for multiple intervals.
#' 
#' @param X \code{data.frame} with at least columns \code{$date1} and \code{$date2}, or \code{$year1} and \code{$year2}. 
#' 
#' @return an decimal year value of the time between two dates.
#' 
#' @export
yint2 <- function(X){
	colnames(X) <- tolower(colnames(X))
	X           <- fakeDates(X)
	yint(   Day1 = unique(X$day1), 
			Month1 = unique(X$month1),
			Year1 = unique(X$year1), 
			Day2 = unique(X$day2), 
			Month2 = unique(X$month2), 
			Year2 = unique(X$year2), 
			reproduce.matlab = FALSE, 
			detect.mid.year = TRUE, 
			detect.start.end = TRUE)
	
}


#'
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