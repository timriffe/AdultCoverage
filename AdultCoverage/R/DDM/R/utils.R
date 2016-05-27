
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


args(grepl)
detectAgeInterval <- function(Dat, MinAge = 5, MaxAge = 70, ageColumn = "Age"){
	
	stopifnot("age" %in% towlower(colnames(Dat)))
	colnames(Dat)[grepl("age",colnames(Dat))] <- "Age"
if ()
	Ages <- with(Dat, unique(Age[Age >= MinAge & Age <= MaxAge]))
	Interval <- unique(diff(sort(Ages)))
	if (length(Interval) > 1){
		warning("You have more than one interval here!")
		return(NA)
	}
	Interval
}
