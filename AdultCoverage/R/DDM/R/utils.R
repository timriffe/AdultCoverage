
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

