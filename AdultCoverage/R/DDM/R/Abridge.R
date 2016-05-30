
# Author: tim
###############################################################################

single2abr <- function(x){
	stopifnot(is.integer(x))
	x5 <- x - x %% 5
	x5[x %in% c(1:4)] <- 1
	x5
}

#' @title Abridge population or death data 
#' @description This function assumes you have columns named Age, Sex, AgeInterval, and Year
#' @importFrom data.table data.table
Abridge <- function(X, value.var = "Deaths"){
	# assume we have Age, AgeInterval, etc.
	
	X <- as.data.table(X)
	
	# this doesn't do regions...
	X[,Age5 := single2abr(Age), by = c(Year, Sex)]
	
	#
	# test: value.var might not pass properly, need live data
	X[,list(V1=sum(value.var),AgeInterval5 = sum(AgeInterval)), by = c(Age5, Year, Sex)]
	
	# rename new abridged column
	setnames(X,"V1", value.var)
	setnames(X,"Age5", "Age")
	
	
	# that's that!
	X <- as.data.frame(X)
	
	X
}

