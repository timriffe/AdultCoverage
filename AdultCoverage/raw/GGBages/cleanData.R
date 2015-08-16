

# this function should receive the input data and sort it into the necessary pieces.

# it should be able to handle a few basic types of data.

head(UFF)
Data <- UFF
cleanData <- function(Data){
	# 1) all colnames to lower:
	colnames(Data) <- tolower(colnames(Data))
	
	# 2) which columns are pop columns:
	popN <- sum(grepl("pop",colnames(Data)))
	# in test case there are four... we should throw an error:
	stopifnot(popN <= 2, "We know how to interpret one or two population columns, but not more. Let's do this one period at a time. Please chop data down.\n")

	# detect age column:
	
}








