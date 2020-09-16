#---------------------------------------
# datasets

#' Example data for Brasil females by federal states, years 1991 to 2000
#'
#' A dataset containing 486 rows and 7 variables: Population counts for 1991 and 2000 in abridged ages 0, 1, 5, ... 75, with an open age of 80. Deaths are given as the average death count per age group over the intercensal period. In total there are 53 states in this dataset.
#'
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{id}{integer an id number for each state}
#'   \item{pop1}{integer the census population count in 1991}
#'   \item{pop2}{integer the census population count in 2000}
#'   \item{deaths}{numeric average deaths between censuses}
#'   \item{date1}{integer 1991}
#'   \item{date2}{integer 2000}
#'   \item{age}{integer lower age bound for each age group}
#'   \item{sex}{character, 'f'}
#' }
#' 
#' @source data downloaded from DATASUS \code{http://www.datasus.gov.br}
"BrasilFemales"


#' Example data for Brasil males by federal states, years 1980 to 1991
#'
#' A dataset containing 486 rows and 7 variables: Population counts for 1980 and 1991 in abridged ages 0, 1, 5, ... 75, with an open age of 80. Deaths are given as the average death count per age group over the intercensal period. In total there are 53 states in this dataset.
#' 
#' @format A data frame with 53940 rows and 10 variables:
#' \describe{
#'   \item{id}{integer an id number for each state}
#'   \item{pop1}{integer the census population count in 1991}
#'   \item{pop2}{integer the census population count in 2000}
#'   \item{deaths}{numeric average deaths between censuses}
#'   \item{date1}{integer 1991}
#'   \item{date2}{integer 2000}
#'   \item{age}{integer lower age bound for each age group}
#'   \item{sex}{character, 'm'}
#' }
#' 
#' @source data downloaded from DATASUS \code{http://www.datasus.gov.br}
"BrasilMales"



#' Example data for Mozambique females 1997-2007
#'
#' A dataset containing 17 rows and 8 variables: Population counts for 1997 and 2007 in quinquennial age groups 0, 5, ... 75, with an open age of 80. Deaths are given as the average of the age-specific deaths in 1997 and 2007.
#'
#' @format A data frame with 17 rows and 8 variables:
#' \describe{
#'   \item{id}{integer a column of 1s}
#'   \item{sex}{character ``f'' for female}
#'   \item{pop1}{integer the census population count in 1997}
#'   \item{pop2}{integer the census population count in 2007}
#'   \item{deaths}{integer average of 1997 and 2007 deaths}
#'   \item{age}{integer lower age bound for each age group}
#'   \item{date1}{integer 1997}
#'   \item{date2}{integer 2007}
#' }
#' 
#' @source Data courtesy of Bernardo Queiroz.
"Moz"

#' Example data for South African Males 1997-2007
#'
#' A dataset containing 18 rows and 8 variables. This includes male population counts from censuses in 2001 and 2007, total deaths between censuses, and total net migration between censuses, all by age.
#'
#' @format A data frame with 17 rows and 8 variables:
#' \describe{
#'   \item{id}{integer a column of 1s}
#'   \item{sex}{character ``m'' for male}
#'   \item{pop1}{integer the census population count by age on 2001-10-10}
#'   \item{pop2}{integer the census population count by age on 2007-02-15}
#'   \item{deaths}{integer total deaths by age in the intercensal period}
#'   \item{mig}{integer total net migration by age in the intercensal period}
#'   \item{age}{integer lower age bound for each age group}
#'   \item{date1}{Date class, 2001-10-10}
#'   \item{date2}{Date class, 2007-02-15}
#' }
#' 
#' @source Data from IUSSP tools for demographic estimation GGB or SEG templates. <http://demographicestimation.iussp.org/sites/demographicestimation.iussp.org/files/AM_GGB_South%20Africa_males_4_0.xlsx>
#' 
"ZA"