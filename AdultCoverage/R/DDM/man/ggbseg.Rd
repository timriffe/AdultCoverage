% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bh.R
\name{ggbseg}
\alias{ggbseg}
\title{estimate death registration coverage using the hybrid generalized growth balance and synthetic extinct generation}
\usage{
ggbseg(
  X,
  minA = 15,
  maxA = 75,
  minAges = 8,
  exact.ages.ggb = NULL,
  exact.ages.seg = NULL,
  eOpen = NULL,
  deaths.summed = FALSE,
  mig.summed = deaths.summed,
  nx.method = 2,
  lm.method = "oldschool",
  opt.method = "r2"
)
}
\arguments{
\item{X}{\code{data.frame} with columns, \code{$pop1}, \code{$pop2}, \code{$deaths}, \code{$date1}, \code{$date2}, \code{$age}, \code{$sex}, and \code{$id} (if there are more than 1 region/sex/intercensal period).}

\item{minA}{the lowest age to be included in search}

\item{maxA}{the highest age to be included in search (the lower bound thereof)}

\item{minAges}{the minimum number of adjacent ages to be used in estimating}

\item{exact.ages.ggb}{optional. A user-specified vector of exact ages to use for coverage estimation in the GGB (first stage) part of the estimation.}

\item{exact.ages.seg}{optional. A user-specified vector of exact ages to use for coverage estimation in the SEG (second stage) part of the estimation.}

\item{eOpen}{optional. A user-specified value for remaining life-expectancy in the open age group.}

\item{deaths.summed}{logical. Is the deaths column given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.}

\item{mig.summed}{logical. Is the (optional) net migration column \code{mig} given as the total per age in the intercensal period (\code{TRUE}). By default we assume \code{FALSE}, i.e. that the average annual was given.}

\item{nx.method}{integer. either 2 or 4. 4 is smoother.}

\item{lm.method}{character, one of:\itemize{
\item{\code{"oldschool"}} default sd ratio operation of still unknown origin
\item{\code{"lm"} or \code{"ols"}} for a simple linear model
\item{\code{"tls"}, \code{"orthogonal"}, or \code{"deming"}} for total least squares
\item{\code{"tukey"}, \code{"resistant"}, or "\code{"median"}} for Tukey's resistant line method
}}

\item{opt.method}{What should we try to minimize when picking ggb age trims? Current options are \code{"RMSE"}, \code{"ORSS"}, \code{"MAE"}, \code{"MAPE"}, or \code{"r2"}. Default \code{"r2"}.}
}
\value{
a \code{data.frame} with columns for the coverage coefficient \code{$coverage}, and the minimum \code{$lower} and maximum \code{$upper} of the age range on which it is based. Rows indicate data partitions, as indicated by the optional \code{$id} variable.
}
\description{
Given two censuses and an average annual number of deaths in each age class between censuses, we can use stable population assumptions to estimate the degree of underregistration of deaths. The method estimates age-specific degrees of coverage. The age pattern of these is assumed to be noisy, so we take the arithmetic mean over some range of ages. One may either specify a particular age-range, or let the age range be determined automatically. If the age-range is found automatically, this is done using the method developed for the generalized growth-balance method. Part of this method relies on a prior value for remaining life expectancy in the open age group. By default, this is estimated using a standard reference to the Coale-Demeny West model life table, although the user may also supply a value. The difference between this method and \code{seg()} is that here we adjust census 1 part way through processing, based on some calculations similar to GGB.
}
\details{
Census dates can be given in a variety of ways: 1) using Date classes, and column names \code{$date1} and \code{$date2} (or an unambiguous character string of the date, like, \code{"1981-05-13"}) or 2) by giving column names \code{"day1","month1","year1","day2","month2","year2"} containing integers. If only \code{year1} and \code{year2} columns are given, then we assume January 1 dates. If year and month are given, then we assume dates on the first of the month. If you want coverage estimates for a variety of intercensal periods/regions/by sex, then stack them, and use a variable called \code{$id} with a unique values for each data chunk. Different values of \code{$id} could indicate sexes, regions, intercensal periods, etc. The \code{$deaths} column should refer to the average annual deaths in each age class in the intercensal period. Sometimes one uses the arithmetic average of recorded deaths in each age, or simply the average of the deaths around the time of census 1 and census 2. To identify an age-range for the GGB method in the traditional visual way, see \code{ggbChooseAges()}, when working with a single year/sex/region of data. A similar visual picker for SEG is not yet implemented. The automatic age-range determination feature of this function tries to implement an intuitive way of picking ages that follows the advice typically given for doing so visually. We minimize the square of the average squared error (RMSE) between the fitted line and right term for the GGB (first stage), and the SEG stage minimizes RMSE for \code{$Cx} against the coverage estimate implied by the age range. Finally, only specify \code{eOpen} when working with a single region/sex/period of data, otherwise the same value will be passed in irrespective of mortality and sex.
}
\examples{
# The Mozambique data
res <- ggbseg(Moz)
res
# The Brasil data
BM <- ggbseg(BrasilMales)
BF <- ggbseg(BrasilFemales)
head(BM)
head(BF)
}
\references{
Hill K. Methods for measuring adult mortality in developing countries: a comparative review. The global burden of disease 2000 in aging populations. Research paper; No. 01.13; 2001.

Hill K, You D, Choi Y.  Death distribution methods for estimating adult mortality: sensitivity analysis with simulated data errors. Demographic Research. 2009; 21:235-254.

Preston, S. H., Coale, A. J., Trussel, J. & Maxine, W. Estimating the completeness of reporting of adult deaths in populations that are approximately stable.  Population Studies, 1980; v.4: 179-202
}
