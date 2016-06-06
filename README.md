AdultCoverage
===============

This repository contains R code for a technical paper in progress, provisionally titled "R implementations of three growth balance methods for estimating adult mortality coverage", with [Everton Lima ](http://www.nepo.unicamp.br/nepo/perfils/everton_lima.html) and [Bernardo Queiroz](https://sites.google.com/site/blanza/). It is likely too early to cite, however you are free to see what we're up to and use (with attribution):

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">"R implementations of three growth balance methods for estimating adult mortality coverage</span> by Everton Lima, Bernardo Queiroz, and <a xmlns:cc="http://creativecommons.org/ns#" href="https://sites.google.com/site/timriffepersonal/" property="cc:attributionName" rel="cc:attributionURL">Timothy L. M. Riffe</a> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

More instructions, and possibly a website, will be made available after publication. Aiming for 100% reproducibility here.


DDM R package
=============
This project has produced a small R package that implements three methods for indirect estimation of death registration coverage. It is still not fully documented, and the methods are still not beta tested. So please use and test away. Package test data coming soon.

A short tutorial 
------------------

To download the most recent version of DDM:

Download the [zip ball](https://github.com/timriffe/LifeTable/zipball/master) or [tar ball](https://github.com/timriffe/LifeTable/tarball/master), decompress and run `R CMD INSTALL` on it in the terminal command line, or use the **devtools** package to install the development version:

```r
# install.packages("devtools")

library(devtools)
install_github("timriffe/AdultCoverage/AdultCoverage/R/DDM")
```

Then you can load the package using:

```r
library(DDM)
```

Your data need to be in this kind of shape:

```r
head(Moz)
cod    pop1    pop2 death Age Sex year1 year2
  1 1388350 1963660 88248   0   f  1997  2007
  1 1113675 1615244 11424   5   f  1997  2007
  1  878429 1183939  5677  10   f  1997  2007
  1  854078  991323  6123  15   f  1997  2007
  1  827614  986526  7280  20   f  1997  2007
  1  654465  841416  7212  25   f  1997  2007
```

Here `cod` indicates the group, a single year, sex, region of data that is to be tested. `pop1` and `pop2` are the first and second census, respectively. `deaths` should contain the average number of deaths in each age group in the intercensal period. This could be a straight arithmetic average, or simply the average of the deaths observed around census 1 or census 2. `Age` should be the lower bound of five-year age groups (incl. age 0-4!). `Sex` is character, either `"f"` or `"m"`. Census dates can be conveyed in a variety of ways. If only `year1` and  `year2` are given, we assume Jan 1. It is best to specify proper date classes and use `date1`, `date2` as column names instead:

```r
cod    pop1    pop2 death Age Sex      date1      date2
  1 1388350 1963660 88248   0   f 1997-08-01 2007-08-01
  1 1113675 1615244 11424   5   f 1997-08-01 2007-08-01
  1  878429 1183939  5677  10   f 1997-08-01 2007-08-01
  1  854078  991323  6123  15   f 1997-08-01 2007-08-01
  1  827614  986526  7280  20   f 1997-08-01 2007-08-01
  1  654465  841416  7212  25   f 1997-08-01 2007-08-01
```

Results are contingent on evaluating results for particular age ranges. In spreadsheets this was typically done visually, which a plot referenced to some cell range that the user could manipulate. Here, we have a function that works similarly, but you need to use it just for one data grouping at a time (`cod`):

```r
my_ages <- ggbChooseAges(x[x$cod==1,])
```

This will open a graphics device, where you can interactively select age ranges by clicking on ages. When you are done, click in the margin to close the device, and it returns the vector of ages. You can use these, or any other vector of ages, to manually specify the age range that each method should use:

```r
ggb(Moz, exact.ages = my_ages)
bh1(Moz, exact.ages = my_ages)
bh2(Moz, exact.ages = my_ages)
```

By default these functions will pick a decent age-range on their own:

```r
ggb(Moz)
bh1(Moz)
bh2(Moz)
```

And the result will depend on the age-range chosen. If left to automatically choose age-ranges, the evaluation methods will pick one independently for each data grouping (`cod`). Let's say your data has a large number of groupings (regions, countries, intercensal periods, whatever). You can get a messy overview of results by running:

```r
Results <- ddm(my.huge.data)
ddmplot(Results)
```

This overview plot also gives the harmonic mean of the coverage estimate given from the three methods provided.

What's missing?
==============

Testing. More documentation. A dataset to provide with the package. A proper vignette. References to the papers where these methods come from (Brass, Bennett-Horiuchi, Hill, etc). We aim to provide a manuscript as an overview of the methods provided here, and once-tested, we'll upload the `DDM` package to CRAN. 
