#' calculate and plot ggb results for each possible age trim
#' 
#' @description Age trim selection is a big open problem for death distribution methods. Here we can see how the `"Mxcoverage"` estiamte varies as a function of lower and upper ages of age trims. Alternatively, one could check other optimal criteria, such as r2, RMS, or ORSS statistics over the same space.
#' 
#' @details The optimal (using the specified `opt.method` criterion) age trim is indicated with a circle. Whether or not `plot` is `TRUE` this function returns a data object that one may wish to analyze.
#' 
#' @inheritParams ggb
#' @param type either `"Mxcoverage"` or `"resid"`
#' @param plot logical. Shall we also render the plot?
#' @import ggplot2
#' @importFrom colorspace scale_fill_continuous_diverging
#' @importFrom colorspace scale_fill_continuous_sequential
#' @export
#' 
#' @examples 
#' \dontrun{
#' }
ggbAgeTrimSensitivity <-
  function(X, 
           minA =5, 
           maxA = 75,
           minAges = 8,
           deaths.summed = FALSE,
           mig.summed = deaths.summed,
           lm.method = "tukey", 
           nx.method = 2,
           opt.method = "r2",
           type = "Mxcoverage",
           plot = TRUE){
    
    codi <- X %>% 
      ggbMakeColumns(minA = minA,
                     maxA = maxA,
                     deaths.summed = deaths.summed,
                     mig.summed = mig.summed,
                     nx.method = nx.method)
    maxAges   <- sum(codi$keep)
    agesUniv  <- codi$age[codi$keep]
    
    FirstAges <- agesUniv[agesUniv < 30]
    
    
    # Create list of all possible age trims
    ind       <- 0
    agesL     <- list()
    # determine ages to test
    for (Nr in maxAges:minAges){ #
      its <- length(agesUniv) - Nr + 1
      for (set in 1:its){ # 
        ind <- ind + 1
        agesL[[ind]] <- agesUniv[set:(set+Nr-1)]
      }
    }
    
    # some reference containers
    res   <- rep(NA,length(agesL))
    mat   <- lapply(agesL,range) %>% do.call("rbind",.)
    colnames(mat) <- c("from","to")
    mat   <- mat %>% 
      as_tibble() %>% 
      mutate(value = NA,
             optimum = FALSE)
    res <- rep(NA,nrow(mat))
    for (i in 1:length(res)){
      mn <- ggb(ZA,
                exact.ages = agesL[[i]],
                deaths.summed =  deaths.summed,
                mig.summed = mig.summed,
                lm.method = lm.method,
                nx.method = nx.method)
      
      if (type == "Mxcoverage"){
        res[i]  <- mn$Mxcoverage
      }
      if (type == "resid"){
        res[i]  <- 
          ggbgetRMS(agesL[[i]],
                    codi,
                    lm.method = lm.method,
                    opt.method = opt.method,
                    scale = 1)
      }
    }
    
    mn <- ggb(ZA,
              minA = minA,
              maxA = maxA,
              minAges = minAges,
              deaths.summed = deaths.summed,
              mig.summed = mig.summed,
              lm.method = lm.method,
              nx.method = nx.method,
              opt.method = opt.method)
    
    sen <-
      mat %>% 
      mutate(value = res,
             optimum = ifelse((from == mn[["lower"]] & to == mn[["upper"]]),TRUE,FALSE) )
    opt <- sen %>% 
      filter(optimum)
    
    
    if (plot){
      if (type == "Mxcoverage"){
        plotit <-
          sen %>% 
          ggplot(aes(x=from,y=to,fill=value)) + 
          geom_tile()+
          scale_fill_continuous_diverging("Blue-Red3",
                                          trans ="log") +
          geom_point(opt,mapping=aes(x=from,y=to),
                     size = 3,
                     fill="black",
                     color='white',
                     pch = 21) + 
          labs(title = "Mxcoverage estimate by age trim",
               subtitle = paste0("opt.method = ",opt.method,"; lm.method = ",lm.method))
        print(plotit)
      }
      if (type == "resid"){
        plotit <- 
          sen %>% 
          ggplot(aes(x=from,y=to,fill=-value)) + 
          geom_tile()+
          scale_fill_continuous_sequential("Burg")+
          geom_point(opt,mapping=aes(x=from,y=to),
                     size = 3,
                     fill="black",
                     color='white',
                     pch = 21) + 
          labs(title = paste0(opt.method, " estimate by age trim"),
               subtitle = paste0("opt.method = ",opt.method,"; lm.method = ",lm.method))
        print(plotit)
      }
    }
    invisible(sen)
  }



# functions related to ggb plotting.


#' @title does a given pairlist of x and y coordinates fall within the plot region?
#' 
#' @description Check to see if a point clicked falls in the plot or outside it. This function is used by \code{ggbChooseAges()}.
#' 
#' @param USR, as given by \code{par("usr")}
#' @param click a pairlist with elements \code{$x} and \code{$y}, as returned by \code{locator(1)}
#' 
#' @return logical. \code{TRUE} if in the plot region.
#' 
#' @export
inUSR <- function(USR, click){
  x <- click$x
  y <- click$y
  if (x < USR[1] | x > USR[2] | y < USR[3] | y > USR[4]){
    return(FALSE)
  }
  TRUE
}


#' @title which age is closest to the point clicked?
#' @description a utility function called by \code{ggbChooseAges()}.
#' 
#' @param xvec \code{$rightterm}, as given by \code{ggbMakeColumns()}
#' @param yvec \code{$lefttterm}, as given by \code{ggbMakeColumns()}
#' @param click a point given by \code{locator(1)}
#' @param age ages present in dataset
#' 
#' @return the age corresponding to the x,y pair of \code{$rightterm}, \code{$lefttterm} closest to the point clicked.
#' 
#' @export
guessage <- function(xvec,yvec,click,age){
  age[which.min(sqrt((xvec - click$x) ^ 2 + (yvec - click$y) ^ 2))]
}


#' @title adjust the range of ages used by \code{ggbChooseAges()}
#' @description a utility function called by \code{ggbChooseAges()}. After clicking a point, this function readjusts the age range
#' 
#' @param a an age specified by the user, as returned by \code{guessage()}
#' @param age ages present in dataset
#' @param agesfit the age range used for calculating the coverage coefficient
#' 
#' @return the adjusted set of ages used for calculating the coverage coefficient
#' 
#' @export
adjustages <- function(a, age, agesfit){
  maxa <- max(agesfit)
  mina <- min(agesfit)
  
  isitamin <- abs(a - mina) < abs(a - maxa)
  if (isitamin){
    mina <- a
  } else {
    maxa <- a
  }
  age[age >= mina & age <= maxa]
}

#' @title get the slope the slope and intercept implied by a set of ages
#' @description Called by \code{ggbFittedFromAges()} and \code{ggbChooseAges()}
#' @param codi \code{data.frame} as produced by \code{ggbMakeColumns()}
#' @param agesfit a set of continuous ages to estimate coverage from 
#' @param lm.method character, one of:\itemize{
#'   \item{\code{"oldschool"}} default sd ratio operation of still unknown origin
#'   \item{\code{"lm"} or \code{"ols"}} for a simple linear model
#'   \item{\code{"tls"}, \code{"orthogonal"}, or \code{"deming"}} for total least squares
#'   \item{\code{"tukey"}, \code{"resistant"}, or "\code{"median"}} for Tukey's resistant line method
#' }
#' @return a pairlist with elements \code{$a} for the intercept and \code{$b} for the slope
#' @importFrom tukeyedar eda_rline
#' @importFrom stats lm
#' @importFrom stats cov
#' @importFrom stats var
#' @importFrom stats sd
#' @export
slopeint <- function(codi, 
                     agesfit, 
                     lm.method = "oldschool"){
  
  lm.method <- match.arg(lm.method, choices = c("oldschool", "lm", "ols", "tls","orthogonal","deming","tukey","resistant","median"))
  stopifnot("leftterm" %in% colnames(codi))
  
  codi <- codi %>% 
    filter(.data$age %in% agesfit)
  
  if (lm.method == "oldschool"){
    #age <- codi$age
    # TODO: find eq numbers to cite here
    slope       <- 	sd(codi$leftterm) /  sd(codi$rightterm)
    
    # PJ: fix https://github.com/timriffe/AdultCoverage/issues/3
    intercept <- mean(codi$leftterm) - mean(codi$rightterm) * slope
    
    #coefs <- with(codi,lm(leftterm~rightterm))$coef
  }
  
  
  if (lm.method %in% c("lm","ols")){
    
    ab        <- lm("leftterm~rightterm", data = codi)$coef
    slope     <- ab[2]
    intercept <- ab[1]
  }
  
  if (lm.method %in% c("tls","orthogonal","deming")){
    y <- codi$leftterm
    x <- codi$rightterm
    delta <- 1 # assumed
    # let's avoid a new package dependency
    # PJ's implementation via Wikipedia formulas
    slope <- (var(y)-delta*var(x)+sqrt((var(y)-delta*var(x))^2+4*delta*cov(x,y)^2)) /
      (2*cov(x,y))
    
    intercept <- mean(y) - slope * mean(x)
  }
  
  if (lm.method %in% c("tukey","resistant","median")){
    ab.etc <- eda_rline(x = codi$rightterm, y = codi$leftterm)
    slope     <- ab.etc$b
    intercept <- ab.etc$a
  }
  
  # PG: return r2 also
  list(a = intercept, b = slope)
}



#' @title interactively determine ages to use for estimating coverage
#' @description In a spreadsheet one would typically set up the GGB method to produce a plot that updates as the user changes the age range. This function implements that kind of work flow. This will be intuitive for spreadsheet users, but it does not scale well. Imagine you have 200 territorial units, then you would not want to repeat this task. \code{ggb()} does the same thing automatically. You can compare the age range you select manually with the one given back by \code{ggb()} as a diagnostic, for instance. To set up the plot device, just give a single year/region/sex of data. By default it will give the RMSE-optimized age range to start with, but you can specify a  vector of exact ages to use as well. All points are plotted, with a fitted line that has been set to a subset of the points, which is plotted in a different color. You can click any point to change the age range, and the plot updates accordingly, up to a maximum of 15 clicks so you don't waste your time. You can stop the plot by either clicking on the graphics device outside the plot area or clicking out the 15 tries (or more if you increase \code{maxit}).
#' @details If you want to send the results of this into \code{ggb()}, you can do so by setting \code{Exact.ages} to \code{seq(lower,upper,by=5)}, where \code{$lower}, and \code{$upper} are the results returned from \code{ggbChooseAges()} after you're done manually determining the age range.
#' 
#' @inheritParams ggb
#' @param maxit up to how many times do you want to let yourself fiddle with the age range?
#' @return \code{data.frame} containing elements \code{$coverage}, \code{$lower}, \code{$upper}, and \code{ages}.
#' 
#' @importFrom grDevices gray
#' @importFrom graphics abline legend locator mtext par plot points rect segments text
#' @export
#' 
#' @examples
#' \dontrun{
#' # for interactive sessions only
#' # *click points to adjus age range used (yellow)
#' # *click in margin to stop and return coverage results
#' ggbChooseAges(Moz)
#' }

ggbChooseAges <- function(
  X, 
  minA = 15, 
  maxA = 75, 
  minAges = 8, 
  exact.ages = NULL, 
  maxit = 15, 
  deaths.summed = FALSE,
  mig.summed = deaths.summed,
  lm.method = "oldschool",
  nx.method = 2){
  # this is the automatic age selection.
  
  # only run if in anteractive r session...
  stopifnot(interactive())
  
  # reset ages if necessary
  if (!is.null(exact.ages) & length(exact.ages) >= 3){
    if (min(exact.ages) < minA){
      minA <- min(exact.ages)
    } 
    if (max(exact.ages) > maxA){
      maxA <- max(exact.ages)
    }
    if (minAges < length(exact.ages)){
      minAges <- length(exact.ages)
    }
  }
  
  codi_in <- X
  codi    <- data.frame(X)           
  
  
  # guess which column is the deaths column, rename it deaths
  # codi    <- guessDeathsColumn(codi)
  
  # start GGB stuff
  codi    <- ggbMakeColumns(codi, 
                            minA = minA, 
                            maxA = maxA, 
                            nx.method = nx.method,
                            deaths.summed = deaths.summed,
                            mig.summed = mig.summed)
  
  # some potential starting ages. either auto or self-supplied
  if (!is.null(exact.ages) & length(exact.ages) >= 3){
    agesfit <- exact.ages
  } else {
    fit.res <- ggbgetAgesFit(codi, 
                             minA = minA, 
                             maxA = maxA, 
                             minAges = minAges,
                             lm.method = lm.method)
    agesfit <- fit.res$agesfit
  }
  
  codi     <- ggbFittedFromAges(codi = codi, 
                                agesfit = agesfit, 
                                lm.method = lm.method)
  
  
  coverage <- ggbcoverageFromYear(codi = codi, 
                                  exact.ages = agesfit, 
                                  minA = minA, 
                                  maxA = maxA,
                                  minAges = minAges,
                                  deaths.summed = deaths.summed,
                                  mig.summed = mig.summed,
                                  lm.method = lm.method,
                                  nx.method = nx.method)
  
  
  age      <- codi$age
  # rm border ages re PJ issue #2
  ages.rm  <- age < minA | age > maxA 
  leftt    <- codi$leftterm
  rightt   <- codi$rightterm
  leftt[ages.rm]  <- NA
  rightt[ages.rm] <- NA
  
  
  # age ranges used for fitting
  amin    <- min(agesfit); amax <- max(agesfit)
  plot(rightt, leftt, asp = 1, pch = 19, col = "#00000050", cex = 1.6,
       xlab = "right term",
       ylab = "left term",
       main = paste0("Age range [",amin,
                     ",",amax,"], \nLine method = ",lm.method,"\nEst. coverage = ",round(coverage$Mxcoverage*100,1)),
       sub = "(optimized age range)")
  # automatically fit line (RMS of ggb)
  abline(a = coverage$a, b = coverage$b, col = "blue")
  # shows points used to fit line
  points(rightt[age %in% agesfit], 
         leftt[age %in% agesfit], col = "#FFFF00", pch = 19, cex = 1.6)
  text(rightt, leftt, age, cex = .6)
  legend("bottomright",lty=1,col="blue",legend="fitted line",bty="n")
  # message to user
  mtext("*Adjusting ages*\nClick any age to change bounds for fitting\nClick in margin to quit",
        side = 3, line = -3, outer = FALSE)
  
  
  USR    <- par()$usr
  agesfit <- agesfit # keep old one...
  for (i in 1:maxit){
    # now enter into interactive loop
    clicki <- locator(1)
    # was the click inside the plot or outside the plot?
    IN     <- inUSR(USR, clicki)
    if (IN){
      
      # if it's inside the plot, then which is the closest age?
      a        <- guessage(rightt,leftt,clicki,age)
      # guess how to adjust ages based on which age was clicked
      agesfit  <- adjustages(a, age, agesfit)
      # new range of ages used for fitting
      amin     <- min(agesfit)
      amax     <- max(agesfit)
      
      coverage <- ggbcoverageFromYear(codi = codi, 
                                      exact.ages = agesfit, 
                                      minA = minA, 
                                      maxA = maxA,
                                      minAges = minAges,
                                      deaths.summed = deaths.summed,
                                      mig.summed = mig.summed,
                                      lm.method = lm.method,
                                      nx.method = nx.method)
      
      # redraw plot
      plot(rightt, 
           leftt, 
           asp = 1, 
           pch = 19, 
           col = "#00000050", 
           cex = 1.6,
           xlab = "right term",
           ylab = "left term",
           main = paste0("Age range [", amin,
                         ",", amax, "], est. coverage = %",round(coverage$Mxcoverage * 100, 1)),
           sub = "(optimized age range)")
      # new fitted slope, intercept
      abline(a=0,b=1,col=gray(.8)) # line of perfection
      #
      abline(a = coverage$a, b = coverage$b, col = "blue")
      # indicate which points used with color
      points(rightt[age %in% agesfit], 
             leftt[age %in% agesfit], col = "#FFFF00", pch = 19, cex = 1.6)
      text(rightt, leftt, age, cex = .6)
      legend("bottomright", lty = 1, col = "blue", legend = "fitted line", bty = "n")
    } else {
      break
    }
  }
  
  # click outside the margin to return results
  out <- ggb(X = codi_in,
             minA = minA, 
             maxA = maxA, 
             minAges = minAges, 
             exact.ages = agesfit, 
             lm.method = lm.method,
             nx.method =  nx.method,
             deaths.summed = deaths.summed,
             mig.summed = deaths.summed)
  out
}

# end
