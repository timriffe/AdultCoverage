
#ITA <- read.table("Data/italy.test.txt", 
#  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#
##Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)
#
## new data for Brasil to run. By regions.
#BR1 <- read.table(file.path("Data","data_Brazil_p1.txt"), 
#		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#BR2 <- read.table(file.path("Data","data_Brazil_p2.txt"), 
#		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#BR3 <- read.table(file.path("Data","data_Brazil_p3.txt"), 
#		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#
library(devtools)

setwd("AdultCoverage/R/DDM")

load_all()
document()
  check()
ddm(Moz,delta=TRUE)
ddm(BrasilMales)
ggb(ZA,
      deaths.summed=TRUE,
      mig.summed =TRUE,
      exact.ages=seq(5,80,by=5),
      lm.method = "oldschool")

load_all()
ggbMakeColumns(ZA) %>% View()
  
seg(ZA,
    deaths.summed=TRUE,
    mig.summed =TRUE,
    exact.ages=seq(25,60,by=5),
    exact.ages.ggb=seq(5,80,by=5),
    delta=TRUE,
    lm.method="oldschool",
    eOpen=4.34676243530911)


segMakeColumns(ZA,
               deaths.summed=TRUE,
               mig.summed =TRUE,
               exact.ages.ggb=seq(5,80,by=5),
               delta=TRUE,
               lm.method="ols",
               eOpen=4.35) %>% 
  View()

segCoverageFromYear(ZA,
               deaths.summed=TRUE,
               mig.summed =TRUE,
               exact.ages.ggb=seq(5,80,by=5),
               delta=TRUE,
               lm.method="tls",
               eOpen=4.35)$codi %>% View()

seg(ZA,
    deaths.summed=TRUE,
    mig.summed =TRUE,
    exact.ages.ggb=seq(5,80,by=5),
    delta=TRUE,
    lm.method="tls",
    eOpen=4.35)

#--------------------
library(dplyr)
library(magrittr)

# change data
ggb(Moz)
seg(Moz,delta=F)
ggbseg(Moz)
res <- ddm(Moz,delta=TRUE)

codi <- Moz %>% ggbMakeColumns() 
ggbcoverageFromYear(codi,lm.method = "ols")
ZA <- readr::read_csv("/home/tim/workspace/AdultCoverage/AdultCoverage/Data/AM_SEG_South Africa_males_3_1.csv")
ZA <- 
ZA %>% 
  mutate(date1 = mdy(date1),
         date2 = mdy(date2))
ggb(ZA, exact.ages=seq(5,80,by=5),lm.method="deming", deaths.summed = TRUE, mig.summed=TRUE)
Moz$id<-1
ggb(Moz, exact.ages=seq(5,75,by=5),lm.method="deming", deaths.summed = FALSE, mig.summed=FALSE)

seg(ZA,maxA=80,minA=15,exact.ages=seq(5,80,by=5), deaths.summed = TRUE, mig.summed=TRUE, eOpen=4.35,exact.ages.ggb=seq(5,80,y=5),delta=TRUE)

seg(Moz, exact.ages=seq(5,75,by=5), deaths.summed = FALSE, mig.summed=FALSE)
seg(ZA)
ggbseg(ZA)
library(DependenciesGraphs)
# Prepare data
dep <- funDependencies("package:DDM","ggbgetRMS")

# visualization
plot(dep)

dep <- envirDependencies("package:DDM")

# visualization
plot(dep)









minA = 15 
maxA = 75 
minAges = 8 
exact.ages.ggb = NULL
exact.ages.seg = NULL 
eOpen = NULL
deaths.summed = FALSE
mig.summed = deaths.summed
lm.method = "tukey"
nx.method = 2
delta = FALSE
opt.method="logRMS"

r_rms <- ggb(BrasilMales, 
             minA=5,
             maxA=75,
             lm.method="tukey",
             opt.method = "RMS",
             minAges=10) %>% 
  mutate(rg = upper - lower) %>% 
  select(upper,lower,rg)

# r_lrms <- ggb(BrasilMales, 
#              minA=5,
#              maxA=75,
#              lm.method="tukey",
#              opt.method = "logRMS",
#              minAges=10) %>% 
#   mutate(rg = upper - lower) %>% 
#   select(upper,lower,rg)
r_ors <- ggb(BrasilMales, 
             minA=5,
             maxA=75,
             lm.method="tukey",
             opt.method = "ORSS",
             minAges=10) %>% 
  mutate(rg = upper - lower) %>% 
  select(upper,lower,rg)
r_hrs <- ggb(BrasilMales, 
             minA=5,
             maxA=75,
             lm.method="tukey",
             opt.method = "hybrid",
             minAges=10) %>% 
  mutate(rg = upper - lower) %>% 
  select(upper,lower,rg)
# r_lors <- ggb(BrasilMales, 
#              minA=5,
#              maxA=75,
#              lm.method="tukey",
#              opt.method = "logORSS",
#              minAges=10) %>% 
#   mutate(rg = upper - lower) %>% 
#   pull(rg) %>% 
#   mean()

codi <- ZA %>% ggbMakeColumns(deaths.summed=TRUE,minA=5,maxA=80)
maxAges   <- sum(codi$keep)
minAges=8
agesUniv  <- codi$age[codi$keep]

FirstAges <- agesUniv[agesUniv < 30]

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

res <- rep(NA,length(agesL))
res2 <- rep(NA,length(agesL))
res3 <- rep(NA,length(agesL))
res4 <- rep(NA,length(agesL))
for (i in 1:length(res)){
  scl<- sqrt(diff(range(agesL[[i]])))
  res[i] <-
    ggb(ZA,
        deaths.summed = TRUE,
        exact.ages = agesL[[i]],
        lm.method = "tukey")$Mxcoverage
  res2[i] <- ggbgetRMS(agesL[[i]],
                       ggbMakeColumns(ZA,deaths.summed=TRUE),
                       lm.method = "tukey",
                       opt.method = "RMS",
                       scale = scl)
  res3[i] <- ggbgetRMS(agesL[[i]],
                       ggbMakeColumns(ZA,deaths.summed=TRUE),
                       lm.method = "tukey",
                       opt.method = "ORSS",
                       scale = scl)
  res4[i] <- ggbgetRMS(agesL[[i]],
                       ggbMakeColumns(ZA,deaths.summed=TRUE),
                       lm.method = "tukey",
                       opt.method = "r2",
                       scale = scl)
}

mat <- lapply(agesL,range) %>% do.call("rbind",.)

froms <- unique(mat[,1]) %>% sort()
tos <- unique(mat[,2]) %>% sort()

MAT  <- matrix(NA,length(froms),length(tos),dimnames=list(from=froms,to=tos))
MAT4 <- MAT3 <- MAT2 <- MAT
for (i in 1:nrow(mat)){
  fc          <- as.character(mat[i,1])
  tc          <- as.character(mat[i,2])
  MAT[fc,tc]  <- res[i]
  MAT2[fc,tc] <- res2[i]
  MAT3[fc,tc] <- res3[i]
  MAT4[fc,tc] <- res4[i]
}
NAmat <- 
rgs           <- outer(rev(tos),rev(froms),"-") + 5
dimnames(rgs) <- dimnames(MAT)

image(froms,tos,MAT,asp=1,xlab="lower",ylab="upper")
contour(froms,tos,MAT,add=TRUE)

image(froms,tos,MAT2,asp=1,xlab="lower",ylab="upper")
contour(froms,tos,MAT2,add=TRUE)

image(froms,tos,MAT3,asp=1,xlab="lower",ylab="upper")
contour(froms,tos,MAT3,add=TRUE)

image(froms,tos,MAT4,asp=1,xlab="lower",ylab="upper")
contour(froms,tos,MAT4,add=TRUE)


image(froms,tos,MAT2,asp=1,xlab="lower",ylab="upper")
image(froms,tos,MAT2/sqrt(rgs),asp=1,xlab="lower",ylab="upper")
image(froms,tos,MAT2/rgs,asp=1,xlab="lower",ylab="upper")
image(froms,tos,MAT2/(rgs^2),asp=1,xlab="lower",ylab="upper")


contour(froms,tos,MAT2,add=TRUE)

image(froms,tos,MAT3/sqrt(rgs),asp=1,xlab="lower",ylab="upper")
contour(froms,tos,MAT3,add=TRUE)

image(froms,tos,MAT4/sqrt(rgs),asp=1,xlab="lower",ylab="upper")
contour(froms,tos,MAT4,add=TRUE)




ggb(ZA,
    deaths.summed = TRUE,
    minA=5,maxA=80,
    lm.method = "tukey",
    opt.method="hybrid")

ggb(ZA,
    deaths.summed = TRUE,
    minA=5,maxA=80,
    lm.method = "tukey",
    opt.method="r2")

ggb(Moz,
    deaths.summed = F,
    minA=5,maxA=80,
    lm.method = "tukey",
    opt.method="r2")

BMtukeyr2 <- 
ggb(BrasilMales,
    minA=5,maxA=75,
    lm.method="tukey",
    opt.method = "r2")

combos <- expand.grid(lm.method = c("oldschool","ols","tls","tukey"),
            opt.method = c("RMS","ORSS","r2"),stringsAsFactors = FALSE)

testBR <- function(.combos,BR){
  .combos <- unlist(.combos)
  ggb(BR,
      deaths.summed = FALSE,
      minA=5,maxA=75,
      lm.method = .combos[1],
      opt.method = .combos[2])
}

BM <-
 combos %>% 
  group_by(lm.method, opt.method) %>% 
  do(testBR(.combos = .data,BR=BrasilMales)) %>% 
  ungroup() %>% 
  select(id,Mxcoverage,lower,upper,lm.method,opt.method,r2)

BM %>% 
  group_by(id) %>% 
  summarize(min=min(Mxcoverage),
            q25=quantile(Mxcoverage,.25),
            q50=quantile(Mxcoverage,.5),
            q75=quantile(Mxcoverage,.75),
            max=max(Mxcoverage)) %>% 
  filter(min < 0)

BM %>% 
  filter(id==22)

BrasilMales %>% 
  filter(id == 22) %>% 
  ggplot(aes(x=age,y=pop2/pop1))+
  geom_line()


ggbSensitivityPlot <-
  function(X, 
           minA =5, 
           MaxA = 75,
           minAges = 8,
           deaths.summed = FALSE,
           mig.summed = deaths.summed,
           lm.method = "tukey", 
           nx.method = 2,
           type = "Mxcoverage"){
    
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
    mat   <- mat %>% 
             as_tibble() %>% 
             mutate(value = NA,
                    optimum = FALSE)
    
    for (i in 1:length(res)){
      #scl<- sqrt(diff(range(agesL[[i]])))
      fc          <- as.character(min(agesL[[i]]))
      tc          <- as.character(max(agesL[[i]]))
      
      mn <- ggb(ZA,
                exact.ages = agesL[[i]],
                deaths.summed =  deaths.summed,
                mig.summed = mig.summed,
                lm.method = lm.method,
                nx.method = nx.method)
      
      if (type == "Mxcoverage"){
        MAT[fc,tc]  <- mn$Mxcoverage
      }
      if (type == "resid"){
        MAT[fc,tc]  <-
          ggbgetRMS(agesL[[i]],
                    codi,
                    lm.method = lm.method,
                    opt.method = opt.method,
                    scale = 1)
      }
    }
    
     
  }

