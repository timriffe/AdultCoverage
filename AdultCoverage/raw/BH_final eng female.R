## B-H original method
head(x)
minA. <- 10
bhtable <-  function(x, minA. = 10, AgeInt. = 5, sex = "f"){        ##  Data
	tab     <- data.frame(x)           ##  Dat in data frame : cod, age, pop1, year1, pop2, year2, death
	dif     <- x$year2[1]-x$year1[1]
	cod     <- x$cod
	tab1    <- split(tab,cod)

	cods    <- unique(tab$cod)
	ages    <- sort(unique(tab$age))
# temp
	codi    <- tab1[[1]]

codi$birthdays            <- 0
# iterate over age groups >= 10

for (j in seq_along(ages)[ages >= minA.]) {
  # take geometric average of p1 pop vs p2 pop within same cohort
  codi$birthdays[j]       <- 
    1 / AgeInt. * sqrt(codi$pop1[j - 1] * codi$pop2[j])
} # end age loop

# age-specific growth
 
codi[["growth"]]	          <-  log(codi$pop2 / codi$pop1) / dif
codi$growth[is.infinite(codi$growth)] <- 0

codi$cumgrowth       <-  0
codi$cumgrowth[1]    <-  2.5*codi$growth[1]

for (j in 2:length(ages)){
  codi$cumgrowth[j]  <-  2.5*codi$growth[j]+5*sum(codi$growth[(j-1):1])
}

# stopped here

  codi$death_tab       <- codi$death * exp(codi$cumgrowth)


  ratio           <- sum(codi$death_tab[ages%in%c(10:39)])/sum(codi$death_tab[ages%in%c(40:59)])

  if (sex == "f"){
    # model lifetable
    # based on Bennett & Horiuchi (1984) 
    # "Mortality Estimateion fro Registered Deaths in Less Developed Countries", Demography
    standardratios <- c(1.376,	1.3,	1.233,	1.171,
                        1.115,	1.062,	1.012,	0.964,
                        0.918,	0.872,	0.827,	0.787,
                        0.729,	0.673,	0.617,	0.56,
                        0.501,	0.438,	0.365,	0.298,
                        0.235,	0.175,	0.117)
   #                   ex <- ....
  }
  if (sex == "m"){
    # need to change this:
    standardratios <- c(1.161,	1.094,	1.034,	0.98,	
						0.93,	0.885,	0.842,	0.802,	
						0.763,	0.725,	0.689,	0.648,
						0.609,	0.57,	0.53,	0.49,	
						0.447,	0.401,	0.352,	0.305,	
						0.255,	0.202,	0.147)
                      
    #                   ex <- ....
  }

  AllLevels <- 3:25
  CDlevel   <- splinefun(AllLevels~standardratios)(ratio)
  # open at age 110 ....
  eOpen     <- splinefun(ex[,ncol(ex)]~1:25)(CDlevel)
  
  N              <- nrow(codi)
  codi$pop_a     <-  0
  codi$pop_a     <-  codi$death[N ] * (exp(codi$aberto[N ] * codi$cresc[N ])-((codi$aberto*codi$cresc[N ])^2/6))

for (i in seq(along = tab1)){
  for(j in 18:1){
    codi$pop_a[j-1] <- codi$pop_a[j]*exp(5*codi$growth[j-1])+codi$death[j-1]*exp(2.5*codi$growth[j-1])
  }
}
for (i in seq(along = tab1)){
  codi$Cx  <-  0
}
for (i in seq(along = tab1)){
  for (j in 4:17){
    codi$Cx[j] <-  round(codi$pop_a[j]/codi$birthdays[j],digits=2)
  }
}
for (i in seq(along = tab1)){
  codi$cob1  <-  round((sum(codi$Cx[4:15])/12), digits=2)                     ## para 10 a 65+
  codi$cob2  <-  round((sum(codi$Cx[5:15])/11), digits=2)                     ## para 15 a 65+
  codi$cob3  <-  round((sum(codi$Cx[6:15])/10), digits=2)                     ## para 20 a 65+
  codi$cob4  <-  round((sum(codi$Cx[7:15])/9), digits=2)                      ## para 25 a 65+
  codi$cob5  <-  round((sum(codi$Cx[8:15])/8), digits=2)                      ## para 30 a 65+
}
tab2      <-  unsplit(tab1,cod)
grau1     <-  tapply(tab2$cob1,tab2$cod,mean)
grau2     <-  tapply(tab2$cob2,tab2$cod,mean)
grau3     <-  tapply(tab2$cob3,tab2$cod,mean)
grau4     <-  tapply(tab2$cob4,tab2$cod,mean)
grau5     <-  tapply(tab2$cob5,tab2$cod,mean)
resultado <-  cbind(grau1,grau2,grau3,grau4,grau5)
results   <-  list(resultado=resultado)
return(results)}

bhtable(banco)

