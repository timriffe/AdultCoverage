if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
	# if I'm on the laptop
	setwd("/home/tim/git/AdultCoverage/AdultCoverage")
} else {
	if (system("hostname",intern=TRUE) == "PC-403478"){
		# on MPIDR PC
		setwd("U://git//AdultCoverage//AdultCoverage")
	} else {
		# in that case I'm on Berkeley system, and other people in the dept can run this too
		setwd(paste0("/data/commons/",system("whoami",intern=TRUE),"/git/AdultCoverage/AdultCoverage"))
	}
}
getwd()
source("raw/GGBages/BH1.R")
source("raw/GGBages/BH2.R")
source("raw/GGBages/GGB.R")
source("R/cdmltw.R")


BR1 <- read.table(file.path("Data","data_Brazil_p1.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BR2 <- read.table(file.path("Data","data_Brazil_p2.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BR3 <- read.table(file.path("Data","data_Brazil_p3.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# functions to group infants. Assumes data come abridged in the standard
# way, 0, 1-4, 5-9, etc.

# this works for one population subset
groupInfants <- function(codi){
	if (all(c(0,1) %in% codi$age)){
		codi[codi$age == 0,c("pop1","death","pop2")] <-
				colSums(codi[codi$age %in% c(0,1),c("pop1","death","pop2")])
		codi <- codi[codi$age != 1, ]
		codi
	}
	codi
}

# for one period, since that's how we have the data organized.
biggroup <- function(X){
	XL        <- split(X, X$cod)
	XLgrouped <- lapply(XL, groupInfants)
	do.call(rbind, XLgrouped)
}

# need to first group age 0 with 1-4...
BR1 <- rbind(biggroup(BR1[BR1$sex == "m", ]),biggroup(BR1[BR1$sex == "f", ]))
BR2 <- rbind(biggroup(BR2[BR2$sex == "m", ]),biggroup(BR2[BR2$sex == "f", ]))
BR3 <- rbind(biggroup(BR3[BR3$sex == "m", ]),biggroup(BR3[BR3$sex == "f", ]))


# BR1
BR1ggb.f       <- ggb(BR1[BR1$sex == "f", ])
BR1bh1.f       <- bh1(BR1[BR1$sex == "f", ], sex = "f")
BR1bh2.f       <- bh2(BR1[BR1$sex == "f", ], sex = "f")
BR1ggb.m       <- ggb(BR1[BR1$sex == "m", ])
BR1bh1.m       <- bh1(BR1[BR1$sex == "m", ], sex = "m")
BR1bh2.m       <- bh2(BR1[BR1$sex == "m", ], sex = "m")
# BR2
BR2ggb.f       <- ggb(BR2[BR2$sex == "f", ])
BR2bh1.f       <- bh1(BR2[BR2$sex == "f", ], sex = "f")
BR2bh2.f       <- bh2(BR2[BR2$sex == "f", ], sex = "f")
BR2ggb.m       <- ggb(BR2[BR2$sex == "m", ])
BR2bh1.m       <- bh1(BR2[BR2$sex == "m", ], sex = "m")
BR2bh2.m       <- bh2(BR2[BR2$sex == "m", ], sex = "m")
# BR3
BR3ggb.f       <- ggb(BR3[BR3$sex == "f", ])
BR3bh1.f       <- bh1(BR3[BR3$sex == "f", ], sex = "f")
BR3bh2.f       <- bh2(BR3[BR3$sex == "f", ], sex = "f")
BR3ggb.m       <- ggb(BR3[BR3$sex == "m", ])
BR3bh1.m       <- bh1(BR3[BR3$sex == "m", ], sex = "m")
BR3bh2.m       <- bh2(BR3[BR3$sex == "m", ], sex = "m")




Results <- rbind(
		# BR1
		data.frame(Data = "BR1", Sex = "f", method = "BH1", cod = names(BR1bh1.f), result = BR1bh1.f ),
		data.frame(Data = "BR1", Sex = "f", method = "BH2", cod = names(BR1bh2.f), result = BR1bh2.f ),
		data.frame(Data = "BR1", Sex = "f", method = "GGB", cod = names(BR1ggb.f), result = BR1ggb.f ),
		data.frame(Data = "BR1", Sex = "m", method = "BH1", cod = names(BR1bh1.m), result = BR1bh1.m ),
		data.frame(Data = "BR1", Sex = "m", method = "BH2", cod = names(BR1bh2.m), result = BR1bh2.m ),
		data.frame(Data = "BR1", Sex = "m", method = "GGB", cod = names(BR1ggb.m), result = BR1ggb.m ),
		# BR2
		data.frame(Data = "BR2", Sex = "f", method = "BH1", cod = names(BR2bh1.f), result = BR2bh1.f ),
		data.frame(Data = "BR2", Sex = "f", method = "BH2", cod = names(BR2bh2.f), result = BR2bh2.f ),
		data.frame(Data = "BR2", Sex = "f", method = "GGB", cod = names(BR2ggb.f), result = BR2ggb.f ),
		data.frame(Data = "BR2", Sex = "m", method = "BH1", cod = names(BR2bh1.m), result = BR2bh1.m ),
		data.frame(Data = "BR2", Sex = "m", method = "BH2", cod = names(BR2bh2.m), result = BR2bh2.m ),
		data.frame(Data = "BR2", Sex = "m", method = "GGB", cod = names(BR2ggb.m), result = BR2ggb.m ),
		# BR3 
		data.frame(Data = "BR3", Sex = "f", method = "BH1", cod = names(BR3bh1.f), result = BR3bh1.f ),
		data.frame(Data = "BR3", Sex = "f", method = "BH2", cod = names(BR3bh2.f), result = BR3bh2.f ),
		data.frame(Data = "BR3", Sex = "f", method = "GGB", cod = names(BR3ggb.f), result = BR3ggb.f ),
		data.frame(Data = "BR3", Sex = "m", method = "BH1", cod = names(BR3bh1.m), result = BR3bh1.m ),
		data.frame(Data = "BR3", Sex = "m", method = "BH2", cod = names(BR3bh2.m), result = BR3bh2.m ),
		data.frame(Data = "BR3", Sex = "m", method = "GGB", cod = names(BR3ggb.m), result = BR3ggb.m )
)

write.table(Results, sep = ",", row.names = FALSE, file = "Data/ResultsGGBages.csv")	


# results for UF datasets

#showperiod <- function(x,codes){
#	N         <- length(x)
#	
#	
#	plot(1:27,x,type="n")
#	text(1:27,x,names(x))
#}


#plot(Bernardo1[names(ggb.m1),"Males"], (ggb.m1  + bh2.m1)/2*100, type = "n")
#text(Bernardo1[names(ggb.m1),"Males"], (ggb.m1  + bh2.m1)/2*100,names(ggb.m1))
#abline(a=0,b=1)
#
#showperiod(ggb.m1 * 100)
#showperiod(bh1.m1 * 100)
#showperiod(bh2.m1 * 100)
#


## males
#ggb.m1       <- ggb(DM1)
#bh1.m1       <- bh1(DM1, sex = "m")
#bh2.m1       <- bh2(DM1, sex = "m")
#
#ggb.m2       <- ggb(DM2)
#bh1.m2       <- bh1(DM2, sex = "m")
#bh2.m2       <- bh2(DM2, sex = "m")
#
#ggb.m3       <- ggb(DM3)
#bh1.m3       <- bh1(DM3, sex = "m")
#bh2.m3       <- bh2(DM3, sex = "m")

##############################################

# Everton's theory is that BH2 internal ages
# are fixed in his code and flexible in ours,
# and that might explain differences...

##############################################


#
#plot(bh1(DM3, sex = "m", exact.ages = seq(30,65,by=5)) / bh1(DM3, sex = "m"))
#abline(h=1)
#
#codi <- DM3[DM3$cod == 14, ]
#
#codi <- DM3[DM3$cod == 53, ]
#codibh1 <- bh1MakeColumns(codi, minA. = 10, AgeInt. = 5, minAges. = 8, ages. = unique(codi$age), sex. = "m")
#codibh1$Cx
#mean(codibh1$Cx[codibh1$age %in% c(seq(15,65,by=5))])
#
##codi <- DM3[DM3$cod == 14, ]
##head(codi)
##codi[1,c("pop1","death","pop2")] <- colSums(codi[1:2,c("pop1","death","pop2")])
##codi <- codi[-2, ]
#write.table(codi, sep = ",", row.names = FALSE, file = "Data/testData53.csv")	
#
#
# codi <- read.csv("Data/testData.csv")
# codi[1,c("pop1","death","pop2")] <- colSums(codi[1:2,c("pop1","death","pop2")])
# codi <- codi[-2, ]

#range(UFF3bh2.f)
#
#plot(UFF1bh2.f,UFF3bh2.f,asp=1)
#abline(a=0,b=1)
#
#head(UFFggb.f)
#
## UFM
#UFMggb.m       <- ggb(BR2[BR2$sex == "m", ])
#UFMbh1.m       <- bh1(BR2[BR2$sex == "m", ], sex = "m")
#UFMbh2.m       <- bh2(BR2[BR2$sex == "m", ], sex = "m")


