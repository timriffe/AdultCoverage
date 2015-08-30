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

UFF <- read.table(file.path("Data","dataUFfemale.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
UFM <- read.table(file.path("Data","dataUFmale.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cnames <- colnames(UFF)
UFF <- UFF[,-1]
UFM <- UFM[,-1]
colnames(UFF) <- colnames(UFM) <- cnames[1:ncol(UFF)]

head(UFF)
head(BR1)

DF1 <- UFF[,c("UF","AGE","POP_80","DEATH_80","POP_91")]
colnames(DF1) <- c("cod","age","pop1","death","pop2")
DF1$year1 <- 1980
DF1$year2 <- 1991
DF1$sex   <- "f"
head(UFF)
DF3 <- UFF[,c("UF","AGE","POP_00","DEATH_00","POP_10")]
colnames(DF3) <- c("cod","age","pop1","death","pop2")
DF3$year1 <- 2000
DF3$year2 <- 2010
DF3$sex   <- "f"

DM3 <- UFM[,c("UF","AGE","POP_00","DEATH_00","POP_10")]
colnames(DM3) <- c("cod","age","pop1","death","pop2")
DM3$year1 <- 2000
DM3$year2 <- 2010
DM3$sex   <- "m"

head(DM3[DM3$cod == 14, ])
# doesn't work well with typical abridged data because
# age groups need to be conformable.
codi <- DM3[DM3$cod == 14, ]
head(DM3)
groupInfants <- function(codi){
	if (all(c(0,1) %in% codi$age)){
		codi[codi$age == 0,c("pop1","death","pop2")] <-
				colSums(codi[codi$age %in% c(0,1),c("pop1","death","pop2")])
		codi <- codi[codi$age != 1, ]
		codi
	}
	codi
}

DM3L <- split(DM3, DM3$cod)

DM3Lgrouped <- lapply(DM3L, groupInfants)
DM3 <- do.call(rbind,DM3Lgrouped)

# UFF
head(DF1)

UFF1ggb.f       <- ggb(DF1)
UFF1bh1.f       <- bh1(DF1, sex = "f")
UFF1bh2.f       <- bh2(DF1, sex = "f")

UFF3ggb.f       <- ggb(DF3)
UFF3bh1.f       <- bh1(DF3, sex = "f")
UFF3bh2.f       <- bh2(DF3, sex = "f")

UFM3ggb.m       <- ggb(DM3)
bh1(DM3, sex = "m")
UFM3bh1.m       <- bh1(DM3, sex = "m", exact.ages = seq(15,65,by=5))
bh1(DM3, sex = "m")
UFM3bh2.m       <- bh2(DM3, sex = "m")

plot(bh1(DM3, sex = "m", exact.ages = seq(30,65,by=5)) / bh1(DM3, sex = "m"))
abline(h=1)

codi <- DM3[DM3$cod == 14, ]

codi <- DM3[DM3$cod == 53, ]
codibh1 <- bh1MakeColumns(codi, minA. = 10, AgeInt. = 5, minAges. = 8, ages. = unique(codi$age), sex. = "m")
codibh1$Cx
mean(codibh1$Cx[codibh1$age %in% c(seq(15,65,by=5))])

#codi <- DM3[DM3$cod == 14, ]
#head(codi)
#codi[1,c("pop1","death","pop2")] <- colSums(codi[1:2,c("pop1","death","pop2")])
#codi <- codi[-2, ]
write.table(codi, sep = ",", row.names = FALSE, file = "Data/testData53.csv")	


 codi <- read.csv("Data/testData.csv")
 codi[1,c("pop1","death","pop2")] <- colSums(codi[1:2,c("pop1","death","pop2")])
 codi <- codi[-2, ]

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


