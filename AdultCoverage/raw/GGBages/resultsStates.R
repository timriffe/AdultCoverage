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

# split?

# females:
# first period
DF1 			<- UFF[, c("UF", "AGE", "POP_80", "POP_91")]
colnames(DF1) 	<- c("cod", "age", "pop1", "pop2")
DF1$year1 		<- 1980
DF1$year2 		<- 1991
DF1$sex   		<- "f"
DF1$death 		<- (UFF$DEATH_80 + UFF$DEATH_91) / 2
# second period
DF2 			<- UFF[, c("UF", "AGE", "POP_91", "POP_00")]
colnames(DF2) 	<- c("cod", "age", "pop1", "pop2")
DF2$year1 		<- 1991
DF2$year2 		<- 2000
DF2$sex   		<- "f"
DF2$death 		<- (UFF$DEATH_91 + UFF$DEATH_00) / 2
# third period
DF3 			<- UFF[, c("UF", "AGE", "POP_00", "POP_10")]
colnames(DF3) 	<- c("cod", "age", "pop1", "pop2")
DF3$year1 		<- 2000
DF3$year2 		<- 2010
DF3$sex   		<- "f"
DF3$death 		<- (UFF$DEATH_00 + UFF$DEATH_10) / 2

# males
# first period
DM1 			<- UFM[, c("UF", "AGE", "POP_80", "POP_91")]
colnames(DM1) 	<- c("cod", "age", "pop1", "pop2")
DM1$year1 		<- 1980
DM1$year2 		<- 1991
DM1$sex   		<- "f"
DM1$death 		<- (UFM$DEATH_80 + UFM$DEATH_91) / 2
# second period
DM2 			<- UFM[,c("UF","AGE","POP_91","POP_00")]
colnames(DM2) 	<- c("cod","age","pop1","pop2")
DM2$year1 		<- 1991
DM2$year2 		<- 2000
DM2$sex   		<- "m"
DM2$death 		<- (UFM$DEATH_91 + UFM$DEATH_00) / 2
# third period
DM3 			<- UFM[,c("UF","AGE","POP_00","POP_10")]
colnames(DM3) 	<- c("cod","age","pop1","pop2")
DM3$year1 		<- 2000
DM3$year2 		<- 2010
DM3$sex   		<- "m"
DM3$death 		<- (UFM$DEATH_00 + UFM$DEATH_10) / 2


groupInfants <- function(codi){
	if (all(c(0,1) %in% codi$age)){
		codi[codi$age == 0,c("pop1","death","pop2")] <-
				colSums(codi[codi$age %in% c(0,1),c("pop1","death","pop2")])
		codi <- codi[codi$age != 1, ]
		codi
	}
	codi
}

biggroup <- function(X){
	XL        <- split(X, X$cod)
	XLgrouped <- lapply(XL, groupInfants)
	do.call(rbind, XLgrouped)
}
DM1 <- biggroup(DM1)
DM2 <- biggroup(DM2)
DM3 <- biggroup(DM3)
DF1 <- biggroup(DF1)
DF2 <- biggroup(DF2)
DF3 <- biggroup(DF3)
# UFF

# females
ggb.f1       <- ggb(DF1)
bh1.f1       <- bh1(DF1, sex = "f")
bh2.f1       <- bh2(DF1, sex = "f")

ggb.f2       <- ggb(DF2)
bh1.f2       <- bh1(DF2, sex = "f")
bh2.f2       <- bh2(DF2, sex = "f")

ggb.f3       <- ggb(DF3)
bh1.f3       <- bh1(DF3, sex = "f")
bh2.f3       <- bh2(DF3, sex = "f")

ResultsStates <- rbind(
		# BR1
		data.frame(Period = 1, Sex = "f", method = "BH1", cod = names(bh1.f1), result = bh1.f1 ),
		data.frame(Period = 1, Sex = "f", method = "BHadj", cod = names(bh2.f1), result = bh2.f1 ),
		data.frame(Period = 1, Sex = "f", method = "GGB", cod = names(ggb.f1), result = ggb.f1 ),
		data.frame(Period = 1, Sex = "m", method = "BH1", cod = names(bh1.m1), result = bh1.m1),
		data.frame(Period = 1, Sex = "m", method = "BHadj", cod = names(bh2.m1), result = bh2.m1),
		data.frame(Period = 1, Sex = "m", method = "GGB", cod = names(ggb.m1), result = ggb.m1),
		# BR2
		data.frame(Period = 2, Sex = "f", method = "BH1", cod = names(bh1.f2), result = bh1.f2 ),
		data.frame(Period = 2, Sex = "f", method = "BHadj", cod = names(bh2.f2), result = bh2.f2 ),
		data.frame(Period = 2, Sex = "f", method = "GGB", cod = names(ggb.f2), result = ggb.f2 ),
		data.frame(Period = 2, Sex = "m", method = "BH1", cod = names(bh1.m2), result = bh1.m2),
		data.frame(Period = 2, Sex = "m", method = "BHadj", cod = names(bh2.m2), result = bh2.m2),
		data.frame(Period = 2, Sex = "m", method = "GGB", cod = names(ggb.m2), result = ggb.m2),
		# BR3 
		data.frame(Period = 3, Sex = "f", method = "BH1", cod = names(bh1.f3), result = bh1.f3 ),
		data.frame(Period = 3, Sex = "f", method = "BHadj", cod = names(bh2.f3), result = bh2.f3 ),
		data.frame(Period = 3, Sex = "f", method = "GGB", cod = names(ggb.f3), result = ggb.f3 ),
		data.frame(Period = 3, Sex = "m", method = "BH1", cod = names(bh1.m3), result = bh1.m3),
		data.frame(Period = 3, Sex = "m", method = "BHadj", cod = names(bh2.m3), result = bh2.m3),
		data.frame(Period = 3, Sex = "m", method = "GGB", cod = names(ggb.m3), result = ggb.m3)
)

write.table(ResultsStates, file = "/home/tim/Dropbox/paper Lima Riffe and Queiroz/results/UF_estimates/StatesAutoAgesGGB.csv",
		row.names = FALSE, col.names = colnames(ResultsStates), sep = ",")

################################################################################
# Comparison with Bernardo's results. 
################################################################################


####################################################
# some data prep necessary in order to match codes #
#################################################### 
# read in data, once, for sake of comparison.
# need code table to make decent comparisons
codes <- as.data.frame(matrix(c(11,'Rondônia',		12,'Acre',
						13,'Amazonas',		14,'Roraima',
						15,'Pará',		16,'Amapá',
						17,'Tocantins',		21,'Maranhão',
						22,'Piauí',		23,'Ceará',
						24,'Rio Grande do Norte',		25,'Paraíba',
						26,'Pernambuco',		27,'Alagoas',
						28,'Sergipe',		29,'Bahia',
						31,'Minas Gerais',		32,'Espírito Santo',
						33,'Rio de Janeiro',		35,'São Paulo',
						41,'Paraná',		42,'Santa Catarina',
						43,'Rio Grande do Sul',		50,'Mato Grosso do Sul',
						51,'Mato Grosso',		52,'Goiás',		53,'Distrito Federal'
				),ncol=2,byrow=TRUE), stringsAsFactors = FALSE)
# Bernardo's results for each of the three periods,
# eyeball age selection, results probably averaged
# between the three methods.
Bernardo1 <- structure(list(UF = c("Brasil", "Acre", "Alagoas", "Amapá"   , 
						"Amazonas", "Bahia", "Ceará"   , "Distrito Federal", "Espírito Santo"   , 
						"Goiás"   , "Maranhão"   , "Mato Grosso", "Mato Grosso do Sul", 
						"Minas Gerais", "Pará"   , "Paraíba"   , "Paraná"   , "Pernambuco", 
						"Piauí"   , "Rio de Janeiro", "Rio Grande do Norte", "Rio Grande do Sul", 
						"Rondônia"   , "Roraima", "Santa Catarina", "São Paulo"   , 
						"Sergipe"), Males = c(88.5, 71.4, 83.3, 71.43, 71.43, 86.96, 
						78.12, 96.15, 93.46, 83.33, 64.52, 80.64, 83.33, 90.91, 62.89, 
						80.64, 90.91, 87.72, 68.49, 93.46, 73.53, 90.91, 64.94, 65.36, 
						94.34, 96.15, 90.91), Females = c(84.03, 66.7, 81.3, 69.44, 69.93, 
						81.97, 72.99, 96.15, 92.56, 84.75, 60.68, 78.74, 79.37, 86.86, 
						62.5, 78.12, 85.47, 82.89, 66.67, 90.09, 68.97, 89.29, 60.98, 
						65.79, 92.59, 95.24, 86.96)), .Names = c("UF", "Males", "Females"
		), class = "data.frame", row.names = c(NA, -27L))
Bernardo2 <- structure(list(UF = c("Brasil", "Acre", "Alagoas", "Amapá"   , 
						"Amazonas", "Bahia", "Ceará"   , "Distrito Federal", "Espírito Santo"   , 
						"Goiás"   , "Maranhão"   , "Mato Grosso", "Mato Grosso do Sul", 
						"Minas Gerais", "Pará"   , "Paraíba"   , "Paraná"   , "Pernambuco", 
						"Piauí"   , "Rio de Janeiro", "Rio Grande do Norte", "Rio Grande do Sul", 
						"Rondônia"   , "Roraima", "Santa Catarina", "São Paulo"   , 
						"Sergipe", "Tocantins"), Males = c(87L, 75L, 87L, 74L, 77L, 79L, 
						72L, 90L, 89L, 77L, 48L, 80L, 90L, 81L, 73L, 70L, 93L, 87L, 60L, 
						93L, 76L, 94L, 74L, 80L, 94L, 99L, 89L, 64L), Females = c(79L, 
						79L, 83L, 79L, 78L, 70L, 66L, 91L, 84L, 78L, 42L, 82L, 85L, 80L, 
						76L, 73L, 90L, 88L, 58L, 91L, 71L, 97L, 68L, 75L, 94L, 96L, 92L, 
						67L)), .Names = c("UF", "Males", "Females"), class = "data.frame", row.names = c(NA, 
				-28L))
Bernardo3 <- structure(list(UF = c("Brasil", "Acre", "Alagoas", "Amapá"   , 
						"Amazonas", "Bahia", "Ceará"   , "Distrito Federal", "Espírito Santo"   , 
						"Goiás"   , "Maranhão"   , "Mato Grosso", "Mato Grosso do Sul", 
						"Minas Gerais", "Pará"   , "Paraíba"   , "Paraná"   , "Pernambuco", 
						"Piauí"   , "Rio de Janeiro", "Rio Grande do Norte", "Rio Grande do Sul", 
						"Rondônia"   , "Roraima", "Santa Catarina", "São Paulo"   , 
						"Sergipe", "Tocantins"), Males = c(98L, 91L, 94L, 81L, 92L, 91L, 
						94L, 100L, 99L, 93L, 77L, 94L, 100L, 96L, 78L, 95L, 100L, 98L, 
						95L, 100L, 90L, 100L, 85L, 82L, 100L, 100L, 97L, 86L), Females = c(96L, 
						83L, 88L, 77L, 92L, 86L, 90L, 100L, 96L, 91L, 64L, 92L, 100L, 
						94L, 76L, 90L, 100L, 95L, 90L, 100L, 85L, 100L, 90L, 90L, 100L, 
						100L, 94L, 82L)), .Names = c("UF", "Males", "Females"), class = "data.frame", row.names = c(NA, 
				-28L))

# remove result for Brazil as whole
Bernardo1 			<- Bernardo1[-1,]
Bernardo2 			<- Bernardo2[-1,]
Bernardo3 			<- Bernardo3[-1,]

# create recode vector, named
recvec         		<- codes$V1
names(recvec)  		<- codes$V2

# create code column (easier for Everton to 
# compare codes than alphabetized names...)
Bernardo1$code 		<- recvec[Bernardo1$UF]
Bernardo2$code 		<- recvec[Bernardo2$UF]
Bernardo3$code 		<- recvec[Bernardo3$UF]

# for downstream selection in graphical function
rownames(Bernardo1) <- Bernardo1$code
rownames(Bernardo2) <- Bernardo2$code
rownames(Bernardo3) <- Bernardo3$code

# a once-off function to compare Bernardo's
# artisanal results with the automatic results
compareWithBernardo <- function(ourvec, Bernardo, sex="Males", name = "bh1", period = 1){
	plot(Bernardo[names(ourvec),sex], ourvec*100, type = "n",
			xlab = "Bernardo's results", ylab = name, main = paste(name,sex,"period",period),
			xlim=c(50,120),ylim=c(50,120))
	text(Bernardo[names(ourvec),sex], ourvec*100,names(ourvec))
	abline(a=0,b=1)
}

#############################################
# graphical comparison
# eyeball things... looks good for bh adjusted.
# bh1 is erratic, which is why bh adjusted was
# created in the first place.
#############################################

# period 1 ggb
compareWithBernardo(ggb.m1,Bernardo1,"Males","ggb",1)
compareWithBernardo(ggb.f1,Bernardo1,"Females","ggb",1)

# period 1 bh2
compareWithBernardo(bh2.m1,Bernardo1,"Males","bhadj",1)
compareWithBernardo(bh2.m1,Bernardo1,"Females","bhadj",1)

# period 2 ggb
compareWithBernardo(ggb.m2,Bernardo2,"Males","ggb",2)
compareWithBernardo(ggb.f2,Bernardo2,"Females","ggb",2)

# period 2 bh2
compareWithBernardo(bh2.m2,Bernardo2,"Males","bhadj",2)
compareWithBernardo(bh2.m2,Bernardo2,"Females","bhadj",2)

# period 3 ggb
compareWithBernardo(ggb.m3,Bernardo3,"Males","ggb",3)
compareWithBernardo(ggb.f3,Bernardo3,"Females","ggb",3)

# period 3 bh2
compareWithBernardo(bh2.m3,Bernardo3,"Males","bhadj",3)
compareWithBernardo(bh2.m3,Bernardo3,"Females","bhadj",3)

# 
# test for averaging deaths over 10 years:
DF <- read.table("clipboard",header=FALSE, sep = "\t")
colnames(DF) <- c("age","death")
DF
codi       <- DM3[DM3$cod == 21, ]
codi$death <- DF$death / 10

plot(DF$death / 10, DM3[DM3$cod == 53, "death"])
abline(a=0,b=1)



