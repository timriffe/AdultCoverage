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
source("raw/BH1.R")
source("raw/BH2.R")
source("raw/GGB.R")
source("R/cdmltw.R")
# need web connection for this hack
# source("http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R")
# doh, website removed. Eliminate parallel dependency for now.
# new data for Brasil to run. By regions.
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
	 
write.table(Results, sep = ",", row.names = FALSE, file = "Data/Results.csv")		
		
		
		
