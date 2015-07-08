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
source("http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R")
# new data for Brasil to run. By regions.
BR1 <- read.table(file.path("Data","data_Brazil_p1.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BR2 <- read.table(file.path("Data","data_Brazil_p2.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BR3 <- read.table(file.path("Data","data_Brazil_p3.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
x <- BR1
BR1ggb.f       <- ggb(BR1[BR1$sex == "f", ])
BR1bh1.f       <- bh1(BR1[BR1$sex == "f", ], sex = "f")
BR1bh2.f       <- bh2(BR1[BR1$sex == "f", ], sex = "f")
BR1ggb.m       <- ggb(BR1[BR1$sex == "m", ])
BR1bh1.m       <- bh1(BR1[BR1$sex == "m", ], sex = "m")
BR1bh2.m       <- bh2(BR1[BR1$sex == "m", ], sex = "m")


