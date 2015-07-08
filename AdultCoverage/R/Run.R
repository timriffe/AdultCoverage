
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

ITA <- read.table("Data/italy.test.txt", 
  header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Data in frame : cod, age, pop1, year1, pop2, year2, death (mean of two periods)

# new data for Brasil to run. By regions.
BR1 <- read.table(file.path("Data","data_Brazil_p1.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BR2 <- read.table(file.path("Data","data_Brazil_p2.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)
BR3 <- read.table(file.path("Data","data_Brazil_p3.txt"), 
		header = TRUE, sep = "\t", stringsAsFactors = FALSE)









