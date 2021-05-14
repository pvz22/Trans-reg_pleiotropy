########################################
###Growth rate analysis from Bioscreener Data
#########################################
#Contents:
#1. Reading in data and general processing ---> Line 14
#2. Fitting growth curves to each well --> line 58
#3. Calculating summary statistics for each mutant strain relative to wild type --> line 123

#General set up
rm(list = ls())
library(growthcurver)
library(dplyr)
growthdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Growthrate" #Must be modified to provide the path to the stored raw bioscreener data.

#########################################
##1. Reading in data from Biotek Synergy and initial processing

##Read in the data
Rep1data <- read.table(paste0(growthdir,"/Rep1GR.txt"), header = TRUE, stringsAsFactors = FALSE)
Rep2data <- read.table(paste0(growthdir,"/Rep2GR.txt"), header = TRUE, stringsAsFactors = FALSE)
Rep3data <- read.table(paste0(growthdir,"/Rep3GR.txt"), header = TRUE, stringsAsFactors = FALSE)
Rep4data <- read.table(paste0(growthdir,"/Rep4GR.txt"), header = TRUE, stringsAsFactors = FALSE)
##Something grew in the third replicate's blank designated well, so I 
#read in the raw non-subtracted file, and am going to subtract background using a different blank, well G9.
for (i in 1:nrow(Rep3data)) {
  for (j in 1:ncol(Rep3data)) {
    Rep3data[i,j] <- Rep3data[i,j] - Rep3data[i,"G9"]
  }
}

#Looking at the raw data
matplot(Rep1data[,2:97], type = "l", xlim = c(0,25))
matplot(Rep2data[,2:97], type = "l", xlim = c(0,45))
matplot(Rep3data[,2:97], type = "l", xlim = c(0,45))
matplot(Rep4data[,2:97], type = "l", xlim = c(0,45))

#Practicing with one sample, and testing different timepoints at which to cut off after the diauxic shift
timepassed <- rep(1:18)
timepassed2 <- rep(1:30)
gc_fit <- SummarizeGrowth(timepassed, Rep4data[1:18,"E9"])
gc_fit$vals
plot(gc_fit)
gc_fit2 <- SummarizeGrowth(timepassed2, Rep4data[1:30,"E9"])
gc_fit2$vals
plot(gc_fit2)

##Transforming so that there are more rows than columns, and cutting down to get rid of the slow increase after the diaxuic shift for a better fit
Rep1datat <- t(Rep1data[1:18,2:97])
Rep2datat <- t(Rep2data[1:18,2:97])
Rep3datat <- t(Rep3data[1:18,2:97])
Rep4datat <- t(Rep4data[1:18,2:97])

###Reading in the sample names for each replicate plate, note that to make things all numerical, 978 is actually the blanks here
StrainID1 <- read.table(paste0(growthdir,"/StrainID1.txt"), sep = "\t")
StrainID2 <- read.table(paste0(growthdir,"/StrainID2.txt"), sep = "\t")
StrainID3 <- read.table(paste0(growthdir,"/StrainID3.txt"), sep = "\t")
StrainID4 <- read.table(paste0(growthdir,"/StrainID4.txt"), sep = "\t")

#Fitting the growth curves to each well
Rep1stats <- data.frame()
for (i in 1:nrow(Rep1datat)) {
  well <- rownames(Rep1datat)[i]
  gc_fit <- SummarizeGrowth(timepassed, Rep1datat[i,])
  Rep1stats[i,"K"] <- gc_fit$vals$k
  Rep1stats[i,"r"] <- gc_fit$vals$r
  Rep1stats[i,"N0"] <- gc_fit$vals$n0
  Rep1stats[i,"t_mid"] <- gc_fit$vals$t_mid
  Rep1stats[i,"sigma"] <- gc_fit$vals$sigma
}

Rep2stats <- data.frame()
for (i in 1:nrow(Rep2datat)) {
  well <- rownames(Rep2datat)[i]
  gc_fit <- SummarizeGrowth(timepassed, Rep2datat[i,])
  Rep2stats[i,"K"] <- gc_fit$vals$k
  Rep2stats[i,"r"] <- gc_fit$vals$r
  Rep2stats[i,"N0"] <- gc_fit$vals$n0
  Rep2stats[i,"t_mid"] <- gc_fit$vals$t_mid
  Rep2stats[i,"sigma"] <- gc_fit$vals$sigma
}

Rep3stats <- data.frame()
for (i in 1:nrow(Rep3datat)) {
  well <- rownames(Rep3datat)[i]
  gc_fit <- SummarizeGrowth(timepassed, Rep3datat[i,])
  Rep3stats[i,"K"] <- gc_fit$vals$k
  Rep3stats[i,"r"] <- gc_fit$vals$r
  Rep3stats[i,"N0"] <- gc_fit$vals$n0
  Rep3stats[i,"t_mid"] <- gc_fit$vals$t_mid
  Rep3stats[i,"sigma"] <- gc_fit$vals$sigma
}

Rep4stats <- data.frame()
for (i in 1:nrow(Rep4datat)) {
  well <- rownames(Rep4datat)[i]
  gc_fit <- SummarizeGrowth(timepassed, Rep4datat[i,])
  Rep4stats[i,"K"] <- gc_fit$vals$k
  Rep4stats[i,"r"] <- gc_fit$vals$r
  Rep4stats[i,"N0"] <- gc_fit$vals$n0
  Rep4stats[i,"t_mid"] <- gc_fit$vals$t_mid
  Rep4stats[i,"sigma"] <- gc_fit$vals$sigma
}

#Adding strain designators to each well
Rep1stats[,"COLLECTION"] <- StrainID1
Rep2stats[,"COLLECTION"] <- StrainID2
Rep3stats[,"COLLECTION"] <- StrainID3
Rep4stats[,"COLLECTION"] <- StrainID4

#Separating the mutants from controls and blanks
Rep1mutants <- Rep1stats[!(Rep1stats$COLLECTION %in% c("3016","978")),]
Rep2mutants <- Rep2stats[!(Rep2stats$COLLECTION %in% c("3016","978")),]
Rep3mutants <- Rep3stats[!(Rep3stats$COLLECTION %in% c("3016","978")),]
Rep4mutants <- Rep4stats[!(Rep4stats$COLLECTION %in% c("3016","978")),]

#Checking that all the controls look good to start normalizing to them
boxplot(Rep1mutants[Rep1mutants$COLLECTION %in% c("30161","30162", "30163"),"r"],Rep2mutants[Rep2mutants$COLLECTION %in% c("30161","30162", "30163"),"r"],Rep3mutants[Rep3mutants$COLLECTION %in% c("30161","30162", "30163"),"r"],Rep4mutants[Rep4mutants$COLLECTION %in% c("30161","30162", "30163"),"r"])
Rep4mutants[Rep4mutants$COLLECTION %in% c("30161","30162", "30163"),]
#All look pretty good EXCEPT for an outlier in Rep4, so I will remove that one, in case it is contaminated
Rep4mutants <- Rep4mutants[!(Rep4mutants$COLLECTION == "30162"),]
boxplot(Rep1mutants[,"r"], Rep2mutants[,"r"], Rep3mutants[,"r"], Rep4mutants[,"r"])

##Now, I'm going to use the average values of these control strains to normalize the rest of the growth rates to, making them each relative to the controls on their plate
##I will then measure mean and sd for replicates across plates.
RELATIVE <- function(x) {
  for (i in 1:nrow(x)) {
    if (!(x[i,"COLLECTION"]) %in% c("30161", "30162", "30163")){
      x[i,"REL.r.u"] <- x[i,"r"]/(mean(x[x$COLLECTION %in% c("30161","30162", "30163"),"r"]))
    }
  }
  return(x)
}

Rep1mutants <- RELATIVE(Rep1mutants)
Rep2mutants <- RELATIVE(Rep2mutants)
Rep3mutants <- RELATIVE(Rep3mutants)
Rep4mutants <- RELATIVE(Rep4mutants)

Rep1mutants <- Rep1mutants[!(Rep1mutants$r >= 1.5),] ##Removing these crazy flocculant ones, CYC8 and SSN2
#Calculing means and standard deviations across replicates for each mutant
Combstats <- data.frame()
for (i in 1:nrow(Rep1mutants)) {
  Combstats[i,"COLLECTION"] <- Rep1mutants[i,"COLLECTION"]
  Combstats[i,"REL.r.u"] <- mean(c(Rep1mutants[i,"REL.r.u"],Rep2mutants[(Rep2mutants$COLLECTION == Combstats[i,"COLLECTION"]),"REL.r.u"],Rep3mutants[(Rep3mutants$COLLECTION == Combstats[i,"COLLECTION"]),"REL.r.u"],Rep4mutants[(Rep4mutants$COLLECTION == Combstats[i,"COLLECTION"]),"REL.r.u"]))
  Combstats[i,"REL.r.sd"] <- sd(c(Rep1mutants[i,"REL.r.u"],Rep2mutants[(Rep2mutants$COLLECTION == Combstats[i,"COLLECTION"]),"REL.r.u"],Rep3mutants[(Rep3mutants$COLLECTION == Combstats[i,"COLLECTION"]),"REL.r.u"],Rep4mutants[(Rep4mutants$COLLECTION == Combstats[i,"COLLECTION"]),"REL.r.u"]))
}

##Manually fixing a few for which there were also outliers, and I just re-calculated the mean and sd manually
Combstats[Combstats$COLLECTION == "3289","REL.r.u"] <- 0.4112339
Combstats[Combstats$COLLECTION == "3289", "REL.r.sd"] <- 0.03202532
Combstats[Combstats$COLLECTION == "2919","REL.r.u"] <- 0.6888549
Combstats[Combstats$COLLECTION == "2919", "REL.r.sd"] <- 0.0456173

#Calculating confidence intervals for plotting
Q <- qnorm(1-0.05/2)
for (i in 1:nrow(Combstats)) {
  Combstats[i,"REL.r.se"] <- Combstats[i,"REL.r.sd"]/2
  Combstats[i,"REL.r.upper"] <- Combstats[i,"REL.r.u"]+Q*Combstats[i,"REL.r.se"]
  Combstats[i,"REL.r.lower"] <- Combstats[i,"REL.r.u"]-Q*Combstats[i,"REL.r.se"]
}
write.table(Combstats, paste0(growthdir,"/Combstats.txt"), sep = "\t")
