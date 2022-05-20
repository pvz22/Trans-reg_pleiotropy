#################################
###Pleiotropy of Cis and Trans-Regulatory Mutations Analysis and Figure Generation
################################
#Contents:
#1. Set up
#2. Comparison of fitness effects of cis and trans regulatory mutations
#3. Comparison of genome-wide effects of cis and trans regulatory mutations
#4. Comparison of downstream effects of chagning TDH3 expression in cis and in trans


##############################
###1. Set up - Files loaded here are generated in "DEseqProcessing.R"
##############################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(gplots)
library(grid)
library(egg)
library(car)
library(cowplot)
library(pheatmap)


#Setting up the environment and loading in data
rm(list = ls())
#Setting the ggplot theme for all figures
THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

#Loading in processed data from "DEseqProcessing.R"
load("~/Documents/Output/Projects/Pleiotropy/DEworkspace/separatenooutliersmatremovpostcontrast.RData")
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Combined/4.0" #Directory for figure output
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Combined/4.0/Resubmission" #Directory for file output

#Defining condition designators for cis strains, and paralog deletions also sequenced in the same experiment, and the trans-regulatory mutants
cissamplesvec <- c("C","B","GG","M","MM")
paralogsamplesvec <- c("F","FF","K")
transvec <- samplesvec[!samplesvec %in% c(cissamplesvec,"BB","TT","A","F","FF","K")] #Taking out cis-reg mutants, controls, and paralogs

###########################
## Fitness consequences of cis and trans-regulatory mutations
###########################
#Plotting the relationship between TDH3 expression and growth rate, as calculated in the "Growthrate.R" script
growthdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Growthrate"
growthstats <- read.table(paste0(growthdir,"/Combstats.txt"), sep = "\t", header = 1) 
#Adding in the two flocculant strains with NAs for later use.
growthstats[46,] <- c(2504, NA, NA, NA, NA, NA)
growthstats[47,] <- c(3092, NA, NA, NA, NA, NA)
growthstats$COLLECTION <- as.character(growthstats$COLLECTION)
DEseqSEs <- c(0.321, 0.285, 0.285, 0.285, 0.370) #Written in from reading out the numbers in the dataset
DESeqvals <- c(-7.22, -2.48, -0.947, -0.200, 0.434)
plotgrowth <- data.frame(DESeqvals, DEseqSEs)
plotgrowth[,"COLLECTION"] <- c("1177","1156","1200","1188","3059") 
plotgrowth <- left_join(plotgrowth, growthstats, by = "COLLECTION")
plotgrowth[6,] <- c(0,0,3016,1,0,0,1,1)
plotgrowth$Shape <- c(rep("Mut",5),"WT")

ggplot(data = plotgrowth, aes(y = REL.r.u, x = (2^DESeqvals))) +
  geom_point(aes(x = (2^DESeqvals), y = REL.r.u, shape = Shape), size = 5) +
  geom_errorbarh(aes(xmin = 2^(DESeqvals - DEseqSEs), xmax = 2^(DESeqvals + DEseqSEs))) +
  geom_errorbar(aes(ymin = REL.r.u - REL.r.se, ymax = REL.r.u + REL.r.se)) +
  labs(y = "Relative Fitness", x = "TDH3 expression \n(% Wild Type)") +
  scale_x_continuous(labels = scales::percent) +
  geom_smooth(method = "loess", color = "black", span = 0.8) +
  THEMEMAIN() +
  theme(legend.position = "none") +
  ylim(0.7,1.05)
ggsave("TDH3vsGrowth.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)


#Looking at the growth rate differentials relative to expectations from TDH3 expression
#Adding in the reference at a growth rate of 1
growthstats[growthstats$COLLECTION == 30161,] <- c(3016,1,0,0,1,1)
growthstats[growthstats$COLLECTION == 30162,] <- c(3016,1,0,0,1,1)
growthstats[growthstats$COLLECTION == 30163,] <- c(3016,1,0,0,1,1)

#Converting strain numbers to the condition indicator used for RNAseq analysis
dir2 <- "/Users/petravandezande/Documents" #Directory containing the document with the Nameskey
Nameskey <- read.table(paste0(dir2,"/Nameskey.txt"), sep = "\t", header = 1, stringsAsFactors = FALSE)
Nameskey$COLLECTION <- as.character(Nameskey$COLLECTION)
for (i in 1:nrow(growthstats)) {
  growthstats[i,"condition"] <- Nameskey[Nameskey$COLLECTION == growthstats[i,"COLLECTION"], "CONDITION"]
  growthstats[i,"nickname"] <- Nameskey[Nameskey$COLLECTION == growthstats[i,"COLLECTION"], "NICKNAME"]
}
for (i in 1:nrow(growthstats)) { #Pulling values of TDH3 expression from the DEseq objects in the environment
  try(growthstats[i,"TDH3level"] <- get(growthstats[i,"condition"])["YGR192C","log2FoldChange"])
  try(growthstats[i,"TDH3SE"] <- get(growthstats[i,"condition"])["YGR192C","lfcSE"])
}
growthstats[growthstats$condition == "A","TDH3level"] <- 0 #Adding in a fold change of 0 for the control


#Calculating the distribution of pleiotropic fitness effects
growthstats$TDH3Ratio <- 2^(growthstats$TDH3level) #Converting to proportion
cisgrowthmod <- loess(growthstats[growthstats$condition %in% c(cissamplesvec,"A"),"REL.r.u"] ~ growthstats[growthstats$condition %in% c(cissamplesvec,"A"),"TDH3Ratio"])

for (i in 1:nrow(growthstats)) {
  growthstats[i,"Predgrowth"] <- predict(cisgrowthmod, growthstats[i,"TDH3Ratio"])
}
plot(growthstats$TDH3Ratio, growthstats$Predgrowth) #Check that predictions follow expectation

#Which predictions are outside the 99% CIs for the trans growth data?
growthstats$GPSIG <- ifelse(growthstats$REL.r.lower < growthstats$Predgrowth & growthstats$Predgrowth < growthstats$REL.r.upper, "NO","YES")
growthstats$Shape <- ifelse(growthstats$condition == "A","WT","Mut")

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = (2^(TDH3level)), y = REL.r.u, label = nickname)) +
  geom_point(data = growthstats[growthstats$condition %in% transvec,], aes(x = (2^(TDH3level)), y = REL.r.u, color = GPSIG), size = 5, alpha = 0.7) +
  scale_color_manual(values = c("grey","#336879")) +
  geom_errorbar(data = growthstats[growthstats$condition %in% transvec,], aes(ymin = (REL.r.u - REL.r.se), ymax = (REL.r.u + REL.r.se), color = GPSIG)) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% transvec,], aes(xmin = (2^(TDH3level - TDH3SE)), xmax = (2^(TDH3level + TDH3SE)), color = GPSIG)) +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = (2^(TDH3level)), y = REL.r.u, shape = Shape), color = "#f47a60", size = 5, alpha = 0.7) +
  geom_errorbar(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(ymin = REL.r.u - REL.r.se, ymax = REL.r.u + REL.r.se), color = "#f47a60") +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(xmin = (2^(TDH3level - TDH3SE)), xmax = (2^(TDH3level + TDH3SE))), color = "#f47a60") +
  geom_smooth(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = (2^(TDH3level)), y = REL.r.u), method = "loess", color = "#f47a60", span = 0.8) +
  THEMEMAIN() +
  theme(legend.position = "none") +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("Relative Fitness") +
  ylim(.20,1.10) + #For use if you want to plot just the cis strains with the same axis values
  xlim(0,1.50) +
  scale_x_continuous(labels = scales::percent)
#geom_text() #To look at strain identities
ggsave("Cistransfit.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

growthstats$Pleiofit <- growthstats$REL.r.u - growthstats$Predgrowth

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = Pleiofit)) +
  geom_density(color = "black", fill = "black", alpha = 0.7) +
  geom_histogram(color = "black", fill = "#58a1b8") +
  THEMEMAIN() +
  xlab("Pleiotropic Fitness Effect") +
  ylab('Frequency') +
  #geom_vline(xintercept = 0, color = "black") +
  xlim(-0.7, 0.2)
ggsave("distpleiofit.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Plotting the residual, or pleiotropic, growth effects vs TDH3 expression
ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = (2^(TDH3level)), y = Pleiofit, label = nickname)) +
  geom_point(aes(color = GPSIG), size = 5, alpha = 0.7) +
  geom_smooth(method = MASS::rlm, color = "black") +
  scale_color_manual(values = c("grey","#336879")) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% transvec,], aes(xmin = (2^(TDH3level - TDH3SE)), xmax = (2^(TDH3level + TDH3SE)), color = GPSIG)) +
  THEMEMAIN() +
  theme(legend.position = "none") +
  xlab("TDH3 expression\n(% Wild Type)") +
  ylab("Pleiotropic Fitness") +
  scale_x_continuous(labels = scales::percent)
#geom_text() #To look at strain identities
ggsave("Transpleiofitvstdh3.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)


#Testing whether positive effects are significantly smaller than negative effects
wilcox.test(abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit > 0,"Pleiofit"]),abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit < 0,"Pleiofit"])) # pvalue: 0.001304
median(abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit > 0,"Pleiofit"]), na.rm = TRUE) #0.02
median(abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit < 0,"Pleiofit"]), na.rm = TRUE) #0.11

#Now trying this with the deletion data
deldir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Friend" #Directory containing a matrix of "M values" from Kettering et al (2014)
Mvalues <- read.table(paste0(deldir,"/Mvalues.txt"), sep = "\t") #Right now these are log2fcs
DeldatTDH3 <- data.frame(Mvalues[,"TDH3"]) #Pulling out each deletion mutant's effect on TDH3
DeldatTDH3$TDH3Ratio <- 2^(DeldatTDH3$Mvalues....TDH3..) #Making them a ratio as the model is
rownames(DeldatTDH3) <- rownames(Mvalues)

for (i in 1:nrow(DeldatTDH3)) {
  DeldatTDH3[i,"Predgrowth"] <- predict(cisgrowthmod, DeldatTDH3[i,"TDH3Ratio"])
}
plot(DeldatTDH3$TDH3Ratio, DeldatTDH3$Predgrowth) #This looks wierd, but look at the y axis

fitdatdir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel" #Directory containing fitness data from MacLean et al
datacheck <- read.table(paste0(fitdatdir,"/dbxref.txt"), fill = NA, stringsAsFactors = FALSE) #Converting strain names to match between the two datafiles
Zhangdat <- read.table(paste0(fitdatdir,"/supplementary_data_6.txt"), header = TRUE, stringsAsFactors = FALSE)
alltrans <- data.frame(Strain = rownames(DeldatTDH3), stringsAsFactors = FALSE)
missingvec <- c()
for (i in 1:nrow(alltrans)) {
  if (!alltrans[i,"Strain"] %in% c(tolower(datacheck$V7))) {
    missingvec <- c(missingvec, alltrans[i,"Strain"])
  }
  else {
    try(
      alltrans[i,"FitnessJZstrain"] <- Zhangdat[Zhangdat$Strain == (datacheck[tolower(datacheck$V7) == alltrans[i,"Strain"],"V5"][1]),"Fitness_relative_to_HO_in_YPD"]
    )
  }
}
table(alltrans$Strain == rownames(DeldatTDH3))
DeldatTDH3 <- cbind(DeldatTDH3, alltrans)
DeldatTDH3 <- DeldatTDH3[!is.na(DeldatTDH3$FitnessJZstrain),]
DeldatTDH3$fitness <- as.numeric(DeldatTDH3$FitnessJZstrain)
write.table(DeldatTDH3, paste0(outputdir,"/DeldatTDH3.txt"), sep = "\t")
ggplot(data = DeldatTDH3, aes(x = TDH3Ratio, y = fitness)) +
  geom_point() +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = 2^(TDH3level), y = REL.r.u, shape = Shape), color = "#f47a60", size = 5) +
  geom_smooth(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = 2^(TDH3level), y = REL.r.u), method = "loess", color = "#f47a60") +
  THEMEMAIN() +
  theme(legend.position = "none") +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("Relative Fitness") +
  scale_x_continuous(labels = scales::percent)
ggsave("Cistransfitdel.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

DeldatTDH3$Pleiofit <- DeldatTDH3$fitness - DeldatTDH3$Predgrowth

ggplot(data = DeldatTDH3, aes(x = Pleiofit)) +
  geom_density(color = "black", fill = "black", alpha = 0.7) +
  #geom_histogram(color = "black", fill = "#58a1b8", alpha = 0.7) +
  THEMEMAIN() +
  xlab("Pleiotropic Fitness Effect") +
  ylab('Density')
ggsave("distpleiofitdeldens.pdf", plot = last_plot(), path = figdir, width = 7, height = 8)

DeldatTDH3 <- DeldatTDH3[!is.na(DeldatTDH3$Pleiofit),]
ggplot(data = DeldatTDH3, aes(x = TDH3Ratio, y = Pleiofit)) +
  geom_point() +
  geom_smooth(method = MASS::rlm, color = "black") +
  THEMEMAIN() +
  theme(legend.position = "none") +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("Pleiotropic Fitness") +
  scale_x_continuous(labels = scales::percent)
ggsave("transfitdelvstdh3.pdf", plot = last_plot(), path = figdir, width = 6, height = 6)

###########################
## Genome-wide effects of cis and trans regulatory mutations
###########################
#Note --> some strains found to still have the URA3 bearing plasmid
#URA3pos <- c("G","L","Y",'E',"U","FF","H","HH")
# transvecclean <- transvec[!transvec %in% URA3pos] #for use to verify removing these strains does not influence conclusions

########################
## Number of DE genes in cis vs trans of similar effect size
for (i in 1:nrow(growthstats)) {
  if (growthstats[i,"condition"] %in% c(transvec, cissamplesvec)) {
    x <- get(growthstats[i,"condition"])
    growthstats[i,"DEGs"] <- nrow(x[x$padj <= 0.1 & !is.na(x$padj),])
    growthstats[i,"MedianFC"] <- median(abs(x$log2FoldChange))
    growthstats[i,"MeanFC"] <- mean(abs(x$log2FoldChange))
  }
}

#Overall relationship between DEGs amd growth rate
growthstats[growthstats$condition %in% c(cissamplesvec, transvec),"type"] <- ifelse(growthstats[growthstats$condition %in% c(cissamplesvec, transvec),"condition"] %in% cissamplesvec, "CIS","TRANS")
growthstats[is.na(growthstats$type),"type"] <- "OTHER"

ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = DEGs, y = REL.r.u, label = nickname)) +
  geom_point(aes(color = type), size = 3) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  geom_smooth(method = "lm", color = "black") +
  THEMEMAIN() +
  xlab("Number of DEGs") +
  ylab("Relative Fitness") +
  theme(legend.position = 'none')
#geom_text()
ggsave("DEGsvsgrowth.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
growthmod <- lm(growthstats[growthstats$condition %in% c(transvec,cissamplesvec),"REL.r.u"] ~ growthstats[growthstats$condition %in% c(transvec,cissamplesvec),"DEGs"])
summary(growthmod) #p-value: 2.1e-9, Rsquared: 0.6351

growthstats[growthstats$condition == "A","DEGs"] <- 0 #Putting in zero DEGs for the reference strain

ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = type, y = (DEGs))) +
  geom_boxplot() +
  geom_point(aes(color = type), size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  THEMEMAIN() +
  xlab("Relative Position") +
  ylab("Number of DE genes") +
  theme(legend.position = "none")
ggsave("cisvstransDEGs.pdf", plot = last_plot(), path = figdir, width = 4, height = 8)


#permutation test
set.seed(12)
permutevec <- c()
permutevecvar <- c()
n <- 0
while (n < 1000) {
  n <- n + 1
  x <- sample(growthstats[growthstats$condition %in% transvec, "DEGs"], 5)
  permutevec[n] <- median(x)
  permutevecvar[n] <- (sd(x))^2
}
ggplot(data = as.data.frame(permutevec), aes(x = permutevec)) +
  geom_histogram() +
  geom_vline(xintercept = 59, color = "red", size = 3) + #Must be changed for the cis-mutants median
  THEMEMAIN() +
  xlab("Median of 5 sampled DEGs") +
  ylab("Frequency")
ggsave("permuteDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevec[permutevec <= 59]) #141 out of 1000

ggplot(data = as.data.frame(permutevecvar), aes(x = permutevecvar)) +
  geom_histogram() +
  geom_vline(xintercept = 7647, color = "red", size = 3) +
  THEMEMAIN() +
  xlab("Variance of 5 sampled DEGs") +
  ylab("Frequency")
ggsave("permuteDEGsvars.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevecvar[permutevecvar <= 7646]) #10 out of 1000

#Converting to percentages for easier plotting
growthstats$TDH3per <- growthstats$TDH3Ratio*100
growthstats$TDH3upr <- (2^(growthstats$TDH3level + growthstats$TDH3SE))*100
growthstats$TDH3lwr <- (2^(growthstats$TDH3level - growthstats$TDH3SE))*100

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = log10(DEGs))) +
  geom_point(color = "#336879", size = 5) +
  geom_errorbarh(aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#336879") +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = log10(DEGs), shape = Shape), color = "#f47a60", size = 5) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#f47a60") +
  geom_line(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = log10(DEGs)), color = "#f47a60") +
  THEMEMAIN() +
  xlab("TDH3 Expression Level\n(% Wild Type)") +
  ylab("Number of DE Genes\n(Log10)") +
  theme(legend.position = "none")
#ylim(0,4) #Again, to use if you want to plot just cis alone
ggsave("TDH3vsDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

############################
#Average log2fold change instead of number of DEGs
#Putting in zero for the wild type
growthstats[growthstats$condition == "A","MeanFC"] <- 0

ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = type, y = MeanFC)) +
  geom_boxplot() +
  geom_point(aes(color = type), size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  THEMEMAIN() +
  xlab("Relative Position") +
  ylab("Mean |Fold Change|") +
  theme(legend.position = "none")
ggsave("cisvstransmeanFC.pdf", plot = last_plot(), path = figdir, width = 4, height = 8)

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = MeanFC)) +
  geom_point(color = "#336879", size = 5) +
  geom_errorbarh(aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#336879") +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = MeanFC, shape = Shape), color = "#f47a60", size = 5) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#f47a60") +
  geom_line(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = MeanFC), color = "#f47a60") +
  THEMEMAIN() +
  xlab("TDH3 Expression Level\n(% Wild Type)") +
  ylab("Mean |Fold Change|") +
  theme(legend.position = "none")
ggsave("TDH3vsMeanFC.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Permutation tests
permutevec <- c()
permutevecvar <- c()
n <- 0
while (n < 1000) {
  n <- n + 1
  x <- sample(growthstats[growthstats$condition %in% transvec, "MeanFC"], 5)
  permutevec[n] <- median(x)
  permutevecvar[n] <- (sd(x))^2
}
ggplot(data = as.data.frame(permutevec), aes(x = permutevec)) +
  geom_histogram() +
  geom_vline(xintercept = 0.21, color = "red", size = 3) + #Must be changed for the cis-mutants median
  THEMEMAIN() +
  xlab("Median of 5 sampled Mean Fold Changes") +
  ylab("Frequency")
ggsave("permutemeanfc.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevec[permutevec <= 0.21]) #514 out of 1000

ggplot(data = as.data.frame(permutevecvar), aes(x = permutevecvar)) +
  geom_histogram() +
  geom_vline(xintercept = 0.002, color = "red", size = 3) +
  THEMEMAIN() +
  xlab("Variance of 5 sampled Mean Fold Changes") +
  ylab("Frequency")
ggsave("permutemeanfcvars.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevecvar[permutevecvar <= 0.002]) #63 out of 1000

#########################
## Looking at euclidean distances over all log2 fold changes present in each background
inallsets <- rownames(C)[which(rownames(C) %in% rownames(MM))]
Cislogfcs <- cbind(C[inallsets,"log2FoldChange"], B[inallsets,"log2FoldChange"], GG[inallsets,"log2FoldChange"], M[inallsets,"log2FoldChange"], MM[inallsets,"log2FoldChange"])
colnames(Cislogfcs) <- c("C", "B", "GG", "M", "MM")
rownames(Cislogfcs) <- inallsets

Translogfcs <- data.frame(D = D$log2FoldChange)
for (i in 2:length(transvec)) {
  Translogfcs[,i] <- get(transvec[i])$log2FoldChange
  colnames(Translogfcs)[i] <- transvec[i]
}
rownames(Translogfcs) <- rownames(D)
Translogfcs <- Translogfcs[rownames(Translogfcs) %in% rownames(Cislogfcs),]
Cislogfcs <- Cislogfcs[rownames(Cislogfcs) %in% rownames(Translogfcs),]
write.table(Cislogfcs, paste0(outputdir,"/Cislogfcs.txt"), sep = "\t")
write.table(Translogfcs, paste0(outputdir,"/Translogfcs.txt"), sep = "\t")

FCtogether <- cbind(Translogfcs, Cislogfcs)
FCtogether <- FCtogether[rownames(FCtogether) != "YGR192C",] #Removing TDH3 itself


for (i in 1:ncol(FCtogether)) {
  growthstats[which(growthstats$condition == colnames(FCtogether)[i]),"Eudist"] <- sqrt(sum((FCtogether[,i])^2))
}
growthstats[growthstats$condition == "A", "Eudist"] <- 0 #Placing the control at the origin

ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = type, y = Eudist)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(color = type), size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  THEMEMAIN() +
  xlab("Relative Position") +
  ylab("Euclidean Distance") +
  theme(legend.position = "none")
ggsave("cisvstranseudist.pdf", plot = last_plot(), path = figdir, width = 4, height = 8)

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = Eudist, label = nickname)) +
  geom_point(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = Eudist), color = "#336879", size = 5) +
  geom_errorbarh(aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#336879") +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = Eudist), color = "#f47a60", size = 5) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#f47a60") +
  geom_line(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = TDH3per, y = Eudist), color = "#f47a60") +
  THEMEMAIN() +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("Euclidean Distance")
#geom_text() #To see which strain is which
ggsave("TransvsCiseuclids.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Permutation test for this one
permutevec <- c()
permutevecvar <- c()
n <- 0
while (n < 1000) {
  n <- n + 1
  x <- sample(growthstats[growthstats$condition %in% transvec, "Eudist"], 5)
  permutevec[n] <- median(x)
  permutevecvar[n] <- (sd(x))^2
}
ggplot(data = as.data.frame(permutevec), aes(x = permutevec)) +
  geom_histogram() +
  geom_vline(xintercept = 25.3, color = "red", size = 3) +
  THEMEMAIN() +
  xlab("Median of 5 sampled \nEuclidean distances") +
  ylab("Frequency")
ggsave("permuteEudists.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevec[permutevec <= 25.3]) #275 out of 1000

#Are their variances different?
ggplot(data = as.data.frame(permutevecvar), aes(x = permutevecvar)) +
  geom_histogram() +
  geom_vline(xintercept = 24.6, color = "red", size = 3) +
  THEMEMAIN() +
  xlab("Variance of 5 sampled\n Euclidean distances") +
  ylab("Frequency")
ggsave("permuteEudistvars.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevecvar[permutevecvar <= 24.6]) #30 out of 1000


#With the deletion data
transsizes <- read.table(paste0(deldir,"/transsizeslim.txt"), sep = "\t", stringsAsFactors = FALSE)
boxexample <- transsizes[transsizes$Focalgene == "htl1",]
boxexample["htl1",] <- c("htl1",111,"htl1",111,111,111)
boxexample$type <- ifelse(boxexample$Strain == "htl1","CIS","TRANS")
boxexample$Size <- as.numeric(boxexample$Size)
ggplot(data = boxexample, aes(x = type, y = (Size))) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(color = type), size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  THEMEMAIN() +
  xlab("Relative Position") +
  ylab("Number of DE genes") +
  theme(legend.position = "none")
ggsave("cisvstransDEGshtl1.pdf", plot = last_plot(), path = figdir, width = 4, height = 8)

#For all focal genes
ggplot(data = transsizes, aes(x = log10(cissize), y = log10(trans_med))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "#f47a60", size = 2) +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Log10)") +
  ylab("Median trans pleiotropy\n(Log10)") +
  xlim(0,3.1) +
  ylim(0,3.1)
ggsave("Scatterplotcisvstransmedlog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#What percentage of these points fall above the x=y line?
nrow(transsizes[transsizes$cissize < transsizes$Size,]) #6877/7254. 95%

###########################
## Analysis of Downstream consequences of changing TDH3 expression in cis and in trans
##########################

##########################
## Downstream effects of perturbing TDH3 in cis
Cislogmelt <- melt(Cislogfcs)

#Looking at things misexpressed deletion of TDH3
DelDEgenes <- C[C$padj < 0.1 & !is.na(C$padj),] #141 genes at 10%, 69 at 5%, 26 at 1%, 320 at 20%
write.table(DelDEgenes, paste0(outputdir,"/DelDEgenes2.txt"), sep = "\t") 
nrow(DelDEgenes[DelDEgenes$log2FoldChange < 0,]) #50
nrow(DelDEgenes[DelDEgenes$log2FoldChange > 0,]) #91
DEgenes1 <- rownames(C[C$padj < 0.01 & !is.na(C$padj),])
DEgenes05 <- rownames(C[C$padj < 0.05 & !is.na(C$padj),])
DEgenes20 <- rownames(C[C$padj < 0.2 & !is.na(C$padj),])

#Are these also DE in the trans-mutants?
DelDEgenenames <- rownames(DelDEgenes)

for (i in 1:nrow(growthstats)) {
  if(growthstats[i,"condition"] %in% c(transvec, cissamplesvec)) {
    x <- get(growthstats[i,"condition"])
    growthstats[i,"downDENum"] <- table(rownames(x[x$padj <= 0.1 & !is.na(x$padj),]) %in% DelDEgenenames)["TRUE"]
  }
}

growthstats$downDENum <- ifelse((growthstats$condition %in% c(cissamplesvec, transvec)) & is.na(growthstats$downDENum), 0, growthstats$downDENum)

growthstatsord <- growthstats[order(growthstats$TDH3Ratio),]
growthstatsord$nickname <- factor(growthstatsord$nickname, levels = unique(growthstatsord$nickname))

ggplot(data = growthstatsord[growthstatsord$condition %in% transvec,], aes(x = nickname, y = downDENum)) +
  geom_col() +
  geom_point(data = growthstatsord[growthstatsord$condition %in% transvec,], aes(x = nickname, y = TDH3Ratio*100), color = 'seagreen3', size = 3) +
  scale_y_continuous(sec.axis = sec_axis(~./100, name = "TDH3 Expression")) +
  xlab("Trans-regulatory mutant strain\n(Gene bearing mutation)") +
  ylab("Number of downstream genes DE") +
  THEMEMAIN() +
  theme(axis.text.x = element_text(angle = 90, size = 15))
ggsave("downDEtranscol.pdf", plot = last_plot(), path = figdir, width = 10, height = 8)

# ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = downDENum)) +
#   geom_point(size = 5) +
#   THEMEMAIN() +
#   xlab("TDH3 Expression\n(% Wild Type)") +
#   ylab("Number of downstream genes DE")
# ggsave("DownDEvsTDH3.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
# 
# ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = DEGs, y = downDENum)) +
#   geom_point(size = 5) +
#   THEMEMAIN() +
#   xlab("Total DEGs") +
#   ylab("Number of downstream genes DE")
# ggsave("DownDEvsDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
# 
# ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = PleioDEGs, y = downDENum)) +
#   geom_point(size = 5) +
#   THEMEMAIN() +
#   xlab("Non-downstream DEGs") +
#   ylab("Number of downstream genes DE")
# ggsave("DownDEvsNondownDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)


#When these are not considered for the trans-mutants, what does the relationship look like btw DEGs and pleiotropic fitness?
growthstats$PleioDEGs <- growthstats$DEGs - growthstats$downDENum
ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = PleioDEGs, y = Pleiofit, label = nickname)) +
  geom_point(color = "#336879",size = 3) +
  #geom_text() +
  geom_smooth(method = "lm", color = "black") +
  THEMEMAIN() +
  xlab("Genes DE in parallel to TDH3") +
  ylab("Pleiotroic fitness effect")
ggsave("pleioDEGsvspleiofit.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

pleiomod <- lm(Pleiofit ~ PleioDEGs,data = growthstats[growthstats$condition %in% transvec,])
summary(pleiomod)
#How are the number of DEGs downstream related to the total DEGs?
ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = PleioDEGs, y = downDENum, label = nickname)) +
  geom_point(size = 3) +
  #geom_text() +
  THEMEMAIN() +
  xlab("Genes DE in parallel to TDH3") +
  ylab("Number of \ndownstream genes DE")
ggsave("downvsotherDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

cor.test(growthstats[growthstats$condition %in% transvec,"PleioDEGs"], growthstats[growthstats$condition %in% transvec,"downDENum"]) #0.95

CisonlyfcDel <- as.data.frame(Cislogfcs[rownames(Cislogfcs) %in% rownames(DelDEgenes),])
TableS2 <- Translogfcs[rownames(Translogfcs) %in% rownames(CisonlyfcDel),]
TableS2 <- cbind(CisonlyfcDel, TableS2)
for (i in 1:ncol(TableS2)) {
  colnames(TableS2)[i] <- Nameskey[Nameskey$CONDITION == colnames(TableS2)[i],"NICKNAME"]
}
write.table(TableS2, paste0(outputdir,"/TableS2.txt"), sep = "\t")

CisonlyfcDel$ID <- c(rownames(CisonlyfcDel))
colnames(CisonlyfcDel) <- c("0","20","50","85","135", "ID")
CisonlyfcDelmelt <- melt(CisonlyfcDel, id.var="ID")
onewayanovaDel <- aov(abs(value) ~ variable, data = CisonlyfcDelmelt)
TukeyHSD(onewayanovaDel) #All significantly different except 50-20, 90-20, 90-50

CisonlyfcDelmelt$variable <- factor(CisonlyfcDelmelt$variable, levels = c("0","20","50","85","135"))
ggplot(data = CisonlyfcDelmelt, aes(y = abs(value), x = variable)) + 
  geom_violin(fill = "lightgrey") +
  THEMEMAIN() +
  xlab("TDH3 expression\n(% Wild Type)") +
  ylab(expression(Absolute~Log[2]~Fold~Change)) +
  stat_summary(fun = median, geom = "point", size = 3) +
  #stat_summary(aes(label=round(..y..,3)),fun = median, geom = "text", size = 5.5, vjust = 1.2) + #Manually placed in illustrator
  theme(legend.position = "none") +
  ggtitle("140 genes significantly DE in TDH3 Null")
ggsave("ViolinDEdel.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Looking at a heatmap of all strains and genes excluding TDH3 itself, for downstream genes and others
heatdf <- FCtogether
for (i in 1:length(colnames(heatdf))) {
  colnames(heatdf)[i] <- Nameskey[Nameskey$CONDITION == colnames(heatdf)[i],"NICKNAME"]
}

# pdf(paste0(figdir,"/heatall.pdf"))
# pheatmap(heatdf, scale = "row", show_rownames = FALSE)
# dev.off()
# heatall <- pheatmap(heatdf, scale = "row", show_rownames = FALSE)

heatdown <- heatdf[rownames(heatdf) %in% DelDEgenenames,]
#heatdown <- heatdown[,heatall$tree_col$order]

pheatmap(heatdown, scale = 'row', show_rownames = FALSE, border_color = NA, filename = paste0(figdir,"/heatdownclust.png"))
pheatmap(heatdown, scale = 'row', show_rownames = FALSE, border_color = NA, cluster_cols = FALSE, filename = paste0(figdir,"/heatdown.pdf"))

heatother <- heatdf[!rownames(heatdf) %in% DelDEgenenames,]
#heatother <- heatother[,heatall$tree_col$order]
#pheatmap(heatother, scale = 'row', show_rownames = FALSE, cluster_cols = FALSE, border_color = NA, filename = paste0(figdir, "/heatother.pdf"))
pheatmap(heatother, scale = 'row', show_rownames = FALSE, border_color = NA, filename = paste0(figdir, "/heatotherclust.pdf"), width = 8, height = 11)

#Building linear models for each of these genes
Percentages <- 2^(CisonlyfcDel[,c("0","20","50","85","135")])*100
CisonlyfcDelt <- as.data.frame(t(CisonlyfcDel[,c("0","20","50","85","135")]), stringsAsFactors = FALSE)
tpercentages <- as.data.frame(t(Percentages), stringsAsFactors = FALSE)
tpercentages['100',] <- c(rep(100, ncol(tpercentages))) ##Add in a point at 100 for everything
CisonlyfcDelt["100",] <- c(rep(0, ncol(CisonlyfcDelt)))

lmlist <- list()
lmstats <- data.frame(Gene = colnames(tpercentages), stringsAsFactors = FALSE)
for (i in 1:ncol(tpercentages)) {
  lmlist[[i]] <- lm(tpercentages[,i] ~ YGR192C, data = tpercentages)
  names(lmlist)[i] <- colnames(tpercentages)[i]
}
for (i in 1:nrow(lmstats)) {
  if(lmstats[i,"Gene"] == names(lmlist)[i]) {
    lmstats[i,"pvalue"] <- summary(lmlist[[i]])$coefficients[2,4]
    lmstats[i,"Coef2"] <- summary(lmlist[[i]])$coefficients[2,1]
    lmstats[i,"SE"] <- summary(lmlist[[i]])$coefficients[2,2]
  }
}

#Multiple testings correction for the DE genes
lmstats$padj <- p.adjust(lmstats$pvalue, method = "BH")
lmstats[order(lmstats$padj),] #corresponds to p of 0.09
nrow(lmstats[lmstats$padj <= 0.1,]) #123 of 141 , 32 at 5%, 1 at 1%, 139 at 20%
table(lmstats$Coef2 > 0) #50 positive, 91 negative total
table(lmstats[lmstats$padj <= 0.1,"Coef2"] > 0) #of sig, 45 positive, 78 negative

CisonlyfcDelmelt$Siglm <- ifelse(CisonlyfcDelmelt$ID %in% lmstats[lmstats$padj <= 0.1,"Gene"], "Yes","No")
CisonlyfcDelmelt$variable <- as.numeric(as.character(CisonlyfcDelmelt$variable))
ggplot(data = CisonlyfcDelmelt[CisonlyfcDelmelt$ID != "YGR192C",], aes(x=variable,y=value,group=ID,colour=Siglm)) + 
  geom_point()+
  geom_line(alpha = 0.7) +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab(expression(Log[2]~Fold~Change)) +
  ggtitle("140 genes significantly DE in TDH3 Null") +
  THEMEMAIN() +
  scale_color_manual(values = c("grey","#f47a60")) +
  theme(legend.position = 'none')
ggsave("SpaghettiDelzoom.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Creating matrices of standard error to use for error bars in plots
CisSEs <- cbind(C[inallsets,"lfcSE"], B[inallsets,"lfcSE"], GG[inallsets,"lfcSE"], M[inallsets,"lfcSE"], MM[inallsets,"lfcSE"])
rownames(CisSEs) <- inallsets
CisonlyseDel <- as.data.frame(CisSEs[rownames(CisSEs) %in% rownames(CisonlyfcDel),])
CisonlyseDel$ID <- c(rownames(CisonlyseDel))
Cisfclower <- CisonlyfcDel[,1:5] - CisonlyseDel[,1:5]
Cisfcupper <- CisonlyfcDel[,1:5] + CisonlyseDel[,1:5]
Percentlower <- 2^(Cisfclower)*100
Percentupper <- 2^(Cisfcupper)*100
tpercentlower <- as.data.frame(t(Percentlower))
tpercentupper <- as.data.frame(t(Percentupper))
tpercentlower['100',] <- c(rep(100, ncol(tpercentlower)))
tpercentupper['100',] <- c(rep(100, ncol(tpercentupper)))


#First, the example of GPD2
ggplot(data = tpercentages, aes(x = YGR192C, y = YOL059W)) +
  geom_point(size = 5, color = "#f47a60") +
  geom_smooth(method = "lm", color = "#f47a60") +
  geom_pointrange(aes(x = YGR192C, ymin = tpercentlower$YOL059W, ymax = tpercentupper$YOL059W), color = "#f47a60") +
  geom_pointrange(aes(y = YOL059W, xmin = tpercentlower$YGR192C, xmax = tpercentupper$YGR192C), color = "#f47a60") +
  THEMEMAIN() +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("GPD2 Expression\n(% Wild Type)")
ggsave("TDH3vsGPD2.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Looking at how these coefficients and pvalues compare to all genes in the cis mutants, not just those DE in the TDH3 null.
lmstatsall <- data.frame(Gene = rownames(Cislogfcs), stringsAsFactors = FALSE)
Percentagesall <- 2^(Cislogfcs)*100
tpercentagesall <- t(Percentagesall)
tpercentagesall <- as.data.frame(tpercentagesall)
tpercentagesall["A",] <- c(rep(100, ncol(tpercentagesall)))
tpercentagesall <- as.data.frame(tpercentagesall)
for (i in 1:ncol(tpercentagesall)) {
  if(colnames(tpercentagesall)[i] == lmstatsall[i,"Gene"]) {
    lmmod <-  lm(tpercentagesall[,i] ~ YGR192C, data = tpercentagesall)
    lmstatsall[i,"pvalue"] <- summary(lmmod)$coefficients[2,4]
    lmstatsall[i,"Coef2"] <- summary(lmmod)$coefficients[2,1]
  }
}

ggplot(data = lmstats, aes(x = Coef2)) +
  geom_histogram(data = lmstatsall, aes(x = Coef2), fill = "grey", alpha = 0.5, bins = 100) +
  geom_histogram(fill = "black", bins = 100) +
  THEMEMAIN() +
  xlab("Regression coefficient") +
  ylab("Count") +
  xlim(-5,5)
ggsave("DelDEregres.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = lmstats, aes(x = Coef2)) +
  geom_histogram(data = lmstatsall[lmstatsall$Gene %in% DEgenes20,], aes(x = Coef2), fill = "green", binwidth = 0.1) +
  geom_histogram(fill = "black", binwidth = 0.1) +
  geom_histogram(data = lmstatsall[lmstatsall$Gene %in% DEgenes05,], aes(x = Coef2), fill = "blue", binwidth = 0.1) +
  geom_histogram(data = lmstatsall[lmstatsall$Gene %in% DEgenes1,], aes(x = Coef2), fill = "red", binwidth = 0.1) +
  THEMEMAIN() +
  xlab("Regression coefficient") +
  ylab("Count")
ggsave("DelDEregreszoom.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = lmstats, aes(x = pvalue)) +
  geom_histogram(data = lmstatsall, aes(x = pvalue), fill = "grey", alpha = 0.5, bins = 100) +
  geom_histogram(fill = "black", bins = 100) +
  THEMEMAIN() +
  xlab("Regression pvalue") +
  ylab("Count")
ggsave("DelDEpval.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = lmstats, aes(x = pvalue)) +
  geom_histogram(data = lmstatsall[lmstatsall$Gene %in% DEgenes20,], aes(x = pvalue), fill = "green", binwidth = 0.01) +
  geom_histogram(fill = "black", binwidth = 0.01) +
  geom_histogram(data = lmstatsall[lmstatsall$Gene %in% DEgenes05,], aes(x = pvalue), fill = "blue", binwidth = 0.01) +
  geom_histogram(data = lmstatsall[lmstatsall$Gene %in% DEgenes1,], aes(x = pvalue), fill = "red", binwidth = 0.01) +
  THEMEMAIN() +
  xlab("Regression pvalue") +
  ylab("Count")
ggsave("DelDEpvalzoom.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

######################
##These gene's expression in trans-regulatory mutants

Translogfcst <- t(Translogfcs)
Transpercent <- as.data.frame(2^(Translogfcst)*100)
#For error bars
Transses <- data.frame(D = D$lfcSE)
for (i in 2:length(transvec)) {
  Transses[,i] <- get(transvec[i])$lfcSE
  colnames(Transses)[i] <- transvec[i]
}
rownames(Transses) <- rownames(D)
Transses <- Transses[rownames(Transses) %in% rownames(Cislogfcs),]
Transfclower <- Translogfcs - Transses
Transfcupper <- Translogfcs + Transses
Transplower <- 2^(Transfclower)*100
Transpupper <- 2^(Transfcupper)*100
ttransplower <- as.data.frame(t(Transplower))
ttranspupper <- as.data.frame(t(Transpupper))


#Making all of these
# for (i in 1:length(rownames(Transpercent))) {
#   rownames(Transpercent)[i] <- Nameskey[Nameskey$CONDITION == rownames(Transpercent)[i],"NICKNAME"]
# }
# for (i in 1:nrow(Transpercent)) {
#   if (rownames(Transpercent)[i] %in% c("RAP1238","RAP154","RAP1484","RAP1357","GCR1162","GCR1037","GCR1241","GCR1281","GCR1339")) {
#     Transpercent[i,"Group"] <- 1
#   }
#   if (rownames(Transpercent)[i] %in% c("ADE2","ADE4","ADE5","ADE6")) {
#     Transpercent[i,"Group"] <- 2
#   }
#   if (rownames(Transpercent)[i] %in% c("CIA2","FRA1","FTR1","NAR1")) {
#     Transpercent[i,"Group"] <- 3
#   }
#   if (rownames(Transpercent)[i] %in% c("MRN1","BRE2","CYC8","TRA1","TUP1","CAF40","SSN2", "TYE7")) {
#     Transpercent[i,"Group"] <- 4
#   }
#   if (rownames(Transpercent)[i] %in% c("RIM8","HXK2","ATP23","SDS23","NAM7","CCC2","IRA2","WWM1","MOD5","PRE7")) {
#     Transpercent[i,"Group"] <- 5
#   }
# }
# Transpercent$Group <- factor(Transpercent$Group)
# 
# for (i in 1:ncol(tpercentages)) {
#   ggplot(data = tpercentages, aes(x = YGR192C, y = tpercentages[,i])) + 
#     geom_point(color = "black", size = 3, shape = 4) +
#     geom_smooth(color = "black", method = "lm") +
#     geom_point(data = Transpercent, aes(x = YGR192C, y = Transpercent[,colnames(tpercentages)[i]], color = Group), size = 3) +
#     theme_bw() +
#     xlab("TDH3 Expression\n(% Wild Type)") +
#     ylab(paste0(i," Expression\n(% Wild Type)"))
#   ggsave(paste0(i,"scatterwcis.pdf"), plot = last_plot(), path = figdir, width = 8, height = 8)
# }


#Calculating predictions for all trans mutants for all strains
lmlistsig <- lmlist[names(lmlist) %in% CisonlyfcDelmelt[CisonlyfcDelmelt$Siglm == "Yes","ID"]]
transpredlist <- list()
for (i in 1:length(lmlistsig)) {
  transpredlist[[i]] <- predict.lm(lmlistsig[[i]], newdata = data.frame(YGR192C = Transpercent$YGR192C), interval = "prediction")
  names(transpredlist)[i] <- names(lmlistsig)[i]
}


#Calculating the residuals
Transpercentdel <- Transpercent[,c(lmstats[lmstats$padj <= 0.1,"Gene"])] 
transresidlist <- list()
for (i in 1:length(transpredlist)) {
  if(names(transpredlist)[i] == colnames(Transpercentdel)[i]) {
    transresidlist[[i]] <- Transpercentdel[,i] - transpredlist[[i]][,"fit"]
  }
}

transresiddf <- data.frame(t(transresidlist[[1]]))
for (i in 1:length(transresidlist)) { 
  transresiddf[i,] <- transresidlist[[i]]
}
colnames(transresiddf) <- transvec
for (i in 1:length(colnames(transresiddf))) {
  colnames(transresiddf)[i] <- Nameskey[Nameskey$CONDITION == colnames(transresiddf)[i],"NICKNAME"]
}
transresiddf$Gene <- names(transpredlist)

#Using SE intervals overlapping prediction intervals
Transplowersig <- Transplower[rownames(Transplower) %in% c(transresiddf$Gene),]
Transpuppersig <- Transpupper[rownames(Transpupper) %in% c(transresiddf$Gene),]
table(rownames(Transplowersig) == names(transpredlist))

transinout <- matrix(ncol = 35, nrow = 123, data = NA)
transinout <- data.frame(transinout)
for (i in 1:length(transresidlist)) { 
  transinout[i,] <- ifelse(Transplowersig[i,] >= transpredlist[[i]][,"upr"] | Transpuppersig[i,] <= transpredlist[[i]][,"lwr"], "OUT", "IN")
}

colnames(transinout) <- transvec
for (i in 1:length(colnames(transinout))) {
  colnames(transinout)[i] <- Nameskey[Nameskey$CONDITION == colnames(transinout)[i],"NICKNAME"]
}

transinout$Gene <- names(transpredlist)
transresidmelt <- melt(transresiddf)
transinoutmelt <- melt(transinout, id.vars = "Gene")
transresidmelt$PREDINT <- transinoutmelt$value

#Looking at whether there is a relationship between model SE and how many are 'out'
transinout <- left_join(transinout, lmstats, by = "Gene")
for (i in 1:nrow(transinout)) {
  transinout[i,"NUMOUT"] <- table(t(transinout[i,1:35]))["OUT"]
}
ggplot(data = transinout, aes(x = SE, y = as.numeric(NUMOUT))) +
  geom_point()

table(transinout$NUMOUT)

#To put them in order of increasing median
for (i in 1:nrow(transresiddf)) {
  transresiddf[i,"Median"] <- median(t(transresiddf[i,1:35]))
}
transinout <- left_join(transinout, transresiddf[,c('Gene','Median')], by = "Gene")
ggplot(data = transinout, aes(x = SE, y = abs(Median))) +
  geom_point()

transresiddf <- transresiddf[order(transresiddf$Median),]
transresidmelt$Gene <- factor(transresidmelt$Gene, levels = transresiddf$Gene)
ggplot(data = transresidmelt, aes(x = Gene, y = (value))) +
  geom_point(aes(color = PREDINT), size = 2) +
  scale_color_manual(values = c("grey","#336879")) +
  stat_summary(fun = median, geom = "point", size = 2) +
  THEMEMAIN() +
  ylab("Pleiotropic Expression Effect\n(Residual)") +
  xlab("Genes (n = 123)") +
  theme(axis.text.x = element_blank(), legend.position = "none")
ggsave("Transresidsgenes.pdf", plot = last_plot(), path = figdir, width = 12, height = 8)

#Example of GPD2
GPD2out <- c("ADE4","GCR1162","RAP154","RAP1238","ADE6","ADE5")
Nameskey[Nameskey$NICKNAME %in% GPD2out,]
GPD2outcond <- c("S","RR","J","PP","QQ","SS")
ggplot(data = tpercentages, aes(x = YGR192C, y = YOL059W)) + 
  geom_point(color = "#f47a60", size = 3) +
  geom_smooth(color = "#f47a60", method = "lm") +
  geom_point(data = Transpercent, aes(x = YGR192C, y = YOL059W), color = "darkgrey", size = 3) +
  geom_point(data = Transpercent[GPD2outcond,], aes(x = YGR192C, y = YOL059W), color = "#336879", size = 3) +
  geom_pointrange(aes(x = YGR192C, ymin = tpercentlower$YOL059W, ymax = tpercentupper$YOL059W), color = "#f47a60") +
  geom_pointrange(aes(y = YOL059W, xmin = tpercentlower$YGR192C, xmax = tpercentupper$YGR192C), color = "#f47a60") +
  geom_pointrange(data = Transpercent,aes(x = YGR192C, ymin = ttransplower$YOL059W, ymax = ttranspupper$YOL059W), color = "darkgrey") +
  geom_pointrange(data = Transpercent,aes(y = YOL059W, xmin = ttransplower$YGR192C, xmax = ttranspupper$YGR192C), color = "darkgrey") +
  geom_pointrange(data = Transpercent[GPD2outcond,],aes(x = YGR192C, ymin = ttransplower[GPD2outcond,"YOL059W"], ymax = ttranspupper[GPD2outcond,"YOL059W"]), color = "#336879") +
  geom_pointrange(data = Transpercent[GPD2outcond,],aes(y = YOL059W, xmin = ttransplower[GPD2outcond,"YGR192C"], xmax = ttranspupper[GPD2outcond,"YGR192C"]), color = "#336879") +
  geom_line(data = Transpercent,aes(x = YGR192C, y = c(data.frame(transpredlist$YOL059W)$lwr)), color = "red") +
  geom_line(data = Transpercent,aes(x = YGR192C, y = c(data.frame(transpredlist$YOL059W)$upr)), color = "red") +
  THEMEMAIN() +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("GPD2 Expression\n(% Wild Type)")
ggsave("TDH3vsGPD2trans.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)


#What do the residuals look like for genes that are differentially expressed in that trans mutant and those that are not?
for (i in 1:nrow(transresidmelt)) {
  x <- get(growthstats[growthstats$nickname == transresidmelt[i,"variable"],"condition"])
  transresidmelt[i,"DEintrans"] <- ifelse(transresidmelt[i,"Gene"] %in% rownames(x[x$padj <= 0.1 & !is.na(x$padj),]), "DE","notDE")
}
table(transresidmelt$PREDINT, transresidmelt$DEintrans) 
#       DE notDE
# IN   493  3462
# OUT  301    49

ggplot(data = transresidmelt, aes(x = DEintrans, y = value)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(aes(color = PREDINT), size = 3, position = position_jitter(width = 0.3), alpha = 0.7) +
  scale_color_manual(values = c("darkgrey","#336879")) +
  THEMEMAIN() +
  xlab("Significantly DE in trans-mutant") +
  ylab("Residual")
ggsave("DEintransvsresid.pdf", plot = last_plot(), path = figdir, width = 6, height = 8)

#Saving final 'growthstats' file for easy access
write.table(growthstats, paste0(outputdir,"/growthstatsfinal.txt"), sep = "\t")