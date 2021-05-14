#################################
###Pleiotropy of Cis and Trans-Regulatory Mutations Analysis and Figure Generation
################################
#Contents:
#1. Set up --> line 10
#2. Comparison of fitness effects of cis and trans regulatory mutations --> line 49
#3. Comparison of genome-wide effects of cis and trans regulatory mutations --> line 146
#4. Comparison of downstream effects of chagning TDH3 expression in cis and in trans --> line 316


##############################
###1. Set up - Files loaded here are generated in "DEseqProcessing.R"
##############################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(data.table)
library(eulerr)
library(UpSetR)
library(gplots)
library(bc3net)
library(ggplot2)
library(grid)
library(egg)
library(car)
library(RColorBrewer)
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
load("~/Documents/Output/Projects/Pleiotropy/DEworkspace/separatenooutlierspostcontrast.RData")
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Combined/4.0" #Directory for figure output
outputdir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Combined/4.0" #Directory for file output

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
DEseqSEs <- c(0.323, 0.286, 0.286, 0.286, 0.373) #Written in from reading out the numbers in the dataset
DESeqvals <- c(-7.22, -2.48, -0.948, -0.200, 0.433)
plotgrowth <- data.frame(DESeqvals, DEseqSEs)
plotgrowth[,"COLLECTION"] <- c("1177","1156","1200","1188","3059") 
plotgrowth <- left_join(plotgrowth, growthstats, by = "COLLECTION")
plotgrowth[6,] <- c(0,0,3016,1,0,0,1,1)

ggplot(data = plotgrowth, aes(y = REL.r.u, x = 2^(DESeqvals))) +
  geom_point(aes(x = 2^(DESeqvals), y = REL.r.u), size = 5) +
  geom_errorbarh(aes(xmin = 2^(DESeqvals - DEseqSEs), xmax = 2^(DESeqvals + DEseqSEs))) +
  geom_errorbar(aes(ymin = REL.r.lower, ymax = REL.r.upper)) +
  labs(y = "Relative Fitness", x = "TDH3 expression \n(Fold change from Wild Type)") +
  geom_smooth(method = "loess", color = "black") +
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

#Which predictions are outside the 95% CIs?
growthstats$GPSIG <- ifelse(growthstats$REL.r.lower < growthstats$Predgrowth & growthstats$Predgrowth < growthstats$REL.r.upper, "YES","NO")

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = 2^(TDH3level), y = REL.r.u, label = nickname)) +
  geom_point(data = growthstats[growthstats$condition %in% transvec,], aes(x = 2^(TDH3level), y = REL.r.u, color = GPSIG), size = 5) +
  scale_color_manual(values = c("#336879", "grey")) +
  geom_errorbar(data = growthstats[growthstats$condition %in% transvec,], aes(ymin = REL.r.lower, ymax = REL.r.upper, color = GPSIG)) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% transvec,], aes(xmin = 2^(TDH3level - TDH3SE), xmax = 2^(TDH3level + TDH3SE), color = GPSIG)) +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = 2^(TDH3level), y = REL.r.u), color = "#f47a60", size = 5) +
  geom_errorbar(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(ymin = REL.r.lower, ymax = REL.r.upper), color = "#f47a60") +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(xmin = 2^(TDH3level - TDH3SE), xmax = 2^(TDH3level + TDH3SE)), color = "#f47a60") +
  geom_smooth(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = 2^(TDH3level), y = REL.r.u), method = "loess", color = "#f47a60") +
  THEMEMAIN() +
  theme(legend.position = "none") +
  xlab("TDH3 Expression\n(Fold Change from Wild Type)") +
  ylab("Relative Fitness") +
  ylim(0.2,1.1) + #For use if you want to plot just the cis strains with the same axis values
  xlim(0,1.5)
#geom_text() #To look at strain identities
ggsave("Cistransfit.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

growthstats$Pleiofit <- growthstats$REL.r.u - growthstats$Predgrowth

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = Pleiofit)) +
  geom_density(color = "black", fill = "black", alpha = 0.7) +
  geom_histogram(color = "black", fill = "#336879") +
  THEMEMAIN() +
  xlab("Pleiotropic Fitness Effect") +
  ylab('Frequency') +
  #geom_vline(xintercept = 0, color = "black") +
  xlim(-0.7, 0.2)
ggsave("distpleiofit.pdf", plot = last_plot(), path = figdir, width = 7, height = 8)

#Testing whether positive effects are significantly smaller than negative effects
wilcox.test(abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit > 0,"Pleiofit"]),abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit < 0,"Pleiofit"]))
median(abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit > 0,"Pleiofit"]), na.rm = TRUE)
median(abs(growthstats[growthstats$condition %in% transvec & growthstats$Pleiofit < 0,"Pleiofit"]), na.rm = TRUE)

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
  }
}

#Overall relationship between DEGs amd growth rate
ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = DEGs, y = REL.r.u)) +
  geom_point(aes(color = type), size = 3) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  geom_smooth(method = "lm", color = "black") +
  THEMEMAIN() +
  xlab("Number of DEGs") +
  ylab("Relative Growth Rate")
ggsave("DEGsvsgrowth.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
growthmod <- lm(growthstats[growthstats$condition %in% c(transvec,cissamplesvec),"REL.r.u"] ~ growthstats[growthstats$condition %in% c(transvec,cissamplesvec),"DEGs"])
summary(growthmod)

growthstats[growthstats$condition %in% c(cissamplesvec, transvec),"type"] <- ifelse(growthstats[growthstats$condition %in% c(cissamplesvec, transvec),"condition"] %in% cissamplesvec, "CIS","TRANS")
growthstats[is.na(growthstats$type),"type"] <- "OTHER"
growthstats[growthstats$condition == "A","DEGs"] <- 0 #Putting in zero DEGs for the reference strain

ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = type, y = log10(DEGs))) +
  geom_boxplot() +
  geom_point(aes(color = type), size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  THEMEMAIN() +
  xlab("Relative Position") +
  ylab("Number of DE genes\n(Log10)") +
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
  geom_vline(xintercept = 64, color = "red") +
  THEMEMAIN() +
  xlab("Median of 5 sampled DEGs") +
  ylab("Frequency")
ggsave("permuteDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevec[permutevec <= 64]) #108 out of 1000

ggplot(data = as.data.frame(permutevecvar), aes(x = permutevecvar)) +
  geom_histogram() +
  geom_vline(xintercept = 10139, color = "red") +
  THEMEMAIN() +
  xlab("Variance of 5 sampled DEGs") +
  ylab("Frequency")
ggsave("permuteDEGsvars.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevecvar[permutevecvar <= 10139]) #12 out of 1000

#Are their variances different?
#Make the log10 of the wild type zero, not -Inf
growthstats[growthstats$condition == "A","DEGs"] <- 1

#Converting to percentages for easier plotting
growthstats$TDH3per <- growthstats$TDH3Ratio*100
growthstats$TDH3upr <- (2^(growthstats$TDH3level + growthstats$TDH3SE))*100
growthstats$TDH3lwr <- (2^(growthstats$TDH3level - growthstats$TDH3SE))*100

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = log10(DEGs))) +
  geom_point(color = "#336879", size = 5) +
  geom_errorbarh(aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#336879") +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = log10(DEGs)), color = "#f47a60", size = 5) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#f47a60") +
  geom_line(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = log10(DEGs)), color = "#f47a60") +
  THEMEMAIN() +
  xlab("TDH3 Expression Level\n(% Wild Type)") +
  ylab("Number of DE Genes\n(Log10)")
  #ylim(0,4) #Again, to use if you want to plot just cis alone
ggsave("TDH3vsDEGs.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

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
FCtogether <- cbind(Translogfcs, Cislogfcs)
FCtogether <- FCtogether[rownames(FCtogether) != "YGR192C",] #Removing TDH3 itself

for (i in 1:ncol(FCtogether)) {
  growthstats[which(growthstats$condition == colnames(FCtogether)[i]),"Eudist"] <- sqrt(sum((FCtogether[,i])^2))
}
growthstats[growthstats$condition == "A", "Eudist"] <- 0 #Placing the control at the origin

ggplot(data = growthstats[growthstats$condition %in% c(transvec,cissamplesvec),], aes(x = type, y = Eudist)) +
  geom_boxplot() +
  geom_point(aes(color = type), size = 3, position = position_jitter(width = 0.3)) +
  scale_color_manual(values = c("#f47a60","#336879")) +
  THEMEMAIN() +
  xlab("Relative Position") +
  ylab("Euclidean Distance\n(To Wild Type)") +
  theme(legend.position = "none")
ggsave("cisvstranseudist.pdf", plot = last_plot(), path = figdir, width = 4, height = 8)

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
  geom_vline(xintercept = 25.8, color = "red") +
  THEMEMAIN() +
  xlab("Median of 5 sampled \nEuclidean distances") +
  ylab("Frequency")
ggsave("permuteEudists.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevec[permutevec <= 25.8]) #254 out of 1000

#Are their variances different?
ggplot(data = as.data.frame(permutevecvar), aes(x = permutevecvar)) +
  geom_histogram() +
  geom_vline(xintercept = 26, color = "red") +
  THEMEMAIN() +
  xlab("Variance of 5 sampled\n Euclidean distances") +
  ylab("Frequency")
ggsave("permuteEudistvars.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
#Calculated p-value from permutations
length(permutevecvar[permutevecvar <= 26]) #21 out of 1000

ggplot(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = Eudist, label = nickname)) +
  geom_point(data = growthstats[growthstats$condition %in% transvec,], aes(x = TDH3per, y = Eudist), color = "#336879", size = 5) +
  geom_errorbarh(aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#336879") +
  geom_point(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(x = TDH3per, y = Eudist), color = "#f47a60", size = 5) +
  geom_errorbarh(data = growthstats[growthstats$condition %in% c(cissamplesvec, "A"),], aes(xmin = TDH3lwr, xmax = TDH3upr), color = "#f47a60") +
  geom_line(data = growthstats[growthstats$condition %in% c(cissamplesvec,"A"),], aes(x = TDH3per, y = Eudist), color = "#f47a60") +
  THEMEMAIN() +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("Euclidean Distance\n(To Wild Type)")
  #geom_text() #To see which strain is which
ggsave("TransvsCiseuclids.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

###########################
## Analysis of Downstream consequences of changing TDH3 expression in cis and in trans
##########################

##########################
## Downstream effects of perturbing TDH3 in cis
Cislogmelt <- melt(Cislogfcs)

#Looking at things misexpressed deletion of TDH3
DelDEgenes <- C[C$padj < 0.1 & !is.na(C$padj),] #154 genes
write.table(DelDEgenes, paste0(outputdir,"/DelDEgenes.txt"), sep = "\t")
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
  ggtitle("Significantly DE in TDH3 Null")
ggsave("ViolinDEdel.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

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
nrow(lmstats[lmstats$padj <= 0.1,]) #132 of 154 
table(lmstats$Coef2 > 0) #55 positive, 99 negative total
table(lmstats[lmstats$padj <= 0.1,"Coef2"] > 0) #of sig, 49 positive, 83 negative

CisonlyfcDelmelt$Siglm <- ifelse(CisonlyfcDelmelt$ID %in% lmstats[lmstats$padj <= 0.1,"Gene"], "Yes","No")
CisonlyfcDelmelt$variable <- as.numeric(as.character(CisonlyfcDelmelt$variable))
ggplot(data = CisonlyfcDelmelt, aes(x=variable,y=value,group=ID,colour=Siglm)) + 
  geom_point()+
  geom_line(alpha = 0.7) +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab(expression(Log[2]~Fold~Change)) +
  ggtitle("Significantly DE in TDH3 Null") +
  THEMEMAIN() +
  scale_color_manual(values = c("grey","#f47a60")) +
  theme(legend.position = 'none') +
  ylim(-2.5,3)
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
lmstatsall <- data.frame(Gene = inallsets)
Percentagesall <- 2^(Cislogfcs)*100
tpercentagesall <- t(Percentagesall)
tpercentagesall <- as.data.frame(tpercentagesall)
tpercentagesall["A",] <- c(rep(100, ncol(tpercentagesall)))
tpercentagesall <- as.data.frame(tpercentagesall)
for (i in 1:ncol(tpercentagesall)) {
  lmmod <-  lm(tpercentagesall[,i] ~ tpercentagesall$YGR192C)
  lmstatsall[i,"pvalue"] <- summary(lmmod)$coefficients[2,4]
  lmstatsall[i,"Coef2"] <- summary(lmmod)$coefficients[2,1]
}

ggplot(data = lmstats, aes(x = Coef2)) +
  geom_histogram(data = lmstatsall, aes(x = Coef2), fill = "grey", alpha = 0.5, bins = 100) +
  geom_histogram(fill = "black", bins = 100) +
  THEMEMAIN() +
  xlab("Regression coefficient") +
  ylab("Density") +
  xlim(-5,5)
ggsave("DelDEregres.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = lmstats, aes(x = Coef2)) +
  geom_histogram(fill = "black", bins = 100) +
  THEMEMAIN() +
  xlab("Regression coefficient") +
  ylab("Density")
ggsave("DelDEregreszoom.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = lmstats, aes(x = pvalue)) +
  geom_histogram(data = lmstatsall, aes(x = pvalue), fill = "grey", alpha = 0.5, bins = 100) +
  geom_histogram(fill = "black", bins = 100) +
  THEMEMAIN() +
  xlab("Regression pvalue") +
  ylab("Density")
ggsave("DelDEpval.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

ggplot(data = lmstats, aes(x = pvalue)) +
  geom_histogram(fill = "black", bins = 100) +
  THEMEMAIN() +
  xlab("Regression pvalue") +
  ylab("Density")
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


#Example of GPD2
ggplot(data = tpercentages, aes(x = YGR192C, y = YOL059W)) + 
  geom_point(color = "#f47a60", size = 3) +
  geom_smooth(color = "#f47a60", method = "lm") +
  geom_point(data = Transpercent, aes(x = YGR192C, y = YOL059W), color = "#336879", size = 3) +
  geom_pointrange(aes(x = YGR192C, ymin = tpercentlower$YOL059W, ymax = tpercentupper$YOL059W), color = "#f47a60") +
  geom_pointrange(aes(y = YOL059W, xmin = tpercentlower$YGR192C, xmax = tpercentupper$YGR192C), color = "#f47a60") +
  geom_pointrange(data = Transpercent,aes(x = YGR192C, ymin = ttransplower$YOL059W, ymax = ttranspupper$YOL059W), color = "#336879") +
  geom_pointrange(data = Transpercent,aes(y = YOL059W, xmin = ttransplower$YGR192C, xmax = ttranspupper$YGR192C), color = "#336879") +
  THEMEMAIN() +
  xlab("TDH3 Expression\n(% Wild Type)") +
  ylab("GPD2 Expression\n(% Wild Type)")
ggsave("TDH3vsGPD2trans.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Making all of these
for (i in 1:length(rownames(Transpercent))) {
  rownames(Transpercent)[i] <- Nameskey[Nameskey$CONDITION == rownames(Transpercent)[i],"NICKNAME"]
}
for (i in 1:nrow(Transpercent)) {
  if (rownames(Transpercent)[i] %in% c("RAP1238","RAP154","RAP1484","RAP1357","GCR1162","GCR1037","GCR1241","GCR1281","GCR1339")) {
    Transpercent[i,"Group"] <- 1
  }
  if (rownames(Transpercent)[i] %in% c("ADE2","ADE4","ADE5","ADE6")) {
    Transpercent[i,"Group"] <- 2
  }
  if (rownames(Transpercent)[i] %in% c("CIA2","FRA1","FTR1","NAR1")) {
    Transpercent[i,"Group"] <- 3
  }
  if (rownames(Transpercent)[i] %in% c("MRN1","BRE2","CYC8","TRA1","TUP1","CAF40","SSN2", "TYE7")) {
    Transpercent[i,"Group"] <- 4
  }
  if (rownames(Transpercent)[i] %in% c("RIM8","HXK2","ATP23","SDS23","NAM7","CCC2","IRA2","WWM1","MOD5","PRE7")) {
    Transpercent[i,"Group"] <- 5
  }
}
Transpercent$Group <- factor(Transpercent$Group)

for (i in 1:ncol(tpercentages)) {
  ggplot(data = tpercentages, aes(x = YGR192C, y = tpercentages[,i])) + 
    geom_point(color = "black", size = 3, shape = 4) +
    geom_smooth(color = "black", method = "lm") +
    geom_point(data = Transpercent, aes(x = YGR192C, y = Transpercent[,colnames(tpercentages)[i]], color = Group), size = 3) +
    theme_bw() +
    xlab("TDH3 Expression\n(% Wild Type)") +
    ylab(paste0(i," Expression\n(% Wild Type)"))
  ggsave(paste0(i,"scatterwcis.pdf"), plot = last_plot(), path = figdir, width = 8, height = 8)
}
dev.off()


#Calculating predictions for all trans mutants for all strains
lmlistsig <- lmlist[names(lmlist) %in% CisonlyfcDelmelt[CisonlyfcDelmelt$Siglm == "Yes","ID"]]
transpredlist <- list()
for (i in 1:length(lmlistsig)) {
  transpredlist[[i]] <- predict.lm(lmlistsig[[i]], newdata = data.frame(YGR192C = Transpercent$YGR192C), interval = "confidence")
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
transinout <- data.frame(t(transresidlist[[1]]))
for (i in 1:length(transresidlist)) {
  transresiddf[i,] <- transresidlist[[i]]
  transinout[i,] <- ifelse(Transpercentdel[,i] <= transpredlist[[i]][,"upr"] & Transpercentdel[,i] >= transpredlist[[i]][,"lwr"], "IN", "OUT")
}

colnames(transresiddf) <- transvec
for (i in 1:length(colnames(transresiddf))) {
  colnames(transresiddf)[i] <- Nameskey[Nameskey$CONDITION == colnames(transresiddf)[i],"NICKNAME"]
}

colnames(transinout) <- transvec
for (i in 1:length(colnames(transinout))) {
  colnames(transinout)[i] <- Nameskey[Nameskey$CONDITION == colnames(transinout)[i],"NICKNAME"]
}

transresiddf$Gene <- names(transpredlist)
transinout$Gene <- names(transpredlist)
transresidmelt <- melt(transresiddf)
transinoutmelt <- melt(transinout, id.vars = "Gene")
transresidmelt$PREDINT <- transinoutmelt$value

#To put them in order of increasing median
for (i in 1:nrow(transresiddf)) {
  transresiddf[i,"Median"] <- median(t(abs(transresiddf[i,1:35])))
}
transresiddf <- transresiddf[order(transresiddf$Median),]
transresidmelt$Gene <- factor(transresidmelt$Gene, levels = transresiddf$Gene)
ggplot(data = transresidmelt, aes(x = Gene, y = abs(value))) +
  geom_point(aes(color = PREDINT)) +
  scale_color_manual(values = c("grey","#336879")) +
  stat_summary(fun = median, geom = "point", size = 4) +
  THEMEMAIN() +
  ylab("Absolute Residual") +
  xlab("Gene (DE in TDH3 Null)") +
  theme(axis.text.x = element_blank(), legend.position = "none") +
  ylim(0,600) #8 large values excluded
ggsave("Transresidsgenes.pdf", plot = last_plot(), path = figdir, width = 15, height = 10)

#Looking at the distribution of residuals mutant-wise

#Manually ordered from largest to smallest absolute residual within each group
transresidmelt$variable <- factor(transresidmelt$variable, levels = c("RAP154","RAP1238","GCR1162","RAP1484","GCR1339","GCR1281","RAP1357","GCR1037","GCR1241",
                                                                      "ADE5","ADE6","ADE4","ADE2",
                                                                      "FTR1","NAR1","CIA2","FRA1",
                                                                      "CAF40","MRN1","BRE2","TRA1","TYE7","SSN2","CYC8","TUP1",
                                                                      "SDS23","IRA2","NAM7","HXK2","PRE7","ATP23","RIM8","WWM1","CCC2","MOD5"))

for (i in 1:nrow(transresidmelt)) {
  if (transresidmelt[i,"variable"] %in% c("RAP1238","RAP154","RAP1484","RAP1357","GCR1162","GCR1037","GCR1241","GCR1281","GCR1339")) {
    transresidmelt[i,"Group"] <- 1
  }
  if (transresidmelt[i,"variable"] %in% c("ADE2","ADE4","ADE5","ADE6")) {
    transresidmelt[i,"Group"] <- 2
  }
  if (transresidmelt[i,"variable"] %in% c("CIA2","FRA1","FTR1","NAR1")) {
    transresidmelt[i,"Group"] <- 3
  }
  if (transresidmelt[i,"variable"] %in% c("MRN1","BRE2","CYC8","TRA1","TUP1","CAF40","SSN2", "TYE7")) {
    transresidmelt[i,"Group"] <- 4
  }
  if (transresidmelt[i,"variable"] %in% c("RIM8","HXK2","ATP23","SDS23","NAM7","CCC2","IRA2","WWM1","MOD5","PRE7")) {
    transresidmelt[i,"Group"] <- 5
  }
}

ggplot(data = transresidmelt, aes(x = variable, y = abs(value))) +
  geom_violin(aes(fill = Group)) +
  stat_summary(fun = median, geom = "point", size = 2, color = "white") +
  #stat_summary(aes(label=round(..y..,3)),fun = median, geom = "text", size = 5.5, vjust = 1.2) + 
  THEMEMAIN() +
  ylab("Absolute Residuals") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  ylim(0,200)
ggsave("Transresidslmabszoom.pdf", plot = last_plot(), path = figdir, width = 15, height = 10)


#Clustering trans mutants according to their residuals to see if they cluster according to category
residslim <- transresiddf[transresiddf$Gene %in% c(lmstats[lmstats$padj <= 0.1,"Gene"]),1:35]
rownames(residslim) <- transresiddf$Gene
trees <- pheatmap(as.matrix(residslim), trace = "none", show_rownames = FALSE)

pdf(paste0(figdir,"/heatmap.pdf"))
trees
dev.off()

ordervec <- as.vector(trees$tree_col$labels)
ordervec <- factor(ordervec, levels = trees$tree_col$labels)
coldf <- data.frame(nickname = ordervec)
for (i in 1:length(trees$tree_col$order)) {
  coldf[i,"nickname"] <- ordervec[trees$tree_col$order[i]]
}
coldf <- left_join(coldf, growthstats[,c("nickname","TDH3level")], by = "nickname")

pdf(paste0(figdir,"/barplot.pdf"), width = 7.5, height = 4)
barplot(coldf$TDH3level)
dev.off()
