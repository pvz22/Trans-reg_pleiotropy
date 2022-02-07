##################################
#Processing datafile from Kemmeren et al (2014) and calculating expression pleiotropy for deletion mutants
##################################
#Loading libraries
library(data.table)
library(ggplot2)
#Setting up ggplot theme
figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Combined/4.0" #Directory for figure output
THEMEMAIN <- function() {
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
#Reading in file
data.dir <- "/Users/petravandezande/Documents/Data/Projects/Pleiodel" #Folder storing datafile from Kemmeren et al (2014)
mutantsexwtvar <- fread(paste0(data.dir,"/deleteome_all_mutants_ex_wt_var_controls.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Processing datafile to generate pleiotropy matrix
mutsvswt <- t(mutantsexwtvar[,4:ncol(mutantsexwtvar)])

Mvalues <- mutsvswt[mutsvswt[,1] == "M",]
Mvalues <- Mvalues[,2:ncol(Mvalues)]
Mvalues <- apply(Mvalues, c(1,2), as.numeric)
Mvalues <- as.data.frame(Mvalues, stringsAsFactors = FALSE)
exp.names.vec <- c(rownames(Mvalues))
temp <- strsplit(exp.names.vec, "-del")
temp2 <- lapply(temp, "[[", 1)
mat <- matrix(unlist(temp2), ncol = 1, byrow = TRUE)
mat.df <- as.data.frame(mat, stringsAsFactors = FALSE)
rownames(Mvalues) <- mat.df$V1
colnames(Mvalues) <- mutantsexwtvar$geneSymbol[2:6124]
write.table(Mvalues, paste0(outputdir,"/Mvalues.txt"), sep = "\t")

p.values <- mutsvswt[mutsvswt[,1] == "p_value",]
p.values <- p.values[,2:ncol(p.values)]
p.values <- apply(p.values, c(1,2), as.numeric)
exp.names.vec <- c(rownames(p.values))
temp <- strsplit(exp.names.vec, "-del")
temp2 <- lapply(temp, "[[", 1)
mat <- matrix(unlist(temp2), ncol = 1, byrow = TRUE)
mat.df <- as.data.frame(mat)
rownames(p.values) <- mat.df$V1
p.values <- as.data.frame(p.values, stringsAsFactors = FALSE)
colnames(p.values) <- mutantsexwtvar$geneSymbol[2:6124]
write.table(p.values, paste0(outputdir,"/p.values.txt"), sep = "\t")


Pleiotropy <- matrix(NA, nrow = nrow(Mvalues), ncol = ncol(Mvalues))
Pleiotropy <- ifelse(abs(Mvalues) >= 0.7655347  & p.values <= 0.05, 1, 0) #Mvalues are actually log2ratios, so FC of 1.7 = 0.7655347
#I end up with 3,000 more than in the Kemmeren publication. This may be due to the eqality or that I didn't filter out 'nonresponsive' experiments.

colnames(Pleiotropy) <- colnames(Mvalues)
rownames(Pleiotropy) <- rownames(Mvalues)
colnames(Pleiotropy) <- tolower(colnames(Pleiotropy))
write.table(Pleiotropy, paste0(outputdir,"/Pleiotropy.txt"), sep = "\t")
Pleiotropy <- read.table(paste0(outputdir,"/Pleiotropy.txt"), sep = "\t")

#Removing focal genes that do not show a significant decrease in expression upon deletion
noreducvec <- c()
for (i in 1:nrow(Pleiotropy)) {
  if(rownames(Pleiotropy)[i] %in% colnames(Pleiotropy)) {
    if(Pleiotropy[i,colnames(Pleiotropy) == rownames(Pleiotropy)[i]] == 0) {
      noreducvec <- c(noreducvec,i)
    } 
  }
}
#Removes 127 genes
Pleiotropy <- Pleiotropy[-noreducvec,]

#Calculating cis and trans pleiotropy
Pleiotropy <- Pleiotropy[rownames(Pleiotropy) %in% colnames(Pleiotropy),] #We can only look at focal genes that we can ID trans regulators for
transsizes <- data.frame()
allsizes <- data.frame("Strain" = rownames(Pleiotropy), "Size" = rowSums(Pleiotropy) - 2) #Trans size includes focal gene, and trans regulator itself, so subtract 2
cissize <- rowSums(Pleiotropy) - 1 #cissize includes the focal gene, so subtract 1
for (i in 1:nrow(Pleiotropy)) {
  focalgene <- rownames(Pleiotropy)[i]
  if(sum(Pleiotropy[,focalgene]) > 1) { #To eliminate where it is only the self loop as the trans regulator
    transregs <- rownames(Pleiotropy[Pleiotropy[,focalgene] == 1,])
    trans <- allsizes[rownames(allsizes) %in% transregs,]
    trans$Focalgene <- rep(focalgene, nrow(trans))
    transsizes <- rbind(transsizes, trans)
  }
}
cissize <- data.frame(cissize)
write.table(transsizes, paste0(outputdir,"/transsizeslim.txt"), sep = "\t")
write.table(cissize, paste0(outputdir,"/cissizelim.txt"), sep = "\t")
transsizes$Focalgene <- as.character(transsizes$Focalgene)
transsizes$Strain <- as.character(transsizes$Strain)
#Removing self-loops as trans-regulatory perturbations - 
#also checking that all have self-loops to begin with!!
table(unique(transsizes$Focalgene) %in% unique(transsizes$Strain))
transsizes <- transsizes[which(transsizes$Strain != transsizes$Focalgene),]

#Creating a single value for all trans-regs for each focal gene
trans_avg <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(avg = mean(Size))

trans_med <- transsizes %>%
  group_by(Focalgene) %>%
  summarize(med = median(Size))


#Scatterplot of cis vs trans
for (i in 1:nrow(transsizes)) {
  transsizes[i,"cissize"] <- cissize[rownames(cissize) == transsizes[i,"Focalgene"],"cissize"]
  transsizes[i,"trans_med"] <- trans_med[trans_med$Focalgene == transsizes[i,"Focalgene"],"med"]
  transsizes[i,"trans_avg"] <- trans_avg[trans_avg$Focalgene == transsizes[i,"Focalgene"],"avg"]
}

#Reading this out to save it
write.table(transsizes, paste0(outputdir,"/transsizeslim.txt"), sep = "\t")

ggplot(data = transsizes, aes(x = log10(cissize), y = log10(trans_med))) +
  geom_point(color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "#f47a60", size = 2) +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Log10)") +
  ylab("Median trans pleiotropy\n(Log10)") +
  xlim(0,3.1) +
  ylim(0,3.1)
ggsave("Scatterplotcisvstransmedlog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#All trans-regulators separately rather than the median
ggplot(data = transsizes, aes(x = log10(cissize), y = log10(Size))) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "#f47a60", size = 2) +
  THEMEMAIN() +
  xlab("Cis pleiotropy\n(Log10)") +
  ylab("Trans pleiotropy\n(Log10)")
ggsave("Scatterplotcisvstranslog.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#Looking at effect sizes for all focal genes in cis and in trans
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsize"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Strain"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}
for (i in 1:nrow(transsizes)) {
  transsizes[i,"Effectsizecis"] <- Mvalues[rownames(Mvalues) == transsizes[i,"Focalgene"],tolower(colnames(Mvalues)) == transsizes[i,"Focalgene"]]
}
ggplot(data = transsizes, aes(x = Focalgene, y = Effectsize)) +
  geom_point(color = "#336879", alpha = 0.7) +
  geom_point(aes(x = Focalgene, y = Effectsizecis), color = "#f47a60", alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_blank(), legend.position = "none",axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm")) +
  ylab("Focal Gene Expression\n(Log2 Fold Change)") +
  xlab("Focal Gene")
ggsave("Cisvstranseffectall.pdf", plot = last_plot(), path = figdir, width = 8, height = 7)

