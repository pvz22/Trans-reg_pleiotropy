########################################
###RNAseq analysis pipeline and QC
#########################################
#Contents:
#1. Import of salmon count matrices and DESeq2 modeling and contrasts  ---> line 10
#2. Comparison of TDH3 values expected based on fluorescence and RNA seq estimates  ---> line 224
#3. Other QC  ---> line 317

#########################################
##1. Import of count data and DESeq2 analysis
rm(list = ls())
options(stringsAsFactors = F)
library(tximport)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
library(DESeq2)
library(ggplot2)

THEMEMAIN <- function() { #Theme for figure generation
  theme_bw() +
    theme(axis.text = element_text(size = 25), axis.title = element_text(size = 25), plot.margin = unit(c(2,1,1,1),"cm"), plot.title = element_text(size = 25, hjust = 0.5))
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

figdir <- "/Users/petravandezande/Documents/Figures/Projects/Pleiotropy/Combined/4.0" #Location for figure output
dir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/YFPcounts" #Path must be modified to indicate the directory containing the transgene quantifications
salmondir <- "/Users/petravandezande/Documents/Output/Projects/Pleiotropy/Salmonquants/" #Must be modified to indicate the directory containing the salmon output files and samples file
samples <- read.table(paste0(salmondir,"samples4fix.txt"), header=TRUE) #read in txt file of samples identities
rownames(samples) <- samples$SampleID #give it rownames

#Splitting samples into different background groups
alphasamples <- samples[!samples$condition %in% c('C','B','GG','M','BB','MM','TT','K'),]
auramsamples <- samples[samples$condition %in% c('C','B','GG','M','TT','K'),]
aurapsamples <- samples[samples$condition %in% c('MM','BB'),]
refssamples <- samples[samples$condition %in% c('A',"TT"),]

###################
#Contrasting alpha and a reference strains
##################
#Reading in the data for the transgenes present
files <- file.path(dir, paste0(refssamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- refssamples$SampleID
all(file.exists(files)) #Check that all files are findable
GENEID <- c("NAT","YFP","KAN","URA3")
TXNAME <- c("NAT","YFP","KAN","URA3")
tx2gene <- data.frame(cbind(GENEID, TXNAME))
txi <- tximport(files, type="salmon", tx2gene=tx2gene) #Do the tximport
countsYFPrefs <- txi$counts #For use for later comparisons
#Prepping to read in the counts for all other genes counted by salmon
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID
files <- file.path(salmondir, paste0("Sample",refssamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- refssamples$SampleID
all(file.exists(files)) 
txi1 <- tximport(files, type="salmon", tx2gene = tx2gene) #Do the tximport
countsallrefs <- txi1$counts #For use for later comparisons
#Merging the two matrices so that transgenes and endogenous genes will be analyzed together
txitotalrefs <- list(rbind(txi1$abundance, txi$abundance), rbind(txi1$counts, txi$counts), rbind(txi1$length, txi$length), txi1$countsFromAbundance)
names(txitotalrefs) <- c("abundance","counts","length","countsFromAbundance")
refssamples$condition <- factor(refssamples$condition, levels = c("A","TT"))
ddsTxirefs <- DESeqDataSetFromTximport(txitotalrefs,
                                        colData = refssamples,
                                        design = ~ condition)
keeprefs <- rowMeans(counts(ddsTxirefs)) >= 10 #Removing transcripts with fewer than 10 counts on average
ddsrefs <- ddsTxirefs[keeprefs,]
ddsrefs <- DESeq(ddsrefs)

refres <- results(ddsrefs)
refres <- data.frame(refres)
refres$sig <- ifelse(refres$padj <=0.1 & !is.na(refres$padj), "sig","ns")
# volcano plot
ggplot(data = refres, aes(x=log2FoldChange, y=-log(pvalue))) +
  geom_point(aes(color = sig)) +
  #geom_text(data = refres[refres$sig == "sig",], aes(x=log2FoldChange, y=-log(pvalue), label = rownames(refres[refres$sig == "sig",]))) +
  ylab("-Log(pvalue)") +
  xlab("Log2 Fold Change") +
  THEMEMAIN()
  #xlim(-10,10) +
  #ylim(0,50)
ggsave("Refvolcano.png", plot = last_plot(), path = figdir, width = 8, height = 8)
refsDEGs <- refres[refres$padj <= 0.1 & !is.na(refres$padj),]
write.table(refsDEGs, paste0(outputdir,"/refsDEGs.txt"), sep = "\t")

###################
##Model for the auramsamples
###################
#Weeding out some outliers
auramsamples <- auramsamples[auramsamples$batch != "G",]
#Reading in the data for the transgenes present
files <- file.path(dir, paste0(auramsamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- auramsamples$SampleID
all(file.exists(files)) #Check that all files are findable
GENEID <- c("NAT","YFP","KAN","URA3")
TXNAME <- c("NAT","YFP","KAN","URA3")
tx2gene <- data.frame(cbind(GENEID, TXNAME))
txi <- tximport(files, type="salmon", tx2gene=tx2gene) #Do the tximport
countsYFPauram <- txi$counts #For use for later comparisons
#Prepping to read in the counts for all other genes counted by salmon
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]  # tx ID, then gene ID
files <- file.path(salmondir, paste0("Sample",auramsamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- auramsamples$SampleID
all(file.exists(files)) 
txi1 <- tximport(files, type="salmon", tx2gene = tx2gene) #Do the tximport
countsallauram <- txi1$counts #For use for later comparisons
#Merging the two matrices so that transgenes and endogenous genes will be analyzed together
txitotalauram <- list(rbind(txi1$abundance, txi$abundance), rbind(txi1$counts, txi$counts), rbind(txi1$length, txi$length), txi1$countsFromAbundance)
names(txitotalauram) <- c("abundance","counts","length","countsFromAbundance")
auramsamples$condition <- factor(auramsamples$condition, levels = c("TT","M","GG","B","C","K"))
ddsTxiauram <- DESeqDataSetFromTximport(txitotalauram,
                                    colData = auramsamples,
                                    design = ~ condition)
keepauram <- rowMeans(counts(ddsTxiauram)) >= 10 #Removing transcripts with fewer than 10 counts on average
ddsauram <- ddsTxiauram[keepauram,]
ddsauram <- DESeq(ddsauram)

#Checking out the PCA
vsd <- vst(ddsauram, blind = TRUE)
PCAstat <- plotPCA(vsd)
pdf(paste0(figdir,"/PCAcis.pdf"))
plotPCA(vsd)
dev.off()

##############################
##Modeling the aurapsamples
#######################
#aurapsamples <- aurapsamples[aurapsamples$batch != "G",]
files <- file.path(dir, paste0(aurapsamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- aurapsamples$SampleID
all(file.exists(files)) #Check that all files are findable
tx2gene <- data.frame(cbind(GENEID, TXNAME))
txi <- tximport(files, type="salmon", tx2gene=tx2gene) #Do the tximport
countsYFPaurap <- txi$counts #For use for later comparisons
#Prepping to read in the counts for all other genes counted by salmon
tx2gene <- df[, 2:1]  # tx ID, then gene ID
files <- file.path(salmondir, paste0("Sample",aurapsamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- aurapsamples$SampleID
all(file.exists(files)) 
txi1 <- tximport(files, type="salmon", tx2gene = tx2gene) #Do the tximport
countsallaurap <- txi1$counts #For use for later comparisons
#Merging the two matrices so that transgenes and endogenous genes will be analyzed together
txitotalaurap <- list(rbind(txi1$abundance, txi$abundance), rbind(txi1$counts, txi$counts), rbind(txi1$length, txi$length), txi1$countsFromAbundance)
names(txitotalaurap) <- c("abundance","counts","length","countsFromAbundance")
aurapsamples$condition <- as.factor(aurapsamples$condition)
ddsTxiaurap <- DESeqDataSetFromTximport(txitotalaurap,
                                   colData = aurapsamples,
                                   design = ~ condition)
keepaurap <- rowMeans(counts(ddsTxiaurap)) >= 10 #Removing transcripts with fewer than 10 counts on average
ddsaurap <- ddsTxiaurap[keepaurap,] #Removed 350 genes
ddsaurap <- DESeq(ddsaurap)

#Checking out the PCA
vsd <- vst(ddsaurap, blind = TRUE) #Doing the transformation, because I am doing it for the pca, I am using blind = TRUE
PCAstat <- plotPCA(vsd)
pdf(paste0(figdir,"/PCAover.pdf"))
plotPCA(vsd)
dev.off()

########################################
##Modeling the alpha strains
##################################
#Getting rid of another strange outlier
alphasamples <- alphasamples[alphasamples$SampleID != 120025,]
files <- file.path(dir, paste0(alphasamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- alphasamples$SampleID
all(file.exists(files)) #Check that all files are findable
tx2gene <- data.frame(cbind(GENEID, TXNAME))
txi <- tximport(files, type="salmon", tx2gene=tx2gene) #Do the tximport
countsYFPalpha <- txi$counts #For use for later comparisons
#Prepping to read in the counts for all other genes counted by salmon
tx2gene <- df[, 2:1]  # tx ID, then gene ID
files <- file.path(salmondir, paste0("Sample",alphasamples$SampleID,"_quant"), "quant.sf") #get paths to each file.
names(files) <- alphasamples$SampleID
all(file.exists(files)) 
txi1 <- tximport(files, type="salmon", tx2gene = tx2gene) #Do the tximport
countsallalpha <- txi1$counts #For use for later comparisons
#Merging the two matrices so that transgenes and endogenous genes will be analyzed together
txitotalalpha <- list(rbind(txi1$abundance, txi$abundance), rbind(txi1$counts, txi$counts), rbind(txi1$length, txi$length), txi1$countsFromAbundance)
names(txitotalalpha) <- c("abundance","counts","length","countsFromAbundance")
alphasamples$condition <- as.factor(alphasamples$condition)
ddsTxialpha <- DESeqDataSetFromTximport(txitotalalpha,
                                        colData = alphasamples,
                                        design = ~ condition)
keepalpha <- rowMeans(counts(ddsTxialpha)) >= 10 #Removing transcripts with fewer than 10 counts on average
ddsalpha <- ddsTxialpha[keepalpha,]
ddsalpha <- DESeq(ddsalpha)

#Checking out the PCA
vsd <- vst(ddsalpha, blind = TRUE) #Doing the transformation, because I am doing it for the pca, I am using blind = TRUE
PCAstat <- plotPCA(vsd)
pdf(paste0(figdir,"/PCAtrans.pdf"))
plotPCA(vsd)
dev.off()

#Saving the environment for later use.
save.image("~/Documents/Output/Projects/Pleiotropy/DEworkspace/separatenooutliers.RData")



##Now doing contrasts between each sample and its respective reference
auramIDs <- c("C","B","GG","M","K")
samplesvec <- as.vector(samples$condition)
samplesvec <- unique(samplesvec)
alphaIDs <- samplesvec[!samplesvec %in% c(auramIDs,"BB","A","TT","MM")]

for (i in 1:length(auramIDs)) {
  print(i)
  x <- results(ddsauram, contrast = c("condition", auramIDs[i],"TT"))
  assign(paste0(auramIDs[i]), x)
}
MM <- results(ddsaurap, contrast = c("condition","MM","BB"))

for (i in 1:length(alphaIDs)) {
  print(i)
  x <- results(ddsalpha, contrast = c("condition", alphaIDs[i],"A"))
  assign(alphaIDs[i], x)
}

save.image("~/Documents/Output/Projects/Pleiotropy/DEworkspace/separatenooutlierspostcontrast.RData")

#Reading out processed files for upload to GEO
countstogether <- cbind(countsallalpha, countsallauram, countsallaurap)
write.table(countstogether, paste0(outputdir,"/counts_matrix_all.txt"),sep = "\t")

for (i in 29:length(samplesvec)) {
  readout <- as.data.frame(get(samplesvec[i]))
  write.table(readout, paste0(outputdir,"/Strain_",unique(samples[samples$condition == samplesvec[i],"Description"]),".txt"), sep = "\t")
}
readout <- as.data.frame(get("Q"))
write.table(readout, paste0(outputdir,"/Strain_3287.txt"), sep = "\t")

##########################################
###2. Comparison of TDH3 values expected based on fluorescence and RNA seq estimates
for (i in 1:length(auramIDs)) {
  print(auramIDs[i])
  print(get(auramIDs[i])["YGR192C","log2FoldChange"])
}
MM["YGR192C","log2FoldChange"]

DESeqvals <- c(-7.22, -2.48, -0.948, -0.200, 0.433) #From above
eLifevalues <- c(0,0.15,0.36,0.72,1.75) #Estimates of the same allele's output by fluorescence of a reporter at HO.

#Calculating the relationship between expression of a fusion protein at the native locus and YFP at HO.
Fusionvalues <- scan()
0.988266704
0.976104363
1.084579001
1.077448279
1.093622433
0.971502365
1.030755442
1.015573263
0.813421542
0.656649097
0.494507532
0.293001368
0.877714474
0.841850311
0.915064541
1.000453979
0.705966421
0.250994934
0.323221904
0.492543039

HOvalues <- scan()
0.978677143
0.980687039
1.112930876
1.126887996
1.131533598
0.948929381
1.055380703
1.035307882
0.706999771
0.558010268
0.373306756
0.152258388
0.829276967
0.771882371
0.875659907
1
0.62285138
0.132185592
0.207003534
0.398723194

cor(HOvalues, Fusionvalues, method = "pearson") #Pearson correlation of 0.997
fusionmod <- lm(Fusionvalues ~ HOvalues)
plot(HOvalues, Fusionvalues)
abline(a=0,b=1)
abline(fusionmod)

eLifevaldf <- as.data.frame(eLifevalues)
colnames(eLifevaldf) <- "HOvalues"
predictfus <- predict(fusionmod, eLifevaldf)
plot(predictfus, 2^DESeqvals)
abline(a=0,b=1)
cor(2^DESeqvals, predictfus,  method = c("pearson")) #Pearson correlation of 0.964

##Plotting this nicely with standard errors
for (i in 1:length(auramIDs)) {
  print(auramIDs[i])
  print(get(auramIDs[i])["YGR192C","lfcSE"])
}
MM["YGR192C","lfcSE"]
DEseqSEs <- c(0.323, 0.286, 0.286, 0.286, 0.373) #Writing in from the output above
eLifeSDs <- c(0.001359901, 0.00437716, 0.006645418, 0.011325647, 0.170158216)

plotdat <- data.frame(DESeqvals, predictfus, DEseqSEs, eLifeSDs)

ggplot(data = plotdat, aes(x = predictfus, y = 2^DESeqvals)) +
  geom_point(aes(y = 2^DESeqvals, x = predictfus),  size = 5) +
  geom_errorbar(aes(ymin = 2^(DESeqvals - DEseqSEs), ymax = 2^(DESeqvals + DEseqSEs))) +
  geom_errorbar(aes(xmin = predictfus - eLifeSDs, xmax = predictfus + eLifeSDs)) +
  labs(x = "TDH3 expression \n(Reporter gene estimates)", y = "TDH3 expression \n(RNAseq estimates)") +
  geom_abline() +
  THEMEMAIN() +
  xlim(0,2) +
  ylim(0,2) +
  theme(legend.position = "none")
ggsave("PredictvsDESeq.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)

#################################
##3. Other QC

#Taking a look at MAplots
for (i in 1:length(alphaIDs)) {
  plotMA(get(alphaIDs[i]))
}

##Making sure that p-value distributions look fine and the multiple testings correction is appropriate
b <- ggplot(data = as.data.frame(B), aes(x=pvalue)) +
  geom_histogram(bins = 100) +
  ggtitle("20%") +
  ylim(0,500)

a <- ggplot(data = as.data.frame(C), aes(x=pvalue)) +
  geom_histogram(bins = 100) +
  ggtitle("0%")+
  ylim(0,500)
c <- ggplot(data = as.data.frame(GG), aes(x=pvalue)) +
  geom_histogram(bins = 100) +
  ggtitle("50%")+
  ylim(0,500)
d <- ggplot(data = as.data.frame(M), aes(x=pvalue)) +
  geom_histogram(bins = 100) +
  ggtitle("90%")+
  ylim(0,500)
e <- ggplot(data = as.data.frame(MM), aes(x=pvalue)) +
  geom_histogram(bins = 100) +
  ggtitle("135%")+
  ylim(0,500)
library(cowplot)
plot_grid(a,b,c,d,e)
ggsave("pvaluehist.pdf", plot = last_plot(), path = figdir, width = 8, height = 8)
