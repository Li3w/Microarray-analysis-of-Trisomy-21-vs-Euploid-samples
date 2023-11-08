if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("limma")
BiocManager::install("annotate")
BiocManager::install("hgu133a.db")
BiocManager::install("biomaRt")
getwd()
setwd("E:/Omnics/MyBioinfData")
library(affy)
library(limma)
library(annotate)
library(hgu133a.db)

#Reading CEL files
MyBioinfData <- ReadAffy()
MyBioinfData
colnames(MyBioinfData)
summary(MyBioinfData)

#Normalization with RMA
eset <- rma(MyBioinfData)
write.exprs(eset, file="RMAnormalized.txt")

#Distinguish Down syndrome from normal 
colors = c(rep("lightpink",11), rep("lightblue",14))
par(mfrow=c(2,2))

#Distinguish Raw and normalized plot
boxraw <- boxplot(MyBioinfData, col=colors, cex.axis=".8", las="3",    main="Boxplot RAW")
boxnorm <- boxplot(exprs(eset), col=colors, cex.axis=".8", las="3", main="Boxplot Normalized")
dev.off()
par(mfrow=c(2,2))
histraw <- hist(MyBioinfData, cex.axis=".8", las="3",  main="Histogram RAW")
histnorm <- plotDensity(exprs(eset), cex.axis=".8", las="3", xlab="log intensity", main="Histogram Normalised")
dev.off()

#Identifying differentially expressed genes with Limma
x <- c(rep(1,11), rep(2,14))
design <- model.matrix(~ -1+factor(x))
colnames(design) <- c("group1", "group2")
design

#Add gene annotation for each probeset
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID,"hgu133a.db")
fData(eset) <- data.frame(ID=ID,Symbol=Symbol)
fit <- lmFit(eset, design)	

#Use eBayes to make empirical bayes adjustment
contrast.matrix <- makeContrasts(group1-group2, levels=design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
limma_complete <- topTable(fit2, adjust="fdr", sort.by="B", number=62000)
#View first 10 results
limma_complete[1:10,]

#Export results and plot volcano plot
write.csv(limma_complete, file="limma_complete.csv")
volcanoplot(fit2, highlight=20,names=fit2$genes$Symbol)


