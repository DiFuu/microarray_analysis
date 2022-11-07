#Important packages
library(GEOquery)
library(affy)
library(affyPLM)
library(genefilter)
library(limma)
library(hgu133plus2.db)


#Auxiliary file 'get_info_kegg.R'.

source('get_info_kegg.R')

#CEL files and phenotypic information

#Work directories
celfiles <- ReadAffy(phenoData = "pheno_data_subset.txt", celfile.path="data")

#When viewing the content of the celfiles object or variable, we see that it indicates that the set has 31 Affymetrix samples or microarrays
#of type hgu133plus2, where 54675 probesets are described. The size of the arrays is 1164*1164. Annotation information is as described by hgu133plus2

celfiles

#Access to the names of the samples or microarrays
sampleNames(celfiles)
#Access to phenotypic information
pData(celfiles)

#Access to the expression values of the probes
eset <- exprs(celfiles)
dim(eset)
#It tells us that we have 1354896 rows, that is, probes

head(eset) #We take a look at the intensity values of some probes

#Raw data quality analysis

#Boxplots
boxplot(celfiles, las=2)

#Histograms
hist(celfiles)

#Pre-processing.
#Objective: to balance the intensity levels of the samples of the whole experiment while maintaining
#the effect due to the treatment under investigation.

#rma performs all 4 stages at once

celfiles.rma <- rma(celfiles)
celfiles.rma

#We access the expression data of the microarrays and check their dimension
eset<-exprs(celfiles.rma)
dim(eset)

head(eset)

#Quality analysis of pre-processed data

boxplot(celfiles.rma,las=2)

hist(celfiles.rma)

#Differential expression analysis
# Objective: to identify those genes whose expression patterns differ significantly according to the phenotype or the experimental conditions.

#In the first place, we are going to carry out a non-specific filtering, prior to the differential expression analysis
#This non-specific filtering consists of eliminating those genes from the chip with very little variability of expression between conditions.
#Objective: (1) to eliminate genes that in all probability are not going to be the answer to our experimental hypothesis and
#(2) reduce the number of tests to perform in differential expression analysis.

#About this non-specific filter:
#  - We eliminate the chip control probesets that are not of interest to the researcher: they start with AFFX (feature.exclude="^AFFX")
#  - As a measure of dispersion for filtering by variance we use IQR (interquartile range) (var.func=IQR) and we choose the 0.5 quartile of the IQR values as threshold (var.cutoff=0.5). Although a less aggressive filter could be applied
#  - We also demand that probesets that do not have an Entrez Gene ID are not eliminated (they can help build the statistical model of differential expression if they meet the above values)

celfiles.rma_filtered<-nsFilter(celfiles.rma, require.entrez=FALSE, var.func=IQR, var.cutoff=0.5, feature.exclude="^AFFX")$eset
exprdata_filtered<-exprs(celfiles.rma_filtered)
dim(exprdata_filtered)

#The matrix has a dimension of 10087 rows and 31 columns. We have removed 54675-10087 = 44588 probesets


#We now proceed with the different steps of differential expression analysis. Remember that we want to identify those genes that
#significantly, their expression levels change between the experimental conditions: post vs. pre.

#Each row represents an array and the value 1 in the column tells us what type of sample the array belongs to. We are interested in 'Group'.

samples <- celfiles$Group
samples

#We store the phenotypic information as a 'factor' data type

samples <- as.factor(samples)

#Design matrix:

design <- model.matrix(~0+samples) 
design

#We provide the header of the columns with more intuitive names
colnames(design) <- c("Pos", "Pre")
design


#Fitting the linear model to each gene of the array set using the pre-processed and filtered data and the design matrix

#We fit a linear model to each gene
fit = lmFit(celfiles.rma_filtered, design)


#Design of the contrast matrix taking into account the contrast to be studied (post - pre) and the design matrix

#We want to study the differential expression between arrays with post-pre conditions.
#To do this, we designed the following matrix of contrasts:

contrast.matrix = makeContrasts(
  pos_pre = Pos - Pre, 
  levels = design)
contrast.matrix


#calculate the estimate of the coefficients and the errors given the fitted linear model and the contrast matrix

fit.cont <- contrasts.fit(fit, contrast.matrix)


#Calculation of the moderate t-statistic and probabilities of differential expression following a Bayesian model given the estimation of the coefficients and errors above

res.limma <- eBayes(fit.cont)

#Generation of the list of the most differentially expressed genes/probesets ordered by p-value and according to the previous restrictions

#The genes/probesets must have an adjusted p-value less than 0.05 (adjustment of the p-values by the FDR Benjamini Hochberg method)
#according to the statistical test used by limma and must be over-expressed or inhibited by at least a factor of 2.

results1 <- topTable(res.limma, p.value=0.05, adjust.method="BH", sort.by="p",number=nrow(res.limma),lfc=1)

dim(results1)
head(results1)

sobreexpr1<-results1[results1$logFC>0,]

inhibidos1<-results1[results1$logFC<0,]

dim(sobreexpr1)
dim(inhibidos1)

#Generation of the list of the most differentially expressed genes/probesets ordered by p-value and according to the restrictions described

#Change the adjustment method of the p-values of the command executed in section 2 above (FDR Benjamini Hochberg) by the Bonferroni adjustment method.

results2 <- topTable(res.limma, p.value=0.05, adjust.method="bonferroni", sort.by="p",number=nrow(res.limma),lfc=1)

dim(results2)
results2

#Generation of the list of the most differentially expressed genes/probesets ordered by p-value and according to the restrictions described

#genes/probesets must have an adjusted p-value less than 0.01 (adjustment of p-values by the FDR Benjamini Hochberg method)
#according to the statistical test used by limma and must be over-expressed or inhibited by at least a factor of 4.

results3 <- topTable(res.limma, p.value=0.01, adjust.method="BH", sort.by="p",number=nrow(res.limma),lfc=2)

dim(results3)
results3

sobreexpr3<-results3[results3$logFC>0,]

inhibidos3<-results3[results3$logFC<0,]

dim(sobreexpr3)
dim(inhibidos3)

#High-level analysis: annotation
#IMPORTANT: make the annotation SEPARATELY
#make an annotation with the list of genes according to section 2 of the problem and another annotation with the list of genes obtained

probesets<-rownames(results1)

#symbol of the gene it represents, using the annotation of the hgu133plus2 chip

symbols_list <- mget(probesets,hgu133plus2SYMBOL)

#We go from list type to character type
symbols<-unlist(symbols_list)

results1 <- cbind(results1, symbols)

probesets<-rownames(results3)

#symbol of the gene it represents, using the annotation of the hgu133plus2 chip

symbols_list <- mget(probesets,hgu133plus2SYMBOL)

#We go from list type to character type
symbols<-unlist(symbols_list)

results3 <- cbind(results3, symbols)

results3


#Annotate each probeset with its gene description

#description of the gene it represents, using the annotation of the hgu133plus2 chip
genename_list <-mget(probesets,hgu133plus2GENENAME)

#To get the name of the gene

genename<-unlist(genename_list)

results1 <- cbind(results1, genename)
head(results1)

#description of the gene it represents, using the annotation of the hgu133plus2 chip
genename_list <-mget(probesets,hgu133plus2GENENAME)

#To get the name of the gene

genename<-unlist(genename_list)

results3 <- cbind(results3, genename)
results3


#Annotation of each probeset with its KEGG PATHWAY using the auxiliary file get_info_kegg.R

#Pathway(s) of the gene it represents, using the annotation of the hgu133plus2 chip

pathway_list <- mget(probesets,hgu133plus2PATH)

#We obtain information from the pathways through the function get_info_kegg (get_info_kegg.R)
Pathway<-sapply(pathway_list,get_info_kegg)

results1 <- cbind(results1, Pathway)
head(results1)

#Pathway(s) of the gene it represents, using the annotation of the hgu133plus2 chip

pathway_list <- mget(probesets,hgu133plus2PATH)

#We obtain information from the pathways through the function get_info_kegg (get_info_kegg.R)
Pathway<-sapply(pathway_list,get_info_kegg)

results3 <- cbind(results3, Pathway)
results3

#Save differentially expressed and annotated genes in two tabular files
## - pos_vs_pre_005.txt: file containing the differentially expressed genes including the annotation information of the post vs pre comparison with the first restrictions indicated
## - pos_vs_pre_001.txt: file containing the differentially expressed genes including the annotation information of the post vs pre comparison with the second restrictions indicated

#pos_vs_pre_005.txt

write.table(results1,file='pos_vs_pre_005.txt',sep='\t',quote=FALSE,col.names=NA)

#pos_vs_pre_001.txt

write.table(results3,file='pos_vs_pre_001.txt',sep='\t',quote=FALSE,col.names=NA)
