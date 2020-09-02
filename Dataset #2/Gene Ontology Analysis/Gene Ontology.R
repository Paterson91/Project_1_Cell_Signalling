#install.packages("readxl")
#install.packages("factoextra")
#install.packages("viridis")
#install.packages("rafalib")
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("gplots")
#BiocManager::install("org.Mm.eg.db")
library(rafalib)
library(readxl)
library("RColorBrewer")
library(limma)
library("gplots")
library(factoextra)
library(dplyr)
library(tidyverse)
library(viridis)
library(genefilter)

library(clusterProfiler)
require("biomaRt")
library(org.Mm.eg.db)

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #1 LuciaCellSignalling 21_2_19/Dataset #2/Gene Ontology Analysis/")
full_table = data.matrix(read.csv("../results_Marucci-3.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE))

metadata = read.table("../MetaData_Exp2.txt", sep = "\t", header = TRUE)
metadata

#Samples to compare

set1 = "EXP4"
set2 = "EXP3"

#P-Value Threshold
pvalue = 0.05

grep1 = full_table[,(grep(set1, colnames(full_table),value = TRUE))]
grep2 = grep1[,(grep(set2, colnames(grep1),value = TRUE))]
final = grep2[,(grep("PValue", colnames(grep2),value = TRUE))]
final = as.data.frame(final)
colnames(final)=c("P-Value")
final = tibble::rownames_to_column(final, "ensembl_gene_id")

sig = final[final$`P-Value`<=pvalue,]

genesonly = sig[,1]

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name","entrezgene"),
  filter="ensembl_gene_id",
  values=genesonly,
  uniqueRows=TRUE)

annotLookup <- data.frame(
  genesonly[match(annotLookup$ensembl_gene_id, genesonly)],
  annotLookup)

ggo <- groupGO(gene     = as.character(annotLookup$entrezgene),
                keyType = "ENTREZID",
                OrgDb    = org.Mm.eg.db,
                ont      = "CC",
                level    = 2,
                readable = TRUE)

barplot(ggo, drop=TRUE, showCategory=12, border = c(10,10))


#ColourStrategy = metadata_t$Sample
#ColourStrategy2 = metadata_t$Condition

order1 = as.character(metadata_t$Original_Sample_Name)
iso_table = full_table[, match(order1, colnames(full_table))]
iso_table[,1:ncol(iso_table)] <- sapply(round(iso_table[,1:ncol(iso_table)], 0), as.integer)
head(iso_table)

####################################
#### Heatmap Generation
####################################

vst = varianceStabilizingTransformation(iso_table, blind = TRUE)
select <- vst[order(rowMeans(vst),decreasing=T)[1:100],]
cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy2))]


