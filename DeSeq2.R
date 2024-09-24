install.packages(c("htmltools","conflicted"))
BiocManager::install(c("DESeq2","GEOquery","knitr","biomaRt","org.Hs.eg.db","AnnotationDbi"))

library(conflicted)
library(GEOquery)
library(DESeq2)
library(ggplot2)
library(htmltools)
library(dplyr)
library(tidyverse)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("C:/Users/Maya/OneDrive - University of Calgary/Maya's folder/Summer students/3featureinvestigation/files")

gse <- getGEO(filename="C:/Users/Maya/OneDrive - University of Calgary/Maya's folder/Summer students/3featureinvestigation/GSE191142_series_matrix.txt.gz", GSEMatrix = FALSE)


listed_files <- list.files("C:/Users/Maya/OneDrive - University of Calgary/Maya's folder/Summer students/3featureinvestigation/files", pattern=".gz")

loaded <- lapply(listed_files, function(x) read.table(x, header=FALSE, sep=" "))

loaded <- lapply(loaded, setNames, (c('geneid', 'count')))

joined <- loaded %>% reduce(inner_join, by='geneid')

joined.drop <- joined[-c(2:7)]

rownames(joined.drop) <- joined$geneid
joined.drop$geneid <-NULL
sample_names=gse$geo_accession[7:34]
colnames(joined.drop) <- c(sample_names) 

mycols=data.frame(gse$geo_accession[7:34], gse$source_name_ch1[7:34])

dds <- DESeqDataSetFromMatrix(
  countData = joined.drop,
  colData = mycols,
  design = ~ gse.source_name_ch1.7.34.)

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)


res <- res[order(res$padj),]
head(res)

par(mfrow=c(2,3))
plotCounts(dds, gene='ENSG00000244734.3', intgroup='gse.source_name_ch1.7.34.')
plotCounts(dds, gene='ENSG00000143546.9', intgroup='gse.source_name_ch1.7.34.')
plotCounts(dds, gene='ENSG00000211592.8', intgroup='gse.source_name_ch1.7.34.')
plotCounts(dds, gene='ENSG00000134352.19', intgroup='gse.source_name_ch1.7.34.')
plotCounts(dds, gene='ENSG00000156467.9', intgroup='gse.source_name_ch1.7.34.')
plotCounts(dds, gene='ENSG00000101825.7', intgroup='gse.source_name_ch1.7.34.')

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup='gse.source_name_ch1.7.34.')

res$fixed_names = sapply(strsplit(rownames(res), ".", fixed=T), function(x) x[1])

IDs <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
             filters = "ensembl_gene_id", values = res$fixed_names,
             mart = mart)

colnames(IDs)[1] = 'fixed_names'

geneidres <- merge(IDs, as.data.frame(res), by='fixed_names', all=TRUE)

colnames(geneidres)[1] = 'ENSEMBL'
colnames(geneidres)[2] = 'gene symbol'

write.csv(geneidres, "C:/Users/Maya/OneDrive - University of Calgary/Maya's folder/Summer students/3featureinvestigation/CTCs.csv", row.names=TRUE)
###



