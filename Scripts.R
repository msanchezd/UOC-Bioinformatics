# All R-scripts in Dynamic_report_TFM.Rmd

# Setting working directory ----
## IN R-STUDIO:
### Session --> Set working directory --> Choose directory 
## WORKING DIRECTORY (THIS FOLDER):
library(rstudioapi)
current_path <- getActiveDocumentContext()$path 
setwd(dirname(current_path))
getwd() # to show the pathway

######## PART 1: Scripts from Dynamic_report_TFM.Rmd ##########

# Installing package if needed ----
list.of.packages <- c("rstudioapi", "dplyr", "xlsx", "rJava", "gplots",
                      "devtools", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Conneting with Bioconductor ----
source("https://bioconductor.org/biocLite.R")

# Installing & loading Bioconductor packages ----
biocLite()
biocLite("biomaRt")
## If biocLite("biomaRt") do not work:
### LINUX COMAND LINE: sudo apt-get install r-bioc-biomart

# Uploading .csv input file ----
numts_coord <- read.csv(file=params$file1, sep = ",", header = TRUE)
str(numts_coord, vec.len = 1)
head(numts_coord)

# The biomaRt package ----
library("biomaRt")
citation("biomaRt") # Package citation

## Preparing Package ‘biomaRt'
gene_mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", 
                    host="grch37.ensembl.org", 
                    path="/biomart/martservice", 
                    dataset="hsapiens_gene_ensembl")

listMarts(gene_mart) # Ensembl version used

## Creating "All_attributes.txt" and "All_filters.txt" with all options:
All_attributes <- listAttributes(gene_mart)
All_filters <- listFilters(gene_mart)

write.table(All_attributes, file = "All_attributes.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)
write.table(All_filters, file = "All_filters.txt", sep = "\t", 
            row.names = FALSE, quote = FALSE)

# Extracting all genes within NUMTs coordinates ----
library(plyr); library(dplyr)

## Adapting coordenates to download attributes
numts_coord$coord_n <- do.call(paste, c(numts_coord[,2:4], sep = ":"))
numts_vector <-as.vector(t(numts_coord$coord_n))
id <-as.vector(t(numts_coord$id))

## Setting attributes and filters
### Our attributes
attributes_gene = c("chromosome_name", "start_position", "end_position", "strand",
                    "hgnc_symbol", "ensembl_gene_id_version","ensembl_gene_id",
                    "transcript_count")

## Getting values: gene_results.txt (loop)

gene_results <- numeric(0)
i <- 1

for (i in 1:length(numts_vector)) {
  b<-i
  
  gene_results_b = getBM(attributes_gene,
                         filters = c("chromosomal_region"),
                         values = list(chromosomal_region=numts_vector[b]), 
                         mart = gene_mart)
  
  if (length(gene_results_b[,1]) == 0) {
    gene_results <- rbind(gene_results, c(rep("", length(attributes_gene)),
                                          do.call(paste, list(numts_coord[b,1]))))
  } else {
    gene_results_b$id <- do.call(paste, list(numts_coord[b,1]))
    gene_results <- rbind(gene_results, gene_results_b)
  }
  
  i <- i + 1
}

gene_results[gene_results==""]  <- NA 

## Reordering columns (gene_results.txt)  

gene_results <- gene_results %>% dplyr::select("id", everything())

## Sorting results (gene_results.txt) 
gene_results <- gene_results[order(gene_results$id, 
                                   gene_results$chromosome_name,
                                   gene_results$start_position),]

## Saving the results (gene_results.txt)
write.table(gene_results, file = "gene_results.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

## Uploading intermediate documents
gene_results <- read.table("gene_results.txt", header = TRUE, sep = "\t")
head(gene_results)
str(gene_results, vec.len = 1)


# Extracting all genes upstream and downstream from the NUMTs coordinates ----
library(plyr); library(dplyr)

## Indicating new coordinates
up_start <- numts_coord[3] - 1000 
up_end <- numts_coord[3] - 100 
down_start <- numts_coord[4] + 100 
down_end <- numts_coord[4] + 1000

up_numts_coord <- data.frame(numts_coord$id, 
                             numts_coord$chr, 
                             up_start, 
                             up_end)
down_numts_coord <- data.frame(numts_coord$id, 
                               numts_coord$chr, 
                               down_start, 
                               down_end)

## Adapting coordenates to download attributes
up_numts_coord$coord_n <- do.call(paste, c(up_numts_coord[,2:4], sep = ":"))
up_numts_vector <-as.vector(t(up_numts_coord$coord_n))
down_numts_coord$coord_n <- do.call(paste, c(down_numts_coord[,2:4], sep = ":"))
down_numts_vector <-as.vector(t(down_numts_coord$coord_n))


## Getting values: up_gene_results.txt (loop)

up_gene_results <- numeric(0)
i <- 1

for (i in 1:length(up_numts_vector)) {
  b<-i
  
  up_gene_results_b = getBM(attributes_gene,
                            filters = c("chromosomal_region"),
                            values = list(chromosomal_region=up_numts_vector[b]), 
                            mart = gene_mart)
  
  if (length(up_gene_results_b[,1]) == 0) {
    up_gene_results <- rbind(up_gene_results, 
                             c(rep("", length(attributes_gene)), 
                               do.call(paste, list(numts_coord[b,1]))))
  } else {
    up_gene_results_b$id <- do.call(paste, list(numts_coord[b,1]))
    up_gene_results <- rbind(up_gene_results, up_gene_results_b)
  }
  
  i <- i + 1
}

up_gene_results[up_gene_results==""]  <- NA 

### Reordering columns (up_gene_results.txt) 
up_gene_results <- up_gene_results %>%  dplyr::select("id", everything())

### Sorting results (up_gene_results.txt)
up_gene_results <- up_gene_results[order(up_gene_results$id,
                                         up_gene_results$chromosome_name,
                                         up_gene_results$start_position),]

### Saving the results (up_gene_results.txt)
write.table(up_gene_results, file = "up_gene_results.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)


## Getting values: down_gene_results.txt (loop)

down_gene_results <- numeric(0)
i <- 1

for (i in 1:length(down_numts_vector)) {
  b<-i
  
  down_gene_results_b = getBM(attributes_gene,
                              filters = c("chromosomal_region"),
                              values = list(chromosomal_region=down_numts_vector[b]), 
                              mart = gene_mart)
  
  if (length(down_gene_results_b[,1]) == 0) {
    down_gene_results <- rbind(down_gene_results, 
                               c(rep("", length(attributes_gene)), 
                                 do.call(paste, list(numts_coord[b,1]))))
  } else {
    down_gene_results_b$id <- do.call(paste, list(numts_coord[b,1]))
    down_gene_results <- rbind(down_gene_results, down_gene_results_b)
  }
  
  i <- i + 1
}

down_gene_results[down_gene_results==""]  <- NA 

### Reordering columns (down_gene_results.txt) 
down_gene_results <- down_gene_results %>% dplyr::select("id", everything())

### Sorting results (down_gene_results.txt)
down_gene_results <- down_gene_results[order(down_gene_results$id, 
                                             down_gene_results$chromosome_name,
                                             down_gene_results$start_position),]

### Saving the results (down_gene_results.txt)
write.table(down_gene_results, file = "down_gene_results.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)


## Uploading intermediate documents
up_gene_results <- read.table("up_gene_results.txt", header = TRUE, sep = "\t")
str(up_gene_results, vec.len = 1)

down_gene_results <- read.table("down_gene_results.txt", header = TRUE, sep = "\t")
str(down_gene_results, vec.len = 1)

# Extracting genes originated by NUMT insertions ----
## Total genes
TOTAL_GENES <- as.character(na.omit(unique(gene_results$ensembl_gene_id_version)))
UP_GENES <- as.character(na.omit(unique(up_gene_results$ensembl_gene_id_version)))
DOWN_GENES <- as.character(na.omit(unique(down_gene_results$ensembl_gene_id_version)))

## Common genes present in all three list of genes  
## (gene_results, upstream and downstream)
library(plyr); library(dplyr)
data_joined <- dplyr::inner_join(up_gene_results, down_gene_results)

## Creating column "Localization"
data_joined$localization <- data_joined$ensembl_gene_id
data_joined$localization[!is.na(data_joined$localization)] <- "intronic"
data_joined$localization[is.na(data_joined$localization)] <- "intergenic"

## Genes in gene_results.txt but not in upper and downstream regions
int_gene_results <- anti_join(gene_results, data_joined)
str(int_gene_results, vec.len = 1)
summary(int_gene_results[c(1,2,6,7,8)])

#### START of exclusive from my dataset (TFM) ####

## Excluding the gene ARHGAP15

for (i in 1:nrow(int_gene_results)) {
  if (!is.na(int_gene_results$hgnc_symbol[i])) {
    if (int_gene_results$hgnc_symbol[i] == "ARHGAP15") {
      int_gene_results[i,c(2:ncol(int_gene_results))] <- NA
    } else {
      int_gene_results$hgnc_symbol[i] <- int_gene_results$hgnc_symbol[i]
    }
    i <- i + 1
  } else {
    int_gene_results$hgnc_symbol[i] <- NA
  }
}
#### END ####

## Saving genes.txt
genes <- as.character(na.omit(unique(int_gene_results$ensembl_gene_id)))

write.table(genes, file="genes.txt",col.names = F, sep="\t",quote=F,row.names=F)

## Uploading intermediate documents
genes <- read.table("genes.txt", header = FALSE, sep = "\t")

table_numts_genes <- gene_results[gene_results$ensembl_gene_id
                                  %in% genes$V1,]

chrom <- (unique(table_numts_genes[c(2,7)]))

# Number of genes per chromosome
table(chrom[1])

# GO terms (go_results.txt) ----
## Setting attributes and filters
### Our attributes
attributes_go = c("hgnc_symbol", "ensembl_gene_id_version",
                  "go_id", "name_1006", "definition_1006")

go_results = getBM(attributes_go,
                   filters = c("ensembl_gene_id"),
                   values = list(ensembl_gene_id=genes$V1), 
                   mart = gene_mart)

go_results[go_results==""]  <- NA 

go_results <- go_results[order(go_results$ensembl_gene_id_version, go_results$go_id),]

go_results <- go_results[complete.cases(go_results$name_1006), ]

write.table(go_results, file = "go_results.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

## Uploading intermediate documents
go_results <-read.table("go_results.txt", header = TRUE, sep = "\t")
go_results

## Plotting GO results 
par(mfrow=c(1,2))
par(mar = c(8.5, 2.5, 2.5, 4), xpd=TRUE)

ensembl_go <- na.omit(unique(go_results[c("hgnc_symbol", "name_1006")]))

ensembl_go <- ensembl_go[complete.cases(ensembl_go), ]

colors = c("aquamarine3","yellow2","azure3", 
           "darkgoldenrod1", "lawngreen", "plum",
           "gray9","deeppink1","cornflowerblue", 
           "antiquewhite3", "slategrey", "tomato")

### Bar plot
barplot(table(ensembl_go$name_1006), las=2, cex.main = 1.2,
        cex.axis = 0.7, cex = 0.7,  col = colors)

### pie chart
counts = table(ensembl_go$name_1006)  ## get counts
labs = paste(levels(ensembl_go$name_1006), counts)  ## create labels
pie(counts, labels = labs, col = colors, cex=0.5)  ## plot
legend("bottom", inset=c(0,-0.6), labs, cex=0.6, fill=colors)

# phenotype_results.txt ----
## Setting attributes and filters
### Our attributes
attributes_phenotype = c("hgnc_symbol", "ensembl_gene_id_version", "transcript_count",
                         "gene_biotype", "description")

## Getting values: phenotype_results ----

phenotype_results = getBM(attributes_phenotype,
                          filters = c("ensembl_gene_id"),
                          values = list(ensembl_gene_id=genes$V1),
                          mart = gene_mart)

# class(phenotype_results) # data.frame
phenotype_results[phenotype_results==""]  <- NA 

phenotype_results <- phenotype_results[order(phenotype_results$ensembl_gene_id_version),]

write.table(phenotype_results, file = "phenotype_results.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

## Uploading intermediate documents
phenotype_results <-read.table("phenotype_results.txt", header = TRUE, sep = "\t")
summary(phenotype_results)

## Plotting Gene biotype results 
par(mfrow=c(1,2))
par(mar = c(6, 2.5, 2.5, 2.5), xpd=TRUE)

ensembl_biotype <- na.omit(unique(phenotype_results[c("ensembl_gene_id_version", "gene_biotype")]))

colors = c("yellow2","orchid3","orangered3", 
           "olivedrab3", "lightskyblue3", "plum")

### Bar plot
barplot(table(ensembl_biotype$gene_biotype), las=2, cex.main = 1.2, 
        cex.axis = 0.7, cex = 0.7,  col = colors)

### pie chart
counts = table(ensembl_biotype$gene_biotype)  ## get counts
labs = paste(levels(ensembl_biotype$gene_biotype), counts)  ## create labels
pie(counts, labels = labs, col = colors, cex=0.5)  ## plot
legend("bottom", inset=c(0, -0.2), labs, cex=0.6, fill=colors)

# EXPRESSION STUDY ----
## Uploading GTEx means in TPM
GTEx_mean_tpm <- read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", 
                            skip = 2, header = TRUE, sep = "\t")

## Creating mean_tpm_GTEx.txt
library(plyr); library(dplyr)

GTEx_tpm <- GTEx_mean_tpm
colnames(GTEx_tpm)[1] <- "GTEx_gene_id_version"

gene_id <- GTEx_mean_tpm$gene_id
GTEx_genes <- numeric(0)
for (i in 1:length(GTEx_mean_tpm$gene_id)){
  x <- unlist(strsplit(as.character(GTEx_mean_tpm[i,1]), split='.', fixed=TRUE))[1]
  GTEx_genes <- rbind(GTEx_genes, x)
}

GTEx_mean_tpm$gene_id <- GTEx_genes

mean_tpm_fromGTEx <- numeric(0)
for (i in 1:nrow(genes)){
  y <- subset(GTEx_mean_tpm, gene_id == genes[i,1])
  mean_tpm_fromGTEx <- rbind(mean_tpm_fromGTEx, y)
}

tissue_means <- rowMeans(mean_tpm_fromGTEx[,3:length(mean_tpm_fromGTEx)])
mean_tpm_fromGTEx$tissue_means <- tissue_means

mean_tpm_fromGTEx$sum <- rowSums(mean_tpm_fromGTEx[,3:length(mean_tpm_fromGTEx)])

mean_tpm_fromGTEx$gene_id <- as.character(mean_tpm_fromGTEx$gene_id)

mean_tpm_GTEx <- dplyr::inner_join(mean_tpm_fromGTEx,
                                   GTEx_tpm)

mean_tpm_GTEx <- mean_tpm_GTEx %>% dplyr::select("gene_id", "GTEx_gene_id_version", everything())
colnames(mean_tpm_GTEx)[1] <- "ensembl_gene_id"
str(mean_tpm_GTEx)

## saving results
write.table(mean_tpm_GTEx, file = "mean_tpm_GTEx.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

## Uploading intermediate documents
mean_tpm_GTEx <- read.table("mean_tpm_GTEx.txt", header = TRUE, 
                            sep = "\t", dec = ".")

# subset_expressed.txt ----
## Creating subset of genes with >= 0.5 TPM (expressed)
GTEx <- mean_tpm_GTEx
subset_expressed <- subset(GTEx, GTEx[4] >= 0.5 | GTEx[5] >= 0.5 | GTEx[6] >= 0.5 |
                             GTEx[7] >= 0.5 | GTEx[8] >= 0.5 | GTEx[9] >= 0.5 | GTEx[10] >= 0.5 |
                             GTEx[11] >= 0.5 | GTEx[12] >= 0.5 | GTEx[13] >= 0.5 | GTEx[14] >= 0.5 |
                             GTEx[15] >= 0.5 | GTEx[16] >= 0.5 | GTEx[17] >= 0.5 | GTEx[18] >= 0.5 |
                             GTEx[19] >= 0.5 | GTEx[20] >= 0.5 | GTEx[21] >= 0.5 | GTEx[22] >= 0.5 |
                             GTEx[23] >= 0.5 | GTEx[24] >= 0.5 | GTEx[25] >= 0.5 | GTEx[26] >= 0.5 |
                             GTEx[27] >= 0.5 | GTEx[28] >= 0.5 | GTEx[29] >= 0.5 | GTEx[30] >= 0.5 |
                             GTEx[31] >= 0.5 | GTEx[32] >= 0.5 | GTEx[33] >= 0.5 | GTEx[34] >= 0.5 | 
                             GTEx[35] >= 0.5 | GTEx[36] >= 0.5 | GTEx[37] >= 0.5 | GTEx[38] >= 0.5 | 
                             GTEx[39] >= 0.5 | GTEx[40] >= 0.5 | GTEx[41] >= 0.5 | GTEx[42] >= 0.5 |
                             GTEx[43] >= 0.5 | GTEx[44] >= 0.5 | GTEx[45] >= 0.5 | GTEx[46] >= 0.5 |
                             GTEx[47] >= 0.5 | GTEx[48] >= 0.5 | GTEx[49] >= 0.5 | GTEx[50] >= 0.5 | 
                             GTEx[51] >= 0.5 | GTEx[52] >= 0.5 | GTEx[53] >= 0.5 | GTEx[54] >= 0.5 | 
                             GTEx[55] >= 0.5 | GTEx[56] >= 0.5 | GTEx[57] >= 0.5 )

## saving results
write.table(subset_expressed, file = "subset_expressed.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)


## Uploading intermediate documents
subset_expressed <- read.table("subset_expressed.txt", header = TRUE, 
                               sep = "\t", dec = ".")
summary(subset_expressed[c(4:ncol(subset_expressed))])

## Heatmap of all expressed genes.
library(gplots)
par(oma=c(10,4,4,2))
subset_mean_tpm <- subset_expressed
subset_mean_tpm[1] <- NULL
rownames(subset_mean_tpm) <- subset_mean_tpm$Description
subset_mean_tpm[1] <- NULL
subset_mean_tpm[1] <- NULL

heatmap.2(data.matrix(subset_mean_tpm[1:53]), trace='none', scale = "row", 
          cexRow=0.6, cexCol = 0.6)

## Graphical representation: expression profile of
## hight expressed genes (>=1000 TMP)
par(mar = c(10.5, 4, 4, 7.5), xpd=TRUE)
GTEx1 <- subset_mean_tpm
subset_mean_tpm2 <- subset(GTEx1, GTEx1[1] >= 1000.0 |GTEx1[2] >= 1000.0 |
                             GTEx1[3] >= 1000.0 | GTEx1[4] >= 1000.0 | 
                             GTEx1[5] >= 1000.0 | GTEx1[6] >= 1000.0 | 
                             GTEx1[7] >= 1000.0 | GTEx1[8] >= 1000.0 | 
                             GTEx1[9] >= 1000.0 | GTEx1[10] >= 1000.0 | 
                             GTEx1[11] >= 1000.0 | GTEx1[12] >= 1000.0 | 
                             GTEx1[13] >= 1000.0 | GTEx1[14] >= 1000.0 |
                             GTEx1[15] >= 1000.0 | GTEx1[16] >= 1000.0 | 
                             GTEx1[17] >= 1000.0 | GTEx1[18] >= 1000.0 | 
                             GTEx1[19] >= 1000.0 | GTEx1[20] >= 1000.0 | 
                             GTEx1[21] >= 1000.0 | GTEx1[22] >= 1000.0 | 
                             GTEx1[23] >= 1000.0 | GTEx1[24] >= 1000.0 | 
                             GTEx1[25] >= 1000.0 | GTEx1[26] >= 1000.0 |
                             GTEx1[27] >= 1000.0 | GTEx1[28] >= 1000.0 | 
                             GTEx1[29] >= 1000.0 | GTEx1[30] >= 1000.0 | 
                             GTEx1[31] >= 1000.0 | GTEx1[32] >= 1000.0 | 
                             GTEx1[33] >= 1000.0 | GTEx1[34] >= 1000.0 | 
                             GTEx1[35] >= 1000.0 | GTEx1[36] >= 1000.0 | 
                             GTEx1[37] >= 1000.0 | GTEx1[38] >= 1000.0 | 
                             GTEx1[39] >= 1000.0 | GTEx1[40] >= 1000.0 | 
                             GTEx1[41] >= 1000.0 | GTEx1[42] >= 1000.0 | 
                             GTEx1[43] >= 1000.0 | GTEx1[44] >= 1000.0 | 
                             GTEx1[45] >= 1000.0 | GTEx1[46] >= 1000.0 | 
                             GTEx1[47] >= 1000.0 | GTEx1[48] >= 1000.0 | 
                             GTEx1[49] >= 1000.0 | GTEx1[50] >= 1000.0 | 
                             GTEx1[51] >= 1000.0 | GTEx1[52] >= 1000.0 | 
                             GTEx1[53] >= 1000.0 )
matplot(t(data.matrix(subset_mean_tpm2[1:53])), type = "b", 
        col = c(1:ncol(subset_mean_tpm2)), 
        cex.main = 1, cex.lab = 0.8,	ylab = "TPM", pch=c(15:18,21:25), 
        main = "TPM/tisse", axes = FALSE)
axis(2, cex.axis=0.7)
axis(side=1,at=1:ncol(subset_mean_tpm2[1:53]), cex.axis=0.6, las = 2,
     labels=colnames(subset_mean_tpm2[1:53]))

legend("right", inset=c(-0.25, 1), legend=rownames(subset_mean_tpm2[1:53]), 
       col=c(1:ncol(subset_mean_tpm2)),pch= c(15:18,21:25), cex = 0.6,
       bg= ("white"), horiz=F)

## Graphical representation: expression profile of
## medium expressed genes (between 10 and 1000 TMP)

subset_mean_tpm3 <- subset(GTEx1, GTEx1[1] >= 10.0 |GTEx1[2] >= 10.0 |
                             GTEx1[3] >= 10.0 | GTEx1[4] >= 10.0 |
                             GTEx1[5] >= 10.0 | GTEx1[6] >= 10.0 | 
                             GTEx1[7] >= 10.0 | GTEx1[8] >= 10.0 | 
                             GTEx1[9] >= 10.0 | GTEx1[10] >= 10.0 | 
                             GTEx1[11] >= 10.0 | GTEx1[12] >= 10.0 | 
                             GTEx1[13] >= 10.0 | GTEx1[14] >= 10.0 |
                             GTEx1[15] >= 10.0 | GTEx1[16] >= 10.0 | 
                             GTEx1[17] >= 10.0 | GTEx1[18] >= 10.0 | 
                             GTEx1[19] >= 10.0 | GTEx1[20] >= 10.0 | 
                             GTEx1[21] >= 10.0 | GTEx1[22] >= 10.0 | 
                             GTEx1[23] >= 10.0 | GTEx1[24] >= 10.0 | 
                             GTEx1[25] >= 10.0 | GTEx1[26] >= 10.0 |
                             GTEx1[27] >= 10.0 | GTEx1[28] >= 10.0 | 
                             GTEx1[29] >= 10.0 | GTEx1[30] >= 10.0 | 
                             GTEx1[31] >= 10.0 | GTEx1[32] >= 10.0 | 
                             GTEx1[33] >= 10.0 | GTEx1[34] >= 10.0 | 
                             GTEx1[35] >= 10.0 | GTEx1[36] >= 10.0 | 
                             GTEx1[37] >= 10.0 | GTEx1[38] >= 10.0 | 
                             GTEx1[39] >= 10.0 | GTEx1[40] >= 10.0 | 
                             GTEx1[41] >= 10.0 | GTEx1[42] >= 10.0 | 
                             GTEx1[43] >= 10.0 | GTEx1[44] >= 10.0 | 
                             GTEx1[45] >= 10.0 | GTEx1[46] >= 10.0 | 
                             GTEx1[47] >= 10.0 | GTEx1[48] >= 10.0 | 
                             GTEx1[49] >= 10.0 | GTEx1[50] >= 10.0 | 
                             GTEx1[51] >= 10.0 | GTEx1[52] >= 10.0 | 
                             GTEx1[53] >= 10.0 ) 
subset_mean_tpm3 <- subset_mean_tpm3[!rownames(subset_mean_tpm3) %in% rownames(subset_mean_tpm2), ]

par(mar = c(10.5, 4, 4, 7.5), xpd=TRUE)
matplot(t(data.matrix(subset_mean_tpm3[1:53])), type = "b", col = c(1:ncol(subset_mean_tpm3)), 
        cex.main = 1, cex.lab = 0.8,	ylab = "TPM", pch=c(15:18,21:25), 
        main = "TPM/tisse", axes = FALSE)
axis(2, cex.axis=0.7)
axis(side=1,at=1:ncol(subset_mean_tpm3[1:53]), cex.axis=0.6, las = 2,
     labels=colnames(subset_mean_tpm3[1:53]))

legend("right", inset=c(-0.25, 1), legend=rownames(subset_mean_tpm3[1:53]), 
       col=c(1:ncol(subset_mean_tpm3)),pch= c(15:18,21:25), cex = 0.6, bg= ("white"), horiz=F)

# FINAL TABLE: FINAL_OUTPUT_TABLE.txt ----
library(plyr); library(dplyr)

## Cheking genes and creating table
table_numts_genes <- gene_results[gene_results$ensembl_gene_id
                                  %in% genes$V1,]

str(table_numts_genes)
length(unique(table_numts_genes$ensembl_gene_id))
summary(table_numts_genes)

table_numts_genes <- dplyr::full_join(numts_coord, 
                                      table_numts_genes)

table_numts_genes <- dplyr::full_join(data_joined[c("id", "localization")], 
                                      table_numts_genes)
table_numts_genes$localization[is.na(table_numts_genes$localization)] <- "partial_gene"
table_numts_genes <- dplyr::full_join(phenotype_results[c("ensembl_gene_id_version", 
                                                          "gene_biotype", "description")], 
                                      table_numts_genes)

table_numts_go <- dplyr::full_join(go_results[c("ensembl_gene_id_version",
                                                "name_1006")], 
                                   table_numts_genes)

table_numts_exp <- dplyr::full_join(mean_tpm_GTEx, 
                                    table_numts_go)

table_numts <- table_numts_exp %>% dplyr::select("id", "localization", "chr",
                                                 "start_n", "end_n", "mt",
                                                 "start_mt", "end_mt",
                                                 "hgnc_symbol", "Description",
                                                 "gene_biotype", "name_1006",
                                                 "ensembl_gene_id_version", 
                                                 "transcript_count",
                                                 "tissue_means", "sum",
                                                 everything())

## Removing columns
table_numts$ensembl_gene_id = NULL

## Sorting results (gene_results.txt) ----

table_numts <- table_numts[order(table_numts$id,
                                 table_numts$chr,
                                 table_numts$start_n,
                                 table_numts$start_position),]

table_numts <- table_numts[!duplicated(table_numts), ]

## Saving results
write.table(table_numts, file="FINAL_OUTPUT_TABLE.txt",
            sep="\t",quote=F,row.names=F)

## Uploading intermediate documents
table_numts <- read.table("FINAL_OUTPUT_TABLE.txt", header = TRUE, sep = "\t")

## Showing 6 first data from FINAL TABLE "table_numts.txt"
head(table_numts)

# ADDITIONAL GRAPHS WITH EXPRESSED DATA ------------
#fig.cap= "Gene biotype and GO term annotated for our list of expressed genes.
library(plyr); library(dplyr)

only_expressed <- table_numts[table_numts$ensembl_gene_id
                              %in% subset_expressed$GTEx_gene_id_version,]
## Plotting Gene biotype results 
par(mfrow=c(1,2))
par(mar = c(6, 2.5, 2.5, 2.5), xpd=TRUE)

ensembl_biotype_ex <- na.omit(unique(only_expressed[c("ensembl_gene_id_version", 
                                                      "gene_biotype")]))

ensembl_go_ex <- na.omit(unique(only_expressed[c("ensembl_gene_id_version", 
                                                 "name_1006")]))
colors = c("yellow2","orchid3","orangered3", 
           "olivedrab3", "lightskyblue3", "plum")

### pie chart
counts = table(ensembl_biotype_ex$gene_biotype)  ## get counts
labs = paste(levels(ensembl_biotype$gene_biotype), counts)  ## create labels
pie(counts, labels = labs, col = colors, cex=0.5)  ## plot
legend("bottom", inset=c(0, -0.2), labs, cex=0.6, fill=colors)

## Plotting GO results 

ensembl_go <- ensembl_go_ex[complete.cases(ensembl_go_ex), ]

colors = c("aquamarine3","yellow2","azure3", 
           "darkgoldenrod1", "lawngreen", "plum",
           "gray9","deeppink1","cornflowerblue", 
           "antiquewhite3", "slategrey", "tomato")

### pie chart
counts = table(ensembl_go$name_1006)  ## get counts
labs = paste(levels(ensembl_go$name_1006), counts)  ## create labels
pie(counts, labels = labs, col = colors, cex=0.5)  ## plot
legend("bottom", inset=c(0, -0.2), labs, cex=0.6, fill=colors)

######## PART 2: Additional scripts for TFM ##########

#install.packages("overlap")
#install.packages("DescTools")
library(DescTools)
library(plyr)
library(dplyr)
library(biomaRt)

# First, we will check general data from our analysis:
table <- gene_results[gene_results$ensembl_gene_id
                      %in% genes$V1,]
# Genes per NUMT
numt_t <- table(table[1])

# NUMTS WITHOUT GENES
sum(numt_t == 0)

x <- numeric(0) 
# Loop for the rest: 
for (i in 0:length(unique(numt_t))){
  x[i]  = paste("NUMTs containing ", i, " genes:",  sum(numt_t == i))
}
x

# NUMTs per gene
numt_g <- table(table_numts_genes[7])

#Genes deleted after filtering:
sum(numt_g == 0)

y <- numeric(0)
# Loop for the rest: 
for (j in 0:length(unique(numt_g))){
  y[j]  = paste("Genes included in ", j, " NUMTs:",  sum(numt_g == j))
}
y

# File 12: "gene_result_mt.txt" -----------
# Input file ----
numts_coord <- read.csv("NUMTs_coord.csv", sep = ",", header = TRUE)
numts_coord$coord_mt <- do.call(paste, c(numts_coord[,5:7], sep = ":"))
numts_vector_mt <-as.vector(t(numts_coord$coord_mt))

attributes_mt = c("chromosome_name", "start_position", "end_position", "strand", 
                  "hgnc_symbol", "ensembl_gene_id", "transcript_count", "gene_biotype")
gene_results_mt <- numeric(0)
i <- 1

for (i in 1:length(numts_vector_mt)) {
  b<-i
  
  gene_results_b = getBM(attributes_mt,
                         filters = c("chromosomal_region"),
                         values = list(chromosomal_region=numts_vector_mt[b]), 
                         mart = gene_mart)
  
  if (length(gene_results_b[,1]) == 0) {
    gene_results_mt <- rbind(gene_results_mt, c(rep("", length(attributes_mt)),
                                                do.call(paste, list(numts_coord[b,1]))))
  } else {
    gene_results_b$id <- do.call(paste, list(numts_coord[b,1]))
    gene_results_mt <- rbind(gene_results_mt, gene_results_b)
  }
  
  i <- i + 1
}

str(gene_results_mt)
class(gene_results_mt) # data.frame

gene_results_mt[gene_results_mt==""]  <- NA 

## Reordering columns ---- 

gene_mt <- gene_results_mt %>% dplyr::select("id", "hgnc_symbol", everything())

## Sorting results  ----
gene_results_mt <- gene_mt[order(gene_results_mt$id,
                                 gene_results_mt$start_position),]
head(gene_results_mt)
ncol(gene_results_mt) # columns: 9
nrow(gene_results_mt) # rows: 3954
length(unique(gene_results_mt$hgnc_symbol)) # num. mito genes: 38

## Saving the results (gene_results_mt.txt) ----
write.table(gene_results_mt, file = "gene_results_mt.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

############################################################

# IDENTIFICATION OF CORRESPONDING MITOCHONDRIAL GENE FOR EACH NEW
# NUCLEAR GENE WITHIN NUMTs
## Ordering
table_mito_genes <- full_join(numts_coord[c("id", "mt", 
                                            "start_mt", "end_mt")], 
                              gene_results_mt)

table_mito <- table_mito_genes %>% dplyr::select("id", "mt", "start_mt", "end_mt", 
                                                 "hgnc_symbol", "ensembl_gene_id",
                                                 "transcript_count", "gene_biotype",
                                                 everything())


# CONDITION: 
## if gene is strand is +1: 
### START: (gene start - NUMT start) and END: (gene end - NUMT start)
## if gene is strand -1: 
###START: (NUMT end - gene end) and END: (NUMT end - gene start)
str(table_mito)
table_mito$start_mt <- as.numeric(table_mito$start_mt)
table_mito$start_position <- as.numeric(table_mito$start_position)
table_mito$end_mt <- as.numeric(table_mito$end_mt)
table_mito$end_position <- as.numeric(table_mito$end_position)

table_mito$start_ref = NULL
table_mito$end_ref = NULL
for (i in 1:nrow(table_mito)) {
  if (!is.na(table_mito$strand[i])) {
    if (table_mito$strand[1] == 1) {
      table_mito$start_ref[i] <- table_mito$start_position[i] - table_mito$start_mt[i] 
      table_mito$end_ref[i] <- table_mito$end_position[i] - table_mito$start_mt[i]
    } else {
      table_mito$start_ref[i] <- table_mito$end_mt[i] - table_mito$end_position[i]
      table_mito$end_ref[i] <- table_mito$end_mt[i] - table_mito$start_position[i]
    }
    i <- i + 1
  } else {
    table_mito$start_ref[i] <- NA
    table_mito$end_ref[i] <- NA
  }
}
table_mito$start_ref <- as.numeric(table_mito$start_ref)
table_mito$end_ref <- as.numeric(table_mito$end_ref)
str(table_mito)
head(table_mito)
table_mito$from_numt <- NULL
table_mito$from_numt <-  substr(table_mito$id, 6,11)
table_mito$from_numt <- gsub("\\.", "", table_mito$from_numt)
table_mito$from_numt <- gsub("\\X", "23", table_mito$from_numt)
table_mito$from_numt <- gsub("\\Y", "24", table_mito$from_numt)
table_mito$from_numt <- as.numeric(table_mito$from_numt)

# In reference interval, replace negative numbers by 0
table_mito$start_ref[table_mito$start_ref<0] <- 0
table_mito$end_ref[table_mito$end_ref<0] <- 0

max(table_mito$start_ref,na.rm=T)
max(table_mito$end_ref,na.rm=T)

# Creating a exclusive numerical code ("comparable coordinates")
# Exclusive for each NUMT to compare mitocondrial and nuclar genes
# included in same NUMT

table_mito$from_numt <- table_mito$from_numt * 100000
table_mito$start_mtnumt <- table_mito$start_ref + table_mito$from_numt
table_mito$end_mtnumt <- table_mito$end_ref + table_mito$from_numt

mt_genes_coord <- data.frame("mt_hgnc_symbol" = table_mito$hgnc_symbol,
                             "mt_start_numt" = table_mito$start_mtnumt,
                             "mt_end_numt" = table_mito$end_mtnumt)
head(table_mito)
ncol(table_mito)
nrow(table_mito)

# Same process with nuclear genes ---------------

table_numts <- read.table("FINAL_OUTPUT_TABLE.txt", header = TRUE, sep = "\t")

length(duplicated(table_numts)[duplicated(table_numts)==TRUE])
head(table_numts)
str(table_numts)

# Removing duplicated rows  
table_numts <- table_numts[!duplicated(table_numts), ]
table_numts$start_n <- as.numeric(table_numts$start_n)
table_numts$start_position <- as.numeric(table_numts$start_position)
table_numts$end_n <- as.numeric(table_numts$end_n)
table_numts$end_position <- as.numeric(table_numts$end_position)

table_numts$start_ref = NULL
table_numts$end_ref = NULL
for (i in 1:nrow(table_numts)) {
  if (!is.na(table_numts$strand[i])) {
    if (table_numts$strand[1] == 1) {
      table_numts$start_ref[i] <- table_numts$start_position[i] - table_numts$start_n[i] 
      table_numts$end_ref[i] <- table_numts$end_position[i] - table_numts$start_n[i]
    } else {
      table_numts$start_ref[i] <- table_numts$end_n[i] - table_numts$end_position[i]
      table_numts$end_ref[i] <- table_numts$end_n[i] - table_numts$start_position[i]
    }
    i <- i + 1
  } else {
    table_numts$start_ref[i] <- NA
    table_numts$end_ref[i] <- NA
  }
}
table_numts$start_ref <- as.numeric(table_numts$start_ref)
table_numts$end_ref <- as.numeric(table_numts$end_ref)

str(table_numts)
head(table_numts)
table_numts$from_numt <- NULL
table_numts$from_numt <-  substr(table_numts$id, 6,11)
table_numts$from_numt <- gsub("\\.", "", table_numts$from_numt)
table_numts$from_numt <- gsub("\\X", "23", table_numts$from_numt)
table_numts$from_numt <- gsub("\\Y", "24", table_numts$from_numt)
table_numts$from_numt <- as.numeric(table_numts$from_numt)

# In reference interval, replace negative numbers by 0
table_numts$start_ref[table_numts$start_ref<0] <- 0
table_numts$end_ref[table_numts$end_ref<0] <- 0

max(table_numts$start_ref,na.rm=T)
max(table_numts$end_ref,na.rm=T)

table_numts$from_numt <- table_numts$from_numt * 100000
table_numts$start_numt <- table_numts$start_ref + table_numts$from_numt
table_numts$end_numt <- table_numts$end_ref + table_numts$from_numt

n_genes_coord <- data.frame("n_hgnc_symbol" = table_numts$hgnc_symbol,
                            "n_start_numt" = table_numts$start_numt,
                            "n_end_numt" = table_numts$end_numt)
head((mt_genes_coord))
head((n_genes_coord))

c(n_genes_coord[1,2], n_genes_coord[1,3]) %overlaps%
  c(mt_genes_coord[1,2], mt_genes_coord[1,3])

c(n_genes_coord[1,2], n_genes_coord[1,3]) %overlaps%
  c(mt_genes_coord[1,2], mt_genes_coord[1,3])

# Searching for overlaping
i <- 1
b <- 1
overlaping <- NULL
for (i in 1:nrow(n_genes_coord)) for (b in 1:nrow(mt_genes_coord)) 
  if (!is.na(n_genes_coord$n_start_numt[i] & n_genes_coord$n_end_numt[i] & 
             mt_genes_coord$mt_start_numt[b] & mt_genes_coord$mt_end_numt[b])){
    if (c(n_genes_coord$n_start_numt[i], n_genes_coord$n_end_numt[i]) %overlaps%
        c(mt_genes_coord$mt_start_numt[b],
          mt_genes_coord$mt_end_numt[b]))
    {
      overlaping_b <- data.frame("n_hgnc_symbol" = n_genes_coord$n_hgnc_symbol[i],
                                 "mt_hgnc_symbol" = mt_genes_coord$mt_hgnc_symbol[b],
                                 "n_start_numt" = n_genes_coord$n_start_numt[i],
                                 "n_end_numt" = n_genes_coord$n_end_numt[i],
                                 "mt_start_numt" = mt_genes_coord$mt_start_numt[b],
                                 "mt_end_numt" = mt_genes_coord$mt_end_numt[b])
      overlaping <- rbind(overlaping, overlaping_b)
      
    }
  }


overlaping <- overlaping[!duplicated(overlaping), ]
head(overlaping)

overlaping$n_lenght <- overlaping$n_end_numt - overlaping$n_start_numt 
overlaping$mt_lenght <- overlaping$mt_end_numt - overlaping$mt_start_numt 

for (i in 1:nrow(overlaping)) {
  overlaping$ov[i] <- (min(c(overlaping$n_end_numt[i], 
                             overlaping$mt_end_numt[i])) 
                       - max(c(overlaping$n_start_numt[i], 
                               overlaping$mt_start_numt[i])))
}

overlaping$percent_n <- overlaping$ov/overlaping$n_lenght
overlaping$percent_mt <- overlaping$ov/overlaping$mt_lenght
head(overlaping)

# To focuss on nuclear genes that are mainly originated from mitochondrial genes
# we filter the output for at least, 70% of representation of mitochondrial gene
# or 70% of nuclear gene originated from a mitochondrial gene

total_overlaping <- overlaping

# CREATING TABLE total_
write.table(total_overlaping, file = "total_overlaping.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)

overlaping <- subset(overlaping, percent_n >= 0.7 | 
                       percent_mt >= 0.7 )


# ORDERING DATA ----

head(table_mito)
mito <- table_mito[c(1,2,3,4,5,6,8,10,11,12,13,14,16,17)]
head(mito)
colnames(mito)[5] <- "mt_hgnc_symbol"
colnames(mito)[6] <- "mt_ensembl_gene_id"
colnames(mito)[7] <- "mt_gene_biotype"
colnames(mito)[8] <- "mt_start_position"
colnames(mito)[9] <- "mt_end_position"
colnames(mito)[10] <- "mt_strand"
colnames(mito)[11] <- "mt_start_ref"
colnames(mito)[12] <- "mt_end_ref"
colnames(mito)[13] <- "mt_start_numt"
colnames(mito)[14] <- "mt_end_numt"

head(mito)
head(overlaping)

nrow(mito)
nrow(overlaping)
table_overlap <- dplyr::inner_join(mito, overlaping)

head(table_numts)

head(data)
data <- table_numts[c(1:13,72:76,78,79,14:70)]
head(data)

colnames(data)[9] <- "n_hgnc_symbol"
colnames(data)[11] <- "n_gene_biotype"
colnames(data)[12] <- "GO_term"
colnames(data)[14] <- "n_start_position"
colnames(data)[15] <- "n_end_position"
colnames(data)[16] <- "n_strand"
colnames(data)[17] <- "n_start_ref"
colnames(data)[18] <- "n_end_ref"
colnames(data)[19] <- "n_start_numt"
colnames(data)[20] <- "n_end_numt"



data <- data[!duplicated(data), ]
nrow(data)

all_data <- dplyr::full_join(data, table_overlap)
all_data <- all_data[!duplicated(all_data), ]
nrow(all_data)

# Ordering

all_data$n_start_numt <- NULL
all_data$n_end_numt <- NULL
all_data$mt_start_numt <- NULL
all_data$mt_end_numt <- NULL
all_data$NUMT_size <- all_data$end_n  - all_data$start_n 
all_data <- all_data %>% dplyr::select("id", "localization", 
                                       "chr", "start_n", "end_n", 
                                       "NUMT_size",
                                       "mt", "start_mt", "end_mt", 
                                       "Description", "n_gene_biotype", 
                                       "mt_hgnc_symbol","mt_gene_biotype", 
                                       "ensembl_gene_id_version",
                                       
                                       "n_strand",
                                       "n_start_position", "n_end_position",
                                       "mt_ensembl_gene_id",
                                       "mt_strand",
                                       "mt_start_position", "mt_end_position",
                                       "n_start_ref", "n_end_ref",
                                       "mt_start_ref", "mt_end_ref",
                                       "n_lenght", "mt_lenght",
                                       "ov", "percent_n", "percent_mt",
                                       "n_hgnc_symbol", 
                                       everything())

# # # # # # CREATING DEFINITIVE TABLE "all_data_70.txt" # # # # # # #

write.table(all_data, file = "all_data_70.txt", sep = "\t",
            quote = FALSE, row.names = FALSE)


############################################################

# # # # # # MITOCHONDRIAL GENES EXPRESSION # # # # # # #

# Saving complete list of mitochondrial genes ( "mito_genes.txt")

gene_results_mt$ensembl_gene_id
mito_genes <- as.character(na.omit(unique(gene_results_mt$ensembl_gene_id)))
length(mito_genes)
write.table(mito_genes, file="mito_genes.txt",col.names = F,
            sep="\t",quote=F,row.names=F)

# Generating expression data "mean_tpm_GTExMITO.txt"

genes <- read.table("mito_genes.txt", header = FALSE, sep = "\t")
GTEx_mean_tpm <- 
  read.table("GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct", 
             skip = 2, header = TRUE, sep = "\t")

library(plyr); library(dplyr)

GTEx_tpm <- GTEx_mean_tpm
colnames(GTEx_tpm)[1] <- "GTEx_gene_id_version"

gene_id <- GTEx_mean_tpm$gene_id
GTEx_genes <- numeric(0)
for (i in 1:length(GTEx_mean_tpm$gene_id)){
  x <- unlist(strsplit(as.character(GTEx_mean_tpm[i,1]), split='.', fixed=TRUE))[1]
  GTEx_genes <- rbind(GTEx_genes, x)
}

GTEx_mean_tpm$gene_id <- GTEx_genes

mean_tpm_fromGTEx <- numeric(0)
for (i in 1:nrow(genes)){
  y <- subset(GTEx_mean_tpm, gene_id == genes[i,1])
  mean_tpm_fromGTEx <- rbind(mean_tpm_fromGTEx, y)
}

tissue_means <- rowMeans(mean_tpm_fromGTEx[,3:length(mean_tpm_fromGTEx)])
mean_tpm_fromGTEx$tissue_means <- tissue_means

mean_tpm_fromGTEx$sum <- rowSums(mean_tpm_fromGTEx[,3:length(mean_tpm_fromGTEx)])

mean_tpm_fromGTEx$gene_id <- as.character(mean_tpm_fromGTEx$gene_id)

mean_tpm_GTExMITO <- inner_join(mean_tpm_fromGTEx,
                                GTEx_tpm)

mean_tpm_GTExMITO <- mean_tpm_GTExMITO %>% 
  select("gene_id", "GTEx_gene_id_version", everything())

colnames(mean_tpm_GTExMITO)[1] <- "ensembl_gene_id"
str(mean_tpm_GTExMITO)

# Saving results

write.table(mean_tpm_GTExMITO, file = "mean_tpm_GTExMITO.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

mean_tpm_GTExMITO <- read.table("mean_tpm_GTExMITO.txt", header = TRUE, 
                                sep = "\t", dec = ".")
nrow(genes)
ncol(genes)
nrow(mean_tpm_GTExMITO)
ncol(mean_tpm_GTExMITO)

head(mean_tpm_GTExMITO[c(1,2,3)])

GTEx <- mean_tpm_GTExMITO
subset_expressed <- GTEx

# fig.cap= "Heat map of all expressed genes."
library(gplots)
par(oma=c(10,4,4,2))
subset_mean_tpm <- subset_expressed
subset_mean_tpm[1] <- NULL
rownames(subset_mean_tpm) <- subset_mean_tpm$Description
subset_mean_tpm[1] <- NULL
subset_mean_tpm[1] <- NULL

heatmap.2(data.matrix(subset_mean_tpm[1:53]), trace='none', scale = "row",
          cexRow=0.6, cexCol = 0.6)

# fig.cap= "Heat map of all expressed genes."
library(gplots)
par(oma=c(10,4,4,2))
subset_mean_tpm <- subset_expressed
subset_mean_tpm[1] <- NULL
rownames(subset_mean_tpm) <- subset_mean_tpm$Description
subset_mean_tpm[1] <- NULL
subset_mean_tpm[1] <- NULL

heatmap.2(data.matrix(subset_mean_tpm[1:53]),
          trace='none', scale = "column",
          cexRow=0.6, cexCol = 0.6)
