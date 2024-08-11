R.version

library(tximport)
packageVersion("tximport")

setwd("/home/RNAseq_tumors/RSEM_align_quant/")

# Metadata file
samples <- read.csv("metadata_samples.csv", header = T)
samples <- samples[1:18,]
head(samples)
dim(samples)

# Make list of files for import. For this analysis, .genes.results files from RSEM alignment are in </home/RNAseq_tumors/RSEM_align_quant/RSEM_quant_results/>. These files can be downloaded from GitHub.
files <- file.path(getwd(), "RSEM_quant_results", list.files("RSEM_quant_results/"))
files

# Name each item in list. These names are used for tximport sample naming, not the full file names.
names(files) <- samples$Sample

# Import as RSEM files. Set txIn = F and txOut = F for gene results. (Isoform results are txIn = T and txOut = T)
txi.rsem <- tximport(files, type = "rsem", txIn = F, txOut = F)

### txi.rsem object contains 3 tables (txi.rsem$counts, txi.rsem$length, txi.rsem$abundance)

## Some genes have "length = 0" in the txi.rsem object. Remove these rows from all tables in txi.rsem

# Determine number of genes before removing any
dim(txi.rsem$counts)

# For each row, determine if the length value is zero
row_sub <- apply(txi.rsem$length, 1, function(row) all(row !=0 ))

# Subset each table in txi.rsem for these rows
txi.rsem$length <- txi.rsem$length[row_sub,]
txi.rsem$abundance <- txi.rsem$abundance[row_sub,]
txi.rsem$counts <- txi.rsem$counts[row_sub,]

# Determine number of genes after removing those with length=0
dim(txi.rsem$counts)

# 55416-51078 = removed ~4000 genes

library(DESeq2)
packageVersion("DESeq2")

### Make sampleTable from metadata table, with some reformatting
sampleTable <- data.frame(Genotype = samples$Genotype, Line = samples$Line, Sex = samples$Sex)
rownames(sampleTable) <- colnames(txi.rsem$counts)

head(sampleTable)

### Import data from txi.rsem into DESeq2
dds_geno <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~ Line + Sex + Genotype)

# The variable of interest is listed LAST in the formula. The variables before Genotype will be controlled for

# Relevel the Genotype factors so that WT is the first/reference level
dds_geno$Genotype <- relevel(dds_geno$Genotype, "WT")

levels(dds_geno$Genotype)

# Filter out genes with 0 or 1 counts across all samples
nrow(dds_geno)

keep <- rowSums(counts(dds_geno)) > 1
dds_geno <- dds_geno[keep,]

nrow(dds_geno)

# DESeq function runs differential expression pipeline on raw counts from dds. Estimates size factors (to control for differences in sequencing depth across samples), estimates dispersion values for each gene, and fits a generalized linear model.
dds_geno_DE <- DESeq(dds_geno)

### Generate results table
# contrast syntax = c("condition", "level_to_compare", "base_level")

contrast_geno <- c("Genotype", "H3KO", "WT")
results_geno <- results(dds_geno_DE, contrast = contrast_geno, alpha = 0.05)

summary(results_geno)

# Select for up-/down-regulated genes with p-adjusted cutoff of 0.05
res.Sig_geno <- subset(results_geno, padj < 0.05)

# Export file with all genes and file with p-adjusted < 0.05. Order by log2-fold change magnitude.

resOrdered <- results_geno[order(results_geno$log2FoldChange), ]
write.csv(resOrdered, file = "DESeq2_all.csv")

resOrdered.Sig <- res.Sig_geno[order(res.Sig_geno$log2FoldChange), ]
write.csv(resOrdered.Sig, file = "DESeq2_padj0.05.csv")
