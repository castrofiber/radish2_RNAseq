# Import libs

library(tidyverse)
library(DESeq2)
library(tximport)
library(stringr)
library(magrittr)

# Read inputs

# Read transcript to gene
ttg <- read_tsv(snakemake@input[["t2g"]], col_names = c('target_id', 'gene_id'))

# Import description table for radish transcripts
description_table <- read_tsv(snakemake@input[["description"]])

# import Arabidopsis BLAST results
blast_results <- read_tsv(snakemake@input[["annotation"]], col_names = F) %>%
  distinct(X1, .keep_all = T) %>%
  transmute(target_id = X1, transcript_id = X2, locus_id = str_extract(transcript_id, "AT[0-9]G[0-9]{5}"))

# combine tables to create annotation table
annotation_table <- left_join(ttg, blast_results, by="target_id") %>%
  left_join(description_table, by="target_id") %>%
  mutate(annotation = ifelse(is.na(transcript_id), gene_id, locus_id))

# Read phenotypes
phenotable <- read_tsv(snakemake@input[["phenotable"]])

# Read abundance from kallisto
files <- snakemake@input[["abundance"]]

txi <- tximport(files, type = 'kallisto', tx2gene = ttg)

# design is hardcoded and "condition" column existence in phenotable is implied

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = phenotable,
                                   design = ~ condition)

# Filter genes by row counts

dds <- ddsTxi[rowSums(counts(ddsTxi)) >= 10,]

# Relevel condition factor

dds$condition <- relevel(dds$condition, ref = "K")

# Running DE

dds <- DESeq(dds)

resLFC <- lfcShrink(dds, coef="condition_O_vs_K", type="apeglm", lfcThreshold = 1) %>%
  as.data.frame %>%
  rownames_to_column("gene") %>%
  arrange(svalue) %>%
  mutate(stat = dds %>% results %>% as.data.frame %>% pull(stat)) %>%  # pull Wald test stat for ranking
  left_join(annotation_table, c("gene" = "gene_id")) %>%  # plug in annotation
  transmute(gene, target_id, annotation, log2FoldChange, svalue, description, stat) %>%
  distinct(gene, .keep_all = T)

norm_expr <- vst(dds, blind=FALSE) %>%
  assay %>%
  set_colnames(phenotable$sample) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

sign_up <- resLFC %>%
  filter(svalue > 0.05, log2FoldChange < 0)

sign_down <- resLFC %>%
  filter(svalue < 0.05, log2FoldChange < 0)

## Return outputs

write_csv(resLFC, snakemake@output[["DEtable"]])

write_csv(norm_expr, snakemake@output[["normalized_expression"]])

write_csv(sign_up, snakemake@output[["up_genes"]])

write_csv(sign_down, snakemake@output[["down_genes"]])
