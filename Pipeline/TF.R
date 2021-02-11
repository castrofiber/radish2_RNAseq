# Import libs

library(tidyverse)
library(pheatmap)
library(ggplot2)

## Read and modify inputs

# Read DE results and filter by significance
sign_results <- read_csv(snakemake@input[["DEtable"]]) %>% 
  filter(svalue < 0.05)

# Read phenotypes
phenotable <- read_tsv(snakemake@input[["phenotable"]])

# Read expression matrix with normalized values
normalized_expression <- read_csv(snakemake@input[["normalized_expression"]])

# TF annotation
TF_annotation_table <- read_tsv(snakemake@input[["AtTFs"]], col_names = F) %>% 
  transmute(TF_loci = toupper(X2))

sign_TFs <- left_join(TF_annotation_table, sign_results, by = c("TF_loci" = "annotation")) %>% 
  drop_na() %>% 
  transmute(gene)

# There are problems with this step: some nonTFs are homologous to TFs and seem like TFs in results: example - AT5G55970
# Manually filtered results can be plugged in via snakemake
# sign_TFs <- read_tsv(snakemake@input[["manual_TFs"]], col_names = T) %>% 
#   dplyr::filter(`is TF?` == 1) %>% 
#   transmute(gene = `Locus ID`)

### Drawing a pheatmap
col_annotation <- phenotable %>% 
  mutate(condition = case_when(condition == 'K' ~ "Root", TRUE ~ "Tumor")) %>% 
  column_to_rownames("sample")

row_annotation <- left_join(TF_annotation_table, sign_results, by = c("TF_loci" = "annotation")) %>% 
  drop_na() %>% 
  transmute(gene, description) %>% 
  column_to_rownames("gene")

plot <- normalized_expression %>% 
  left_join(sign_TFs, ., by = "gene") %>%
  column_to_rownames("gene") %>%
  pheatmap(cluster_rows = T,
           cluster_cols = T,
           annotation_col = col_annotation,
           labels_row = row_annotation$description,
           main = "Differentially expressed TFs")

ggsave(snakemake@output[["figure"]], plot=plot, height=25, width=20)
