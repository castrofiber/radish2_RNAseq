# Import libs

library(biomaRt)
library(tidyverse)
library(clusterProfiler)
library(fgsea)
library(org.At.tair.db)
library(stringr)
library(cowplot)
# library(DOSE)

# Read DE results and arrange by significance
DEtable <- read_csv(snakemake@input[["DEtable"]]) %>% 
  arrange(svalue)

# GSEA

ranks <- DEtable %>% 
  distinct(annotation, .keep_all = T)

ranks_vector <- setNames(ranks$log2FoldChange, ranks$annotation)

## Custom function to read plantGSEA
custom_gmt <- function(gmt.file, category, include = T)
{
  lines <- str_split(readLines(gmt.file), "\t")
  if (include)
  {
    lines <- Filter(function(x) x[2] == category, lines)
  } else
    lines <- Filter(function(x) x[2] != category, lines)
  pathways <- lapply(lines, function(x) tail(x, 1) %>% str_split(",") %>% unlist)
  names(pathways) <- sapply(lines, head, 1)
  pathways
}

gmt <- custom_gmt(snakemake@input[["PlantGSEA"]], category = "Literature", include = F)

gsea_res <- fgsea(pathways = gmt,
                  stats = ranks_vector,
                  minSize = 15,
                  maxSize = 500,
                  eps = 0)

gsea_sign <- gsea_res %>% 
  filter(padj < 0.05)

gsea_sign_up <- gsea_sign %>% 
  filter(NES > 0) %>% 
  arrange(padj)

gsea_sign_down <- gsea_sign %>% 
  filter(NES < 0) %>% 
  arrange(padj)

gsea_sign_up$leadingEdge <- gsea_sign_up$leadingEdge %>% 
  lapply(function(x) x %>% paste(collapse = "/")) %>% unlist
gsea_sign_up %>% write_tsv(snakemake@output[["GSEA_up"]])

gsea_sign_down$leadingEdge <- gsea_sign_down$leadingEdge %>% 
  lapply(function(x) x %>% paste(collapse = "/")) %>% unlist
gsea_sign_down %>% write_tsv(snakemake@output[["GSEA_down"]])


## Custom function to make a volcano plot but with GSEA results
plot_gsea_volcano <- function(gsea_results)
{
  # gsea_results - data frame from fgsea function output
  gsea_results <- gsea_results %>% 
    mutate(sign = padj < 0.05) %>% 
    arrange(NES) %>% 
    mutate(top = ifelse(row_number() < 2 | row_number() > (n() - 1), T, F))
  
  ggplot(gsea_results, aes(NES, -log10(padj)))+
    geom_point(aes(color = sign))+
    scale_color_manual(values = c("black", "red"))+
    geom_hline(yintercept = -log10(0.05), color = "red")+
    geom_text(aes(label = ifelse(top, pathway, "")))+
    theme_bw()
}

GSEA_volcano <- plot_gsea_volcano(gsea_res) + ylim(0, 5)

ggsave(snakemake@output[["GSEA_volcano"]], plot=GSEA_volcano, height=25, width=20)

### GO and KEGG enrichment with clusterProfiler

mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "athaliana_eg_gene",
                         host = 'plants.ensembl.org')

mart_annotation <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "tair_locus",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "entrezgene_id"),
  mart = mart) %>% 
  dplyr::select(c('tair_locus', 'entrezgene_id'))

top <- ranks %>% 
  left_join(mart_annotation, by = c("annotation" = "tair_locus")) %>% 
  dplyr::filter(svalue < 0.05)

top_up <- top %>% 
  dplyr::filter(log2FoldChange > 0)

top_down <- top %>% 
  dplyr::filter(log2FoldChange < 0)

## Remove NA from enrich_results
remove_NA_from_enrichment <- function(enrichment_object)
{
  enrichment_object@result$geneID <- enrichment_object@result$geneID %>% str_replace_all("NA/|NA$", "")
  enrichment_object
}

GO_down <- enrichGO(gene = top_down$entrezgene_id,
                    pAdjustMethod = "BH",
                    OrgDb = org.At.tair.db,
                    ont = "ALL",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05,
                    readable      = TRUE) %>% 
  remove_NA_from_enrichment()

GO_up <- enrichGO(gene = top_up$entrezgene_id,
                  pAdjustMethod = "BH",
                  OrgDb = org.At.tair.db,
                  ont = "ALL",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE) %>% 
  remove_NA_from_enrichment()

## Inspect results

GO_down@result  %>% write_tsv(snakemake@output[["GO_down"]])
GO_up@result %>% write_tsv(snakemake@output[["GO_up"]])

GO_dotplot <- plot_grid(
  dotplot(GO_up, showCategory=15) + ggtitle("GO_up"),
  dotplot(GO_down, showCategory=15) + ggtitle("GO_down"),
  ncol = 2
)

save_plot(snakemake@output[["GO_dotplot"]], GO_dotplot, base_height=20, base_width=25)

# KEGG enrichment

KEGG_down <- enrichKEGG(gene = top_down$annotation,
                        organism     = 'ath',
                        pvalueCutoff = 0.05) %>%
  setReadable('org.At.tair.db', keyType = "TAIR") %>% 
  remove_NA_from_enrichment() 

KEGG_up <- enrichKEGG(gene = top_up$annotation,
                      organism     = 'ath',
                      pvalueCutoff = 0.05) %>%
  setReadable('org.At.tair.db', keyType = "TAIR") %>% 
  remove_NA_from_enrichment() 


KEGG_dotplot <- plot_grid(
  dotplot(KEGG_up, showCategory=15) + ggtitle("KEGG_up"),
  dotplot(KEGG_down, showCategory=15) + ggtitle("KEGG_down"),
  ncol = 2
)

save_plot(snakemake@output[["KEGG_dotplot"]], KEGG_dotplot, base_height=20, base_width=25)



