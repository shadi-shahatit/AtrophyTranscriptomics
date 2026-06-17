#### Transcriptomics Analysis Training - Al-Zghoul Lab, JUST
#### Shadi Shahatit - RA, JUST 2026
# libraries ---------------------------------------------------------------

## Install required packages
## CRAN packages are installed via install.packages(); Bioconductor packages via BiocManager::install()
## ONLY instal ONCE

pkgs <- c(
  "tidyverse", "forcats", "stringr", "tibble", "readxl", "writexl",
  "tximport", "txdbmaker", "AnnotationDbi", "org.Gg.eg.db", "DESeq2",
  "apeglm", "clusterProfiler", "enrichplot", "topGO", "GOSemSim",
  "simplifyEnrichment", "ggrepel", "ggsci", "viridis", "viridisLite",
  "RColorBrewer", "ggnewscale", "patchwork", "gridExtra", "UpSetR",
  "ggvenn", "VennDiagram")

bioc_pkgs <- c(
  "tximport", "txdbmaker", "AnnotationDbi", "org.Gg.eg.db", "DESeq2",
  "apeglm", "clusterProfiler", "enrichplot", "topGO", "GOSemSim",
  "simplifyEnrichment")

cran_pkgs <- setdiff(pkgs, bioc_pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(cran_pkgs[!cran_pkgs %in% installed.packages()[,"Package"]])
BiocManager::install(bioc_pkgs[!bioc_pkgs %in% installed.packages()[,"Package"]])

## Load all packages at the start of each session using the library()

library(tidyverse)
library(forcats)
library(stringr)
library(tibble)
library(readxl)
library(writexl)
library(tximport)
library(txdbmaker)
library(AnnotationDbi)
library(org.Gg.eg.db)
library(DESeq2)
library(apeglm)
library(clusterProfiler)
library(enrichplot)
library(topGO)
library(GOSemSim)
library(simplifyEnrichment)
library(ggrepel)
library(ggsci)
library(viridis)
library(viridisLite)
library(RColorBrewer)
library(ggnewscale)
library(patchwork)
library(gridExtra)
library(UpSetR)
library(ggvenn)
library(VennDiagram)

sys_dir <- "/home/shadi/Desktop/transcriptome_alzghoul_SS/TMBroilers_Transcriptomics/"

# Gene-level quantification - tximport -----------------------------------

## tximport

## tximport to get gene-level quantification from Salmon transcript-level output

## make sure you have installed a reference GTF file for your organism in the correct directory (Ggallus_Refs/gtf_r115)
## in bash, use wget https://ftp.ensembl.org/pub/release-115/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.chr.gtf.gz        	# ref anno - gtf release 115
## convert GTF to TxDb

# TxDb object from Ensembl release 115 GTF for use in genomic feature annotation and transcript-level analyses

## A GTF (Gene Transfer Format) file includes genomic features (e.g., genes, transcripts, exons, CDS) along with their chromosomal coordinates and metadata
## A TxDb object is a Bioconductor database that indexes these features for efficient use in R

gtf_dir <- file.path(sys_dir, "transcriptomics_training/Ggallus_Refs/gtf_r115/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.115.chr.gtf.gz")
txdb_Ggallus <- txdbmaker::makeTxDbFromGFF(gtf_dir,
                                           format="gtf",
                                           dataSource="Ensembl",
                                           organism="Gallus_gallus",
                                           taxonomyId=9031)

## extract transcript-to-gene mapping from TxDb object
## TXNAME (transcript ID) is mapped to its parent gene ID (GENEID)

k <- keys(txdb_Ggallus, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb_Ggallus, k, "GENEID", "TXNAME")

length(unique(tx2gene$GENEID))
length(unique(tx2gene$TXNAME))

## upload salmon for tximport and sum into gene level quants

## your quant files should be your last results from the first stage
## To run an example
## download quant files for 4 conditions in chicken liver D22 (3 samples each) from https://github.com/shadi-shahatit/AtrophyTranscriptomics/tree/main/scripts_training/quant_r115

quants_dir <- file.path(sys_dir, "transcriptomics_training/quant_r115")
files <- list.files(quants_dir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
names(files) <- basename(dirname(files))
TxiSalmon_liver <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)

head(TxiSalmon_liver$counts)     # raw_counts
head(TxiSalmon_liver$abundance)  # TPM
head(TxiSalmon_liver$length)     # gene_lengths (bps)
ncol(TxiSalmon_liver$counts)     # 12 is the number of our samples in this example 

## get exp matrix at the level of transcripts as well

TxiTranscripts_liver <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE,
                                 txOut = TRUE)
transcript_counts_liver <- as.data.frame(TxiTranscripts_liver$counts)

## plot the distribution of transcript counts per gene.

# to remove transcript versions
transcript_counts_liver$TXNAME <- gsub("\\.\\d+$", "", rownames(transcript_counts_liver))
transcript_counts_liver_mapped <- merge(transcript_counts_liver, tx2gene, by = "TXNAME")
transcripts_per_gene_liver <- transcript_counts_liver_mapped %>%
  dplyr::select(TXNAME, GENEID) %>%
  distinct() %>%
  group_by(GENEID) %>%
  summarise(transcript_count = n())

ggplot(subset(transcripts_per_gene_liver,transcript_count<100), aes(x = transcript_count)) +
  geom_histogram(binwidth = 1, fill = "#323232") +
  labs(title = "", x = "Transcripts per gene", y = "Gene counts") +
  theme_classic()

summary(transcripts_per_gene_liver$transcript_count)
mean(transcripts_per_gene_liver$transcript_count)
sd(transcripts_per_gene_liver$transcript_count)
nrow(transcripts_per_gene_liver)

# Metadata and ExpMat preparation  ---------------------------------------------------

## prep your metadata - liver

## create a table with consistent naming with txi salmon files

## VERY specific to each example (could also be loaded from a file)

colnames(TxiSalmon_liver$counts)

metadata_liver_D22 <- tibble(sample_id = colnames(TxiSalmon_liver$counts)) %>%
  mutate(id = sample_id,
         lower = tolower(sample_id),
         tissue = "liver",
         day = case_when(
           str_detect(lower, "d22|a22|at22|ct22") ~ "D22",
           TRUE ~ NA_character_),
         condition = case_when(
           str_detect(lower, "^a") & str_detect(lower, "(^|_)con(_|$)") ~ "Con_AHS",
           str_detect(lower, "^a") & str_detect(lower, "(^|_)tm(_|$)") ~ "TM_AHS",
           ! str_detect(lower, "^a") & str_detect(lower, "(^|_)con(_|$)") ~ "Con_N",
           ! str_detect(lower, "^a") & str_detect(lower, "(^|_)tm(_|$)") ~ "TM_N",
           TRUE ~ NA_character_)) %>%
  transmute(id,
            sample_id,
            tissue,
            day = factor(day, levels = c("D22")),
            condition = factor(condition, levels = c("Con_N", "TM_N", "Con_AHS", "TM_AHS"))) %>%
  column_to_rownames("sample_id")

stopifnot(identical(rownames(metadata_liver_D22), colnames(TxiSalmon_liver$counts)))

unique(metadata_liver_D22$condition)

## The analysis of Con_N, TM_N, Con_AHS, and TM_AHS in Day 22
## Con_N:     thermal neutral control
## TM-N:      thermally manipulated under neutral conditions
## Con-AHS:   control subjected to acute heat stress
## TM-AHS:    thermally manipulated subjected to acute heat stress

TxiSalmon_liver_D22 <- TxiSalmon_liver

## QC

## QC - overall distribution of counts

## plot the data

## when plotting expression levels across conditions one must use normalized counts for scalability 

TxiSalmon_liver_TPM_tidy <- as.data.frame(TxiSalmon_liver_D22$abundance) %>% ## TPM
  rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to="id", values_to="TPM") %>%
  left_join(metadata_liver_D22, by="id")

condition_colors <- c(Con_N = "#D62728",
                      TM_N = "#1F77B4",
                      Con_AHS="#2CA02C",
                      TM_AHS="#9467BD"
)

ggplot(TxiSalmon_liver_TPM_tidy, aes(x=reorder(id, as.numeric(condition)), y=log10(TPM+1), fill=condition))+
  geom_boxplot(outlier.size=0.5)+
  scale_fill_manual(values=condition_colors, name="Conditions")+
  labs(x="Sample",y="log10 TPM")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggplot(TxiSalmon_liver_TPM_tidy,aes(x=condition,y=log10(TPM+1),fill=condition))+
  geom_boxplot(outlier.size=0.5)+
  scale_fill_manual(values=condition_colors, name="Condition")+
  labs(x="Condition",y="log10 TPM")+
  theme_classic()

## QC - PCA

## get a pca from DESeqDataSet from Salmon counts, filter lowly expressed genes, apply VST for variance stabilization,
## and extract PCA coordinates with variance explained by PC1 and PC2

dds_liver_D22_pca_00 <- DESeqDataSetFromTximport(
  txi = TxiSalmon_liver_D22,
  colData = metadata_liver_D22,
  design = ~ condition) # model does not matter in PCA as blind is True

dds_liver_D22_pca_01 <- dds_liver_D22_pca_00[rowSums(counts(dds_liver_D22_pca_00)) >= 10, ]
vsd_liver_D22_all <- vst(dds_liver_D22_pca_01, blind = TRUE)
pca_liver_D22_df  <- plotPCA(vsd_liver_D22_all, intgroup = c("condition","day"), returnData = TRUE, 
                             ntop = Inf) # top 500, 2000, or Inf genes; we need info from all genes
pv_liver_D22 <- round(100 * attr(pca_liver_D22_df, "percentVar")[1:2], 1)

ggplot(pca_liver_D22_df, aes(x = PC1, y = PC2, color = condition, shape = day)) +
  geom_point(size = 3) +
  scale_color_manual(values = condition_colors) +
  scale_shape_discrete(drop = FALSE) +
  coord_equal() +
  labs(x = paste0("PC1 (", pv_liver_D22[1], "%)"),
       y = paste0("PC2 (", pv_liver_D22[2], "%)"),
       color = "Condition",
       shape = "Day") +
  theme_bw()

# DESeq2 modeling -----------------------------------------------------------

## you need the following two var to compare conditions  
TxiSalmon_liver_D22
metadata_liver_D22

## several stats model can be considered:
## (1) Per condition pairwise comp      - Wald

## create a DESeq object
dds_liver_D22_00 <- DESeqDataSetFromTximport(
  txi = TxiSalmon_liver_D22,
  colData = metadata_liver_D22,
  design = ~ condition) # Pairwise comp of all condition at D22; Wald

## Pre-filter low-expression genes
dds_liver_D22_01 <- dds_liver_D22_00[rowSums(counts(dds_liver_D22_00)) >= 10, ]
# dds$day <- droplevels(dds$day)
# dds$condition <- droplevels(dds$condition)

## DESeq modeling
dds_liver_D22 <- DESeq(dds_liver_D22_01)

# estimating size factors
# using 'avgTxLength' from assays(dds), correcting for library size
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(counts(dds_liver_D22))
colData(dds_liver_D22)
summary(dds_liver_D22)
levels(dds_liver_D22$condition)

plotDispEsts(dds_liver_D22)

# This is a DESeq2 dispersion (degree of variability) estimate plot used to quality check the fit of the model (negative binomial)
# The red line is the mean dispersion trend fitted across all genes,
# The classic funnel shape with high dispersion at low counts converging toward stability at high counts is expected
# and confirms the model is fitting the data appropriately

## extract results (AEG_full_info_liver_D22; DEG_full_info_liver_D22)

# AEG_full_info_liver_D22: All Expressed Genes
# DEG_full_info_liver_D22: Differentially Expressed Genes

## main function is results(dds_liver_D22, contrast = c("condition","TM_N","Con_N"))

padj_thresh <- 0.05
lfc_thresh  <- 1

pairs_target <- list(
  c("TM_N","Con_N"),
  c("Con_AHS","Con_N"),
  c("Con_AHS","TM_N"),
  c("TM_AHS","Con_N"),
  c("TM_AHS","TM_N"),
  c("TM_AHS","Con_AHS")
)

AEG_full_info_liver_D22_list <- lapply(seq_along(pairs_target), function(i){
  g1 <- pairs_target[[i]][1]; g2 <- pairs_target[[i]][2]
  labs <- paste(g1, "vs", g2)
  res <- try(results(dds_liver_D22, contrast = c("condition", g1, g2)), silent = TRUE)
  tibble::as_tibble(res, rownames = "Gene") %>% mutate(contrast = labs)
})

AEG_full_info_liver_D22 <- dplyr::bind_rows(AEG_full_info_liver_D22_list)

AEG_full_info_liver_D22$contrast <- factor(AEG_full_info_liver_D22$contrast,
                                           levels = vapply(pairs_target, function(p) paste(p[1],"vs",p[2]), ""))

DEG_full_info_liver_D22 <- AEG_full_info_liver_D22 %>%
  filter(!is.na(padj),
         padj < padj_thresh,
         abs(log2FoldChange) > lfc_thresh) %>%
  mutate(direction = if_else(log2FoldChange > 0, "upregulated", "downregulated"))

length(unique(AEG_full_info_liver_D22$Gene)) < length(unique(rownames(TxiSalmon_liver_D22$counts)))
nrow(DEG_full_info_liver_D22)

## append different gene ids and download your final list of DEGs

## OrgDb are the annotation packages for particular organisms
## using org.Gg.eg.db, we map Ensembl IDs to gene symbols, Entrez IDs, and RefSeq accessions,
## collapsing many-to-one mapping with semicolons, and then add annotations back to the DEGs table

collapse_unique <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

IDmap_Gg <- AnnotationDbi::select(
  org.Gg.eg.db,
  keys    = unique(DEG_full_info_liver_D22$Gene),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL","SYMBOL","ENTREZID","REFSEQ")) %>%
  as_tibble() %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL   = collapse_unique(SYMBOL),
    ENTREZID = collapse_unique(ENTREZID),
    REFSEQ   = collapse_unique(REFSEQ),
    .groups  = "drop")

DEG_full_info_liver_D22_all_anno <- DEG_full_info_liver_D22 %>%
  left_join(IDmap_Gg, by = c(Gene = "ENSEMBL")) %>%
  mutate(SYMBOL = if_else(is.na(SYMBOL) | SYMBOL == "", Gene, SYMBOL))

write_xlsx(
  DEG_full_info_liver_D22_all_anno,
  file.path(sys_dir, "transcriptomics_training", "TMomics_DEG_liver_D22_allanno_v1.xlsx"))

## lfcShrink methods

# res_1 <- results(dds_liver_D22, contrast = c("condition","TM_N","Con_N"))
# res_1_apeglm <- lfcShrink(dds_liver_D22, coef = "condition_TM_N_vs_Con_N", type = "apeglm")
# DESeq2::plotMA(res_1)
# DESeq2::plotMA(res_1_apeglm)

# DEGs visualizations (Bar + Volcano + Venn diagrams + Upset plots) -----------------------------------------------------------

GeneSet_PairSet_target <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    "Con_AHS vs TM_N",
    "TM_AHS vs Con_N",
    "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"
  )
)

## Bar plots

## count and pivot DEGs per contrast

DEG_liver_D22_sum <- DEG_full_info_liver_D22 %>%
  group_by(contrast, direction) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
  mutate(total = upregulated + downregulated)

DEG_liver_D22_sum_long <- DEG_liver_D22_sum %>%
  pivot_longer(cols = c("upregulated", "downregulated", "total"),
               names_to = "category", values_to = "count") %>% 
  mutate(category = factor(category, levels = c("downregulated", "upregulated", "total")))

ggplot(subset(DEG_liver_D22_sum_long), aes(x = contrast, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  theme_classic() +
  labs(title = "",
       x = "Conditions",
       y = "DEGs counts",
       fill = "DEGs category") +
  scale_fill_manual(values = c("#FFD662FF", "#00539CFF", "#949398FF")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

## Volcano plots - v2

collapse_unique <- function(x) {
  x <- unique(na.omit(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

IDmap_Gg <- AnnotationDbi::select(
  org.Gg.eg.db,
  keys    = unique(DEG_full_info_liver_D22$Gene),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL","SYMBOL")) %>%
  as_tibble() %>%
  group_by(ENSEMBL) %>%
  summarise(
    SYMBOL   = collapse_unique(SYMBOL),
    .groups  = "drop")

DEG_full_info_liver_D22_all_anno <- DEG_full_info_liver_D22 %>%
  left_join(IDmap_Gg, by = c(Gene = "ENSEMBL")) %>%
  mutate(SYMBOL = if_else(is.na(SYMBOL) | SYMBOL == "" | grepl(";", SYMBOL), NA, SYMBOL))

volc_cols <- c("Upregulated" = "#22A884FF", "Downregulated" = "#440154FF", "Not Significant" = "#B3B3B3")

VolcanoPlot_label <- function(df, df_DEG, ttl){
  df2 <- df %>% mutate(direction = case_when(
    log2FoldChange >  lfc_thresh & padj < padj_thresh ~ "Upregulated",
    log2FoldChange < -lfc_thresh & padj < padj_thresh ~ "Downregulated",
    TRUE ~ "Not Significant")) %>%
    left_join(df_DEG %>% select(Gene, contrast, SYMBOL), by = c("Gene", "contrast"))
  
  # df_labels <- df2 %>% filter(direction != "Not Significant")
  df_labels <- df2 %>% 
    filter(direction != "Not Significant") %>%
    group_by(contrast, direction) %>%
    slice_max(order_by = log2FoldChange, n = 10) %>%
    ungroup()
  
  upc <- sum(df2$direction=="Upregulated", na.rm=TRUE)
  dnc <- sum(df2$direction=="Downregulated", na.rm=TRUE)
  ymax <- max(-log10(df2$padj), na.rm=TRUE); if(!is.finite(ymax)) ymax <- 1
  
  ggplot(df2, aes(log2FoldChange, -log10(padj), color=direction)) +
    geom_point(size=1.6, alpha=0.8) +
    scale_color_manual(values=volc_cols) +
    geom_hline(yintercept=-log10(padj_thresh), linetype="dashed") +
    geom_vline(xintercept=c(-lfc_thresh, lfc_thresh), linetype="dashed") +
    annotate("text", x=-5.5, y=ymax, label=dnc, size=3.8, fontface="bold") +
    annotate("text", x= 5.5, y=ymax, label=upc, size=3.8, fontface="bold") +
    xlim(-6,6) + labs(title=ttl, x="log2 FC", y="-log10(q)") +
    geom_text_repel(data = df_labels, aes(label = SYMBOL, color = direction),
                    size = 3, max.overlaps = Inf,
                    box.padding = 0.5, point.padding = 0.25,
                    show.legend = FALSE) +
    theme_minimal(base_size=11) +
    theme(panel.grid=element_blank(),
          panel.border=element_rect(color="black", fill=NA, linewidth=0.6),
          legend.position="none",
          plot.title=element_text(hjust=0.5, face="bold"))
}

AEG_full_info_liver_D22_set1 <- AEG_full_info_liver_D22 %>%
  filter(contrast %in% GeneSet_PairSet_target$set1) %>%
  droplevels()

VolcPlot_liver_D22_set1 <- lapply(
  split(AEG_full_info_liver_D22_set1, AEG_full_info_liver_D22_set1$contrast),
  function(d) VolcanoPlot_label(d %>% select(Gene, log2FoldChange, padj, contrast),
                                DEG_full_info_liver_D22_all_anno,
                                ttl = unique(d$contrast)[1]))

wrap_plots(VolcPlot_liver_D22_set1, ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 20, face = "bold", family = "Times New Roman",color = "white"))

## Upset plots + Venn diagrams

GeneSet_PairSet_target <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    "Con_AHS vs TM_N",
    # "TM_AHS vs Con_N",
    # "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"
  )
)

GeneSet_liver_D22 <- DEG_full_info_liver_D22_all_anno %>%
  group_by(contrast) %>%
  summarize(Gene_Symbols = list(unique(Gene)), .groups = "drop") %>%
  deframe()

GeneSet_liver_D22_1 <- GeneSet_liver_D22[GeneSet_PairSet_target$set1]

## Upset

GeneSet_Matrix_liver_D22_1 <- fromList(GeneSet_liver_D22_1)

UpSetR::upset(GeneSet_Matrix_liver_D22_1,
              main.bar.color = "#404080",
              sets.bar.color = "#69b3a2",
              order.by = "freq",
              sets = names(GeneSet_liver_D22_1),
              set_size.show = TRUE)

## Venn diagrams

ggvenn(
  GeneSet_liver_D22_1, 
  fill_color = viridis(length(names(GeneSet_liver_D22_1))),
  stroke_size = 0.5, 
  set_name_size = 4, 
  text_size = 3)

# Functional enrichment ---------------------------------------------------------------------

## using these dfs
## AEG_full_info_liver_D22
## DEG_full_info_liver_D22
## dds_liver_D22

## using helper functions to perform analysis and plot
## perform_ORA_GO_BP
## ORA_PublishFEBarplot

## over-representation analysis ORA

## helper function for ordered panel FC bar plots: ORA table (Description, Fold Enrichment, p-value, contrast)
ORA_PublishFEBarplot <- function(ora_df, n_top = 10, ncol = 3, title_col = "contrast", order_levels = NULL) {
  dfp <- ora_df %>%
    mutate(
      log10p = -log10(p.adjust),
      Desc   = str_wrap(Description, width = 40),
      !!title_col := if (is.null(order_levels)) .data[[title_col]]
      else factor(.data[[title_col]], levels = order_levels)
    ) %>%
    group_by(.data[[title_col]]) %>%
    slice_min(order_by = p.adjust, n = n_top, with_ties = FALSE) %>%
    ungroup()
  
  levs <- if (is.null(order_levels)) levels(factor(dfp[[title_col]])) else order_levels
  plots <- lapply(levs, function(lbl) {
    dd <- dfp %>% filter(.data[[title_col]] == lbl)
    ggplot(dd, aes(x = FoldEnrichment, y = fct_reorder(Desc, FoldEnrichment), fill = log10p)) +
      geom_bar(stat = "identity", width = 0.7) +
      scale_fill_gradient(low = "firebrick", high = "darkgreen", name = expression(-log[10](q))) +
      labs(x = "Fold enrichment", y = NULL, title = lbl) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.y        = element_text(size = 11, color = "black"),
        axis.text.x        = element_text(size = 11, color = "black"),
        axis.title.x       = element_text(size = 13, color = "black"),
        plot.title         = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title       = element_text(size = 13, color = "black"),
        legend.text        = element_text(size = 12, color = "black"),
        panel.grid.major.y = element_blank()
      )
  })
  wrap_plots(plots, ncol = ncol)
}

## make sure to use adj. p-values

## DEGs list per contrast with ENSEMBL IDs and a specific pairwise sets

degs_sets_liver_D22 <-
  DEG_full_info_liver_D22 %>%
  distinct(contrast, Gene) %>%
  group_by(contrast) %>%
  summarise(genes = list(unique(Gene)), .groups = "drop") %>%
  { setNames(.$genes, .$contrast) }

## in case of targeted pairwise comp

GeneSet_PairSet_target <- list(
  set1 = c(
    "TM_N vs Con_N",
    "Con_AHS vs Con_N",
    "Con_AHS vs TM_N",
    "TM_AHS vs Con_N",
    "TM_AHS vs TM_N",
    "TM_AHS vs Con_AHS"
  )
)

## Background universe genes with ENSEMBL IDs

background_genes_liver_D22 <- rownames(dds_liver_D22)

## over-representation analysis ORA (GO:BP)

perform_ORA_GO_BP <- function(gene_set, label) {
  enrichGO(
    gene          = gene_set,
    OrgDb         = org.Gg.eg.db,
    keyType       = "ENSEMBL",
    ont           = "BP",
    universe      = background_genes_liver_D22,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = FALSE
  ) |>
    as.data.frame() |>
    dplyr::mutate(contrast = label)
}

ORA_GO_liver_D22_list <- lapply(names(degs_sets_liver_D22), function(lbl) {
  perform_ORA_GO_BP(degs_sets_liver_D22[[lbl]], lbl)
})

ORA_GO_liver_D22 <- bind_rows(ORA_GO_liver_D22_list) %>%
  mutate(level = "all_DEGs") 

ORA_GO_liver_D22_PlotmePanel <- ORA_PublishFEBarplot(ORA_GO_liver_D22, n_top = 10, ncol = 3, title_col = "contrast")

ORA_PublishFEBarplot(ORA_GO_liver_D22, n_top = 10, ncol = 3, title_col = "contrast",
                     order_levels = GeneSet_PairSet_target$set1)

## EXTRA

## one can remove redundant enriched GO terms using similarity matrix analysis via simplifyEnrichment::simplify()
## this mainly involve clustering redundant terms by semantic similarity, then retaining only the most significant term per cluster per contrast
## which reduces GO redundancy before plotting the top 10 enriched terms per contrast

sem_bp_gga <- GOSemSim::godata(annoDb = "org.Gg.eg.db", ont = "BP", computeIC = FALSE)
clustermap  <- simplifyEnrichment::simplifyGO(unique(ORA_GO_liver_D22$ID),
                                              ont = "BP", sem_data = sem_bp_gga, cutoff = 0.7, plot = FALSE)
ORA_GO_liver_D22_simp <- ORA_GO_liver_D22 %>%
  dplyr::left_join(clustermap, by = c("ID" = "id")) %>%
  dplyr::group_by(contrast, cluster) %>%
  dplyr::slice_min(order_by = p.adjust, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(level = "all_DEGs")
ORA_GO_liver_D22_PlotmePanel_simp <- ORA_PublishFEBarplot(ORA_GO_liver_D22_simp, n_top = 10, ncol = 3, title_col = "contrast")
ORA_PublishFEBarplot(ORA_GO_liver_D22_simp, n_top = 10, ncol = 3, title_col = "contrast",
                     order_levels = GeneSet_PairSet_target$set1)

# Supp_tables -------------------------------------------------------------

# For your supp materials, ensure all output files, tables, and figures are systematically named and reproted for better computational reproducibility

## convert all Supp_tables to proper data.frame class, if needed

TxiSalmon_liver_abund <- as.data.frame(TxiSalmon_liver_D22$abundance) %>%
  rownames_to_column("isoform_id") %>%
  rename_with(~ paste0(.x, "_TPM"), -isoform_id)
TxiSalmon_liver_counts <- as.data.frame(TxiSalmon_liver_D22$counts) %>%
  rownames_to_column("isoform_id") %>%
  rename_with(~ paste0(.x, "_count"), -isoform_id)
TxiSalmon_liver_ExpMatrix <- inner_join(TxiSalmon_liver_abund, TxiSalmon_liver_counts, by = "isoform_id")
rm(TxiSalmon_liver_abund,TxiSalmon_liver_counts)

TxiSalmon_liver_ExpMatrix                                          <- as.data.frame(TxiSalmon_liver_ExpMatrix)
DEG_full_info_liver_D22_all_anno                                   <- as.data.frame(DEG_full_info_liver_D22_all_anno)
ORA_GO_liver_D22                                                   <- as.data.frame(ORA_GO_liver_D22)

## create the excel file table of content

Rdataframes_titles <- data.frame(
  Sheet = c(
    "Sheet1",
    "Sheet2",
    "Sheet3"
  ),
  
  DataFrameName = c(
    "TxiSalmon_liver_ExpMatrix",
    "DEG_full_info_liver_D22_all_anno",
    "ORA_GO_liver_D22"
  ),
  
  DataFrameNote = c(
    "Salmon expression matrix results at the level of transcripts with raw and TPM counts",
    "Differentially Expressed Genes (DEGs) for 6 group comparisons - Day 22 post-hatch",
    "ORA for GO terms results (Biological Process) - 6 group comparisons - Day 22 post-hatch"
  )
)

## create the excel file

write_xlsx(list(
  
  ContentTable = Rdataframes_titles,
  
  Sheet1  = TxiSalmon_liver_ExpMatrix_v2,
  Sheet2  = DEG_full_info_liver_v2_Bias_D22_lfcShrink_normal_all_anno,
  Sheet3  = ORA_GO_liver_v2_Bias_D22_lfcShrink_normal
),

file.path(sys_dir,"transcriptomics_training/supp_tables/SuppTables_Transcriptomics_Liver_v1.xlsx")
)


