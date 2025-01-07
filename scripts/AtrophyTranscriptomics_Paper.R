#### Transcriptomics analysis - Muscle Atrophy Project
#### Shadi Shahatit - RA, JUST 2024-2025
# libraries ---------------------------------------------------------------

library(DT)
library(data.table)
library(RcisTarget)
library(cowplot)
library(ggsci)
library(pathview)
library(writexl)
library(gridExtra)
library(DOSE)
library(broom)
library(ggvenn)
library(WGCNA)
library(igraph)
library(BiocManager)
library(DESeq2)
library(edgeR)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(readxl)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ReactomePA)
library(enrichplot)
library(ggraph)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Summary statistics - Macrogen --------------------------------------------

## define your stats data for raw, trimmed, and mapped reads

Raw_data_stats <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6),
  Sample_id = c("Control1", "treatment1-1", "treatment1-2", "treatment2-1", "treatment2-3", "treatment3-2"),
  Total_read_bases = c(6082470278, 4998486162, 4483018926, 4940785872, 4769793478, 4639114426),
  Total_reads = c(60222478, 49489962, 44386326, 48918672, 47225678, 45931826),
  GC_percent = c(46.99, 46.38, 47.15, 46.50, 46.44, 46.40),
  Q20_percent = c(97.73, 98.10, 98.01, 98.11, 98.32, 98.39),
  Q30_percent = c(94.08, 95.08, 94.50, 95.17, 95.24, 95.53))

Raw_data_stats$Sample_id <- gsub("Control1", "Control_1", Raw_data_stats$Sample_id)
Raw_data_stats$Sample_id <- gsub("treatment1-1", "Atrophy_1", Raw_data_stats$Sample_id)
Raw_data_stats$Sample_id <- gsub("treatment1-2", "Atrophy_2", Raw_data_stats$Sample_id)
Raw_data_stats$Sample_id <- gsub("treatment2-1", "Atrophy_Exercise_1", Raw_data_stats$Sample_id)
Raw_data_stats$Sample_id <- gsub("treatment2-3", "Atrophy_Exercise_2", Raw_data_stats$Sample_id)
Raw_data_stats$Sample_id <- gsub("treatment3-2", "Exercise_1", Raw_data_stats$Sample_id)

Trimming_data_stats <- data.frame(
  Sample_id = c("Control1", "treatment1-1", "treatment1-2", "treatment2-1", "treatment2-3", "treatment3-2"),
  Total_read_bases = c(5888398645, 4849899891, 4374154281, 4790122852, 4666328481, 4532874044),
  Total_reads = c(58667774, 48269198, 43572456, 47680060, 46434878, 45088010),
  GC_percent = c(46.98, 46.38, 47.15, 46.51, 46.44, 46.41),
  Q20_percent = c(98.60, 98.92, 98.65, 98.95, 98.88, 98.99),
  Q30_percent = c(95.20, 96.15, 95.35, 96.26, 96.01, 96.33))

Trimming_data_stats$Sample_id <- gsub("Control1", "Control_1", Raw_data_stats$Sample_id)
Trimming_data_stats$Sample_id <- gsub("treatment1-1", "Atrophy_1", Trimming_data_stats$Sample_id)
Trimming_data_stats$Sample_id <- gsub("treatment1-2", "Atrophy_2", Trimming_data_stats$Sample_id)
Trimming_data_stats$Sample_id <- gsub("treatment2-1", "Atrophy_Exercise_1", Trimming_data_stats$Sample_id)
Trimming_data_stats$Sample_id <- gsub("treatment2-3", "Atrophy_Exercise_2", Trimming_data_stats$Sample_id)
Trimming_data_stats$Sample_id <- gsub("treatment3-2", "Exercise_1", Trimming_data_stats$Sample_id)

Mapping_data_stats <- data.frame(
  Index = c(1, 2, 3, 4, 5, 6),
  Sample_id = c("Control1", "treatment1-1", "treatment1-2", "treatment2-1", "treatment2-3", "treatment3-2"),
  Processed_reads = c(58667774, 48269198, 43572456, 47680060, 46434878, 45088010),
  Mapped_reads = c(57045527, 46307713, 42354910, 45806114, 44850572, 43588553),
  Mapped_percent = c(97.23, 95.94, 97.21, 96.07, 96.59, 96.67),
  Unmapped_reads = c(1622247, 1961485, 1217546, 1873946, 1584306, 1499457),
  Unmapped_percent = c(2.77, 4.06, 2.79, 3.93, 3.41, 3.33))

Mapping_data_stats$Sample_id <- gsub("Control1", "Control_1", Mapping_data_stats$Sample_id)
Mapping_data_stats$Sample_id <- gsub("treatment1-1", "Atrophy_1", Mapping_data_stats$Sample_id)
Mapping_data_stats$Sample_id <- gsub("treatment1-2", "Atrophy_2", Mapping_data_stats$Sample_id)
Mapping_data_stats$Sample_id <- gsub("treatment2-1", "Atrophy_Exercise_1", Mapping_data_stats$Sample_id)
Mapping_data_stats$Sample_id <- gsub("treatment2-3", "Atrophy_Exercise_2", Mapping_data_stats$Sample_id)
Mapping_data_stats$Sample_id <- gsub("treatment3-2", "Exercise_1", Mapping_data_stats$Sample_id)

Raw_Trim_data <- merge(Raw_data_stats, Trimming_data_stats, by = "Sample_id", suffixes = c("_raw", "_trimmed"))
Raw_Trim_Map_data <- merge(Raw_Trim_data, Mapping_data_stats, by = "Sample_id", suffixes = c("", "_mapped"))

Sample_id_order <- c("Control_1", 
                     "Atrophy_1","Atrophy_2",
                     "Atrophy_Exercise_1","Atrophy_Exercise_2",
                     "Exercise_1")

## get the read counts in GB

Raw_Trim_Map_data$Total_read_bases_raw_GB <- Raw_Trim_Map_data$Total_read_bases_raw/1e9
Raw_Trim_Map_data$Total_read_bases_trimmed_GB <- Raw_Trim_Map_data$Total_read_bases_trimmed/1e9
Raw_Trim_Map_data$Processed_reads_GB <- Raw_Trim_Map_data$Processed_reads/1e9
Raw_Trim_Map_data$Mapped_reads_GB <- Raw_Trim_Map_data$Mapped_reads/1e9

## restructure the data for ggplot and get the bar plots
## filter for target samples

Sample_id_order <- c("Control_1","Atrophy_1","Atrophy_2")

Raw_Trim_Map_data <- Raw_Trim_Map_data %>%
  filter(Sample_id %in% c("Control_1", "Atrophy_1", "Atrophy_2"))

throughput_Raw_Trim_Map <- rbind(
  data.frame(Sample_id = Raw_Trim_Map_data$Sample_id, Throughput = Raw_Trim_Map_data$Total_read_bases_raw_GB, Reads = "Raw"),
  data.frame(Sample_id = Raw_Trim_Map_data$Sample_id, Throughput = Raw_Trim_Map_data$Total_read_bases_trimmed_GB, Reads = "Trimmed"))

ggplot(throughput_Raw_Trim_Map, aes(y = Sample_id, x = Throughput, fill = Reads)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Throughput)), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) +
  scale_fill_manual(values =c("#4a94cf","#f4814a")) +
  labs(title = "", y = "Samples", x = "Throughput (Gb)") +
  scale_y_discrete(limits = rev(Sample_id_order)) +
  theme_classic()

q30_Raw_Trim_Map <- rbind(
  data.frame(Sample_id = Raw_Trim_Map_data$Sample_id, Q30_percent = Raw_Trim_Map_data$Q30_percent_raw, Reads = "Raw"),
  data.frame(Sample_id = Raw_Trim_Map_data$Sample_id, Q30_percent = Raw_Trim_Map_data$Q30_percent_trimmed, Reads = "Trimmed"))

ggplot(q30_Raw_Trim_Map, aes(y = Sample_id, x = Q30_percent, fill = Reads)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f", Q30_percent)), position = position_dodge(width = 0.9), hjust = -0.2, size = 3) +
  scale_fill_manual(values =c("#4a94cf","#f4814a")) +
  labs(title = "", y = "Samples", x = "Q30 (%)") +
  scale_y_discrete(limits = rev(Sample_id_order)) +
  theme_classic()

ggplot(Raw_Trim_Map_data, aes(y = Sample_id, x = Mapped_percent)) +
  geom_bar(stat = "identity", fill = "#6a7f90") +
  geom_text(aes(label = sprintf("%.1f", Mapped_percent)), hjust = -0.2, size = 3) +
  labs(title = "", y = "Samples", x = "Mapping Ratio (%)") + 
  scale_y_discrete(limits = rev(Sample_id_order)) +
  theme_classic()

# Final - DEGs_visualizations_Enrichment ----------------------------------

## define your RNEseq count results

file_path <- "/home/shadi/Desktop/transcriptome_alzghoul_SS/2023 transcriptome muscle analysis/transcriptome file analysis/result_RNAseq_excel/Expression_profile/StringTie/Expression_Profile.mm10.gene.xlsx"
expression_data <- read_xlsx(file_path)
head(expression_data)

rawcount_data <- expression_data %>%
  dplyr::select(starts_with("Control1") & ends_with("Read_Count") | starts_with("treatment") & ends_with("Read_Count")) %>%
  as.data.frame()
rownames(rawcount_data) <- expression_data$Gene_Symbol ## replace with Gene_Symbol if you need

## create your metadata

sample_ids <- colnames(rawcount_data)
group_labels <- c("Control", "Treatment1", "Treatment1","Treatment2", "Treatment2", "Treatment3")
metadata <- data.frame(
  sample_id = sample_ids,
  group = factor(group_labels, levels = c("Control", "Treatment1", "Treatment2", "Treatment3")))

metadata$group <- factor(metadata$group, levels = c("Control", "Treatment1", "Treatment2", "Treatment3"))
head(rawcount_data)
head(metadata)

## create a DESeq object

dds <- DESeqDataSetFromMatrix(
  countData = rawcount_data,
  colData = metadata,
  design = ~ group
)

# Pre-filter low-expression genes
dds <- dds[rowSums(counts(dds)) > 10, ]

## DESeq modeling

# Perform DESeq with LRT for global analysis
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)

## get DEGs for pair wise comparisons

DEG_full_info_list <- list()
AEG_full_info_list <- list()
groups <- levels(metadata$group)
reference <- "Control"
pairwise_comparisons <- combn(groups, 2, simplify = FALSE)
all_comparisons <- list()
for (pair in pairwise_comparisons) {
  group1 <- pair[1]
  group2 <- pair[2]
  all_comparisons <- c(all_comparisons, list(c(group1, group2), c(group2, group1)))
}
all_comparisons_filtered <- all_comparisons[!sapply(all_comparisons, function(x) x[1] == "Control")]

for (pair in all_comparisons_filtered) {
  group1 <- pair[1]
  group2 <- pair[2]
  res <- results(dds, contrast = c("group", group1, group2))
  # create differential expressed genes df for each comp
  degs <- as.data.frame(res) %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  degs$Gene_Symbol <- rownames(degs)
  # degs$comparison <- paste(group2, "vs", group1)
  DEG_full_info_list[[paste(group1, "vs", group2)]] <- degs
  # create all expressed genes df for each comp
  aegs <- as.data.frame(res) 
  aegs$Gene_Symbol <- rownames(aegs)
  AEG_full_info_list[[paste(group1, "vs", group2)]] <- aegs
}
DEG_full_info_df <- do.call(rbind, lapply(names(DEG_full_info_list), function(group) {
  cbind(group = group, DEG_full_info_list[[group]])
}))
AEG_full_info_df <- do.call(rbind, lapply(names(AEG_full_info_list), function(group) {
  cbind(group = group, AEG_full_info_list[[group]])
}))

DEG_full_info_df$group <- factor(DEG_full_info_df$group, levels = c("Treatment1 vs Control",
                                                                    "Treatment1 vs Treatment2",
                                                                    "Treatment1 vs Treatment3",
                                                                    "Treatment2 vs Control",
                                                                    "Treatment2 vs Treatment1",
                                                                    "Treatment2 vs Treatment3",
                                                                    "Treatment3 vs Control",
                                                                    "Treatment3 vs Treatment1",
                                                                    "Treatment3 vs Treatment2"))

# add directionality
DEG_full_info_df <- DEG_full_info_df %>%
  mutate(direction = ifelse(log2FoldChange > 0, "upregulated", "downregulated"))

DEG_full_info_df$group <- gsub("Control", "Control", DEG_full_info_df$group)
DEG_full_info_df$group <- gsub("Treatment1", "Atrophy", DEG_full_info_df$group)
DEG_full_info_df$group <- gsub("Treatment2", "Atrophy_Exercise", DEG_full_info_df$group)
DEG_full_info_df$group <- gsub("Treatment3", "Exercise", DEG_full_info_df$group)

DEG_full_info_df$group <- factor(DEG_full_info_df$group, levels = c("Atrophy vs Control",
                                                                    "Atrophy vs Atrophy_Exercise",
                                                                    "Atrophy vs Exercise",
                                                                    "Atrophy_Exercise vs Control",
                                                                    "Atrophy_Exercise vs Atrophy",
                                                                    "Atrophy_Exercise vs Exercise",
                                                                    "Exercise vs Control",
                                                                    "Exercise vs Atrophy",
                                                                    "Exercise vs Atrophy_Exercise"))

DEG_full_info_df_6comp <- DEG_full_info_df %>% filter(group == "Atrophy vs Control" | 
                                                        group == "Exercise vs Control" |
                                                        group == "Atrophy_Exercise vs Control" |
                                                        group == "Atrophy_Exercise vs Atrophy" |
                                                        group == "Atrophy_Exercise vs Exercise" |
                                                        group == "Atrophy vs Exercise")

DEG_df_sum_6comp <- DEG_full_info_df_6comp %>%
  group_by(group, direction) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
  mutate(total = upregulated + downregulated)

AEG_full_info_df <- AEG_full_info_df %>%
  mutate(direction = ifelse(log2FoldChange > 0, "upregulated", "downregulated"))
AEG_full_info_df$group <- gsub("Control", "Control", AEG_full_info_df$group)
AEG_full_info_df$group <- gsub("Treatment1", "Atrophy", AEG_full_info_df$group)
AEG_full_info_df$group <- gsub("Treatment2", "Atrophy_Exercise", AEG_full_info_df$group)
AEG_full_info_df$group <- gsub("Treatment3", "Exercise", AEG_full_info_df$group)
AEG_full_info_df_6comp <- AEG_full_info_df %>% filter(group == "Atrophy vs Control" | 
                                                        group == "Exercise vs Control" |
                                                        group == "Atrophy_Exercise vs Control" |
                                                        group == "Atrophy_Exercise vs Atrophy" |
                                                        group == "Atrophy_Exercise vs Exercise" |
                                                        group == "Atrophy vs Exercise")

## Bar Plots

DEG_sum_6comp_long <- DEG_df_sum_6comp %>%
  pivot_longer(cols = c("upregulated", "downregulated", "total"),
               names_to = "category", values_to = "count") %>% 
  mutate(category = factor(category, levels = c("downregulated", "upregulated", "total")))

ggplot(subset(DEG_sum_6comp_long), aes(x = group, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  theme_classic() +
  labs(title = "",
       x = "Conditions",
       y = "DEGs counts",
       fill = "DEGs category") +
  scale_fill_manual(values = c("#FFD662FF", "#00539CFF", "#949398FF")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# or all at once

DEG_sum_long <- DEG_sum_long %>%
  mutate(shared_agent = case_when(
    grepl("Control", group) ~ "Control",
    grepl("Atrophy vs", group) ~ "Atrophy",
    grepl("Atrophy_Exercise vs", group) ~ "Atrophy_Exercise",
    grepl("Exercise vs", group) ~ "Exercise",
    TRUE ~ "Other"))

shared_agents <- unique(DEG_sum_long$shared_agent)
plots_list <- list()

for (agent in shared_agents) {
  p <- ggplot(subset(DEG_sum_long, shared_agent == agent), aes(x = group, y = count, fill = category)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
    theme_classic() +
    labs(x = "",
         y = "DEGs counts",
         fill = "DEGs category") +
    scale_fill_manual(values = c("#FFD662FF", "#00539CFF", "#949398FF")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(
      # plot.title = element_text(size = 16, face = "bold"),
      # axis.title.x = element_text(size = 14, face = "bold"),
      # axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "none"  # Remove the legend
    )
  plots_list[[agent]] <- p
}
grid.arrange(plots_list[[1]], plots_list[[2]], plots_list[[3]], plots_list[[4]], ncol = 4)

## volcano plots

padj_threshold <- 0.05
logFC_threshold <- 1
generate_nice_volcano_plot <- function(data, treatment) {
  data$Significance <- ifelse(data$log2FoldChange > logFC_threshold & data$padj < padj_threshold, 
                              "Upregulated",
                              ifelse(data$log2FoldChange < -logFC_threshold & data$padj < padj_threshold, 
                                     "Downregulated", "Not Significant"))
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Not Significant" = "#949398FF", "Upregulated" = "#00539CFF", "Downregulated" = "#FFD662FF")) +
    theme_minimal() +
    labs(title = paste(treatment), x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          axis.title = element_blank()) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
}
results_names <- unique(AEG_full_info_df_6comp$group)
plots <- list()
for (cond in results_names) {
  cond_data <- AEG_full_info_df_6comp %>% filter(group == cond)
  cond_plot <- generate_nice_volcano_plot(cond_data, cond)
  plots[[cond]] <- cond_plot
}

combined_volcano_plot <- wrap_plots(plots) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(hjust = 0.5))) &
  theme(
    legend.position = "top",
    plot.margin = margin(5, 5, 5, 5),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12))
combined_volcano_plot <- combined_volcano_plot &
  labs(x = "Log2 Fold Change", y = "-Log10 P-value adj")
print(combined_volcano_plot)

## Venn and upset diagrams

## upset

gene_sets_final_upset <- DEG_full_info_df_6comp %>%
  group_by(group) %>%
  summarize(Gene_Symbols = list(unique(Gene_Symbol)), .groups = "drop") %>%
  deframe()

gene_sets_matrix <- fromList(gene_sets_final_upset)

upset(gene_sets_matrix, 
      main.bar.color = "#404080", 
      sets.bar.color = "#69b3a2",
      order.by = "freq", 
      sets = names(gene_sets_final),
      set_size.show = TRUE)

## Venn

gene_sets_final <- DEG_full_info_df_6comp %>%
  group_by(group) %>%
  summarize(Gene_Symbols = list(unique(Gene_Symbol)), .groups = "drop") %>%
  deframe()

ggVennDiagram(gene_sets_final, label = "none",
              set_color = viridis(6))

names(gene_sets_final)
gene_sets_final_2 <- gene_sets_final[-6]

venn <- venn.diagram(
  x = gene_sets_final_2,
  category.names = names(gene_sets_final_2),
  filename = NULL,
  fill = viridis(5),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.fontface = "bold",
  margin = 0.1)

grid.newpage()
grid.draw(venn)

ggvenn(
  gene_sets_final, 
  fill_color = viridis(6),
  stroke_size = 0.5, 
  set_name_size = 4, 
  text_size = 3)

## unique and intersections

DEG_full_info_df_6comp$group <- as.character(DEG_full_info_df_6comp$group)
conditions <- unique(DEG_full_info_df_6comp$group)
DEGs_lists <- lapply(conditions, function(condition) {
  DEG_full_info_df_6comp$Gene_Symbol[DEG_full_info_df_6comp$group == condition]
})
names(DEGs_lists) <- conditions
generate_combinations <- function(conditions, DEGs_lists) {
  result_list <- list()
  
  for (n in 2:length(conditions)) {  # Start from 2 to ensure combinations are of size >=2
    combinations <- combn(conditions, n, simplify = FALSE)
    for (comb in combinations) {
      selected_genes <- DEGs_lists[comb]
      shared_genes <- Reduce(intersect, selected_genes)
      
      if (length(shared_genes) > 0) {  # Only keep combinations with shared genes
        result_list[[paste(comb, collapse = " ∩ ")]] <- list(
          shared_genes = shared_genes
        )
      }
    }
  }
  return(result_list)
}
intersection_results_list <- generate_combinations(conditions, DEGs_lists)
intersection_df <- data.frame(
  group_pair = character(),
  shared_genes = character(),
  DEGs_count = integer(),
  num_intersections = integer(),
  stringsAsFactors = FALSE)
for (group_pair in names(intersection_results_list)) {
  shared_genes <- intersection_results_list[[group_pair]]$shared_genes
  shared_genes_str <- paste(shared_genes, collapse = ", ")
  DEGs_count <- length(shared_genes)
  num_intersections <- str_count(group_pair, "∩") + 1
  intersection_df <- rbind(intersection_df, data.frame(
    group_pair = group_pair,
    shared_genes = shared_genes_str,
    DEGs_count = DEGs_count,
    num_intersections = num_intersections,
    stringsAsFactors = FALSE
  ))
}
unique_DEGs_results <- list()
for (condition in conditions) {
  current_genes <- DEGs_lists[[condition]]
  other_genes <- unlist(DEGs_lists[setdiff(conditions, condition)])
  unique_genes <- setdiff(current_genes, other_genes)
  unique_DEGs_results[[condition]] <- unique_genes
}
unique_df <- data.frame(
  group = names(unique_DEGs_results),
  unique_genes = sapply(unique_DEGs_results, function(x) paste(x, collapse = ", ")),
  DEGs_count = sapply(unique_DEGs_results, length),
  stringsAsFactors = FALSE)

intersection_DEGs_df <- intersection_df
unique_DEGsdf <- unique_df

## enrichment

## over-representation analysis ORA (GO)

## all DEGs

perform_ORA_GO <- function(gene_set, treatment_name, ontologies = c("BP", "MF", "CC")) {
  results <- lapply(ontologies, function(ont) {
    enrich_result <- enrichGO(
      gene = gene_set,
      OrgDb = org.Mm.eg.db,
      ont = ont, 
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL", # depends on your gene IDs
      universe = rownames(rawcount_data) # all genes expressed as background
    )
    as.data.frame(enrich_result) %>%
      mutate(ontology = ont, treatment = treatment_name)
  })
  bind_rows(results)
}

DEG_full_info_df_6comp$group <- as.character(DEG_full_info_df_6comp$group)
deg_sets <- DEG_full_info_df_6comp %>%
  group_by(group) %>%
  summarize(Gene_Symbols = list(unique(Gene_Symbol)), .groups = "drop") %>%
  deframe()
ORA_GO_alldeg_level <- lapply(names(deg_sets), function(treatment) {
  perform_ORA_GO(deg_sets[[treatment]], treatment)
})
ORA_GO_alldeg_level_df <- bind_rows(ORA_GO_alldeg_level)
ORA_GO_alldeg_level_df$level <- "all_DEGs"

## compareCluster to compare enrichment pathways at once

ORA_GO_alldeg_level_df_bycond_list <- ORA_GO_alldeg_level_df %>%
  group_by(treatment) %>%
  summarise(
    gene_symbols = list(unique(unlist(strsplit(geneID, "/")))), 
    .groups = "drop"
  ) %>%
  deframe()

ck <- compareCluster(geneCluster = ORA_GO_alldeg_level_df_bycond_list,
                     fun = enrichGO,
                     OrgDb = org.Mm.eg.db,
                     keyType="SYMBOL")
dotplot(ck)
cnetplot(ck)

## unique DEGs

unique_DEG_sets <- unique_DEGs_results
ORA_GO_unique_level <- lapply(names(unique_DEG_sets), function(treatment) {
  perform_ORA_GO(unique_DEG_sets[[treatment]], treatment)
})
ORA_GO_unique_level_df <- bind_rows(ORA_GO_unique_level)
ORA_GO_unique_level_df$level <- "unique_DEGs"

## pairwise intersected DEGs

pairwise_intersected_DEG_sets_1 <- intersection_results_list[
  sapply(names(intersection_results_list), function(name) {
    sum(strsplit(name, "∩")[[1]] != "") == 2  # Condition: exactly one intersection (∩)
  })
]
pairwise_intersected_DEG_sets <- setNames(
  lapply(pairwise_intersected_DEG_sets_1, function(item) item$shared_genes),
  names(pairwise_intersected_DEG_sets_1))
ORA_GO_pairwise_intersected_level <- lapply(names(pairwise_intersected_DEG_sets), function(treatment) {
  perform_ORA_GO(pairwise_intersected_DEG_sets[[treatment]], treatment)
})
ORA_GO_pairwise_intersected_level_df <- bind_rows(ORA_GO_pairwise_intersected_level)
ORA_GO_pairwise_intersected_level_df$level <- "pairwise_intersected_DEGs"

## high order intersected DEGs

highorder_intersected_DEG_sets_1 <- intersection_results_list[
  sapply(names(intersection_results_list), function(name) {
    sum(strsplit(name, "∩")[[1]] != "") > 1  # Condition: more than one intersection (∩)
  })
]
highorder_intersected_DEG_sets <- setNames(
  lapply(highorder_intersected_DEG_sets_1, function(item) item$shared_genes),
  names(highorder_intersected_DEG_sets_1))
ORA_GO_highorder_intersected_level <- lapply(names(highorder_intersected_DEG_sets), function(treatment) {
  perform_ORA_GO(highorder_intersected_DEG_sets[[treatment]], treatment)
})
ORA_GO_highorder_intersected_level_df <- bind_rows(ORA_GO_highorder_intersected_level)
ORA_GO_highorder_intersected_level_df$level <- "highorder_intersected_DEGs"

## filter and organize results

ORA_GO_combined_df <- bind_rows(ORA_GO_alldeg_level_df, ORA_GO_unique_level_df,
                                ORA_GO_pairwise_intersected_level_df, ORA_GO_highorder_intersected_level_df)

ORA_GO_combined_top5 <- ORA_GO_combined_df %>%
  group_by(treatment, ontology) %>%
  top_n(5, wt = -log10(pvalue)) %>%
  ungroup()

ORA_GO_summarized <- ORA_GO_combined_top5 %>%
  group_by(treatment, ontology) %>%
  summarise(pathways = paste(Description[order(pvalue)], collapse = ", "), .groups = "drop") %>%
  pivot_wider(names_from = ontology, values_from = pathways)

write_xlsx(ORA_GO_summarized, path = "ORA_GO_summarized.xlsx")

# View the result

ORA_top5_plotme <- ORA_GO_combined_top5 %>% filter(level == "all_DEGs")

ggplot(ORA_top5_plotme, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)), fill = -log10(pvalue))) + 
  geom_bar(stat = "identity", show.legend = FALSE, width = 0.8) + 
  scale_fill_gradient(low = "blue", high = "firebrick") +
  facet_grid(ontology ~ treatment, scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(title = "", x = "-log10(p-value)", y = "GO Term") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    strip.text = element_text(size = 4, color = "black"), 
    strip.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16), 
    axis.title = element_text(size = 14), 
    panel.grid.minor = element_blank()
  )

## gene set enrichment analysis GSEA (GO + KEGG)

## GO pathways

## get the ranked expression matrix for each treatment versus control and then perform the defined functions

rank_and_clean_genes <- function(genes_full_info_df, results_names) {
  ranked_gene_lists <- lapply(results_names, function(treatment) {
    genes_full_info_df %>%
      filter(group == treatment) %>%
      arrange(desc(log2FoldChange)) %>%
      pull(log2FoldChange, name = Gene_Symbol)  # Named vector: log2FoldChange as values, Gene_Symbol as names
  })
  names(ranked_gene_lists) <- results_names
  cleaned_ranked_gene_lists <- lapply(ranked_gene_lists, function(gene_list) {
    gene_list <- gene_list[!is.na(gene_list)]  # Remove NAs
    sort(gene_list, decreasing = TRUE)  # Sort by descending log2FoldChange
  })
  return(cleaned_ranked_gene_lists)
}
perform_GSEA_GO <- function(ranked_gene_lists, ontologies = c("BP", "MF", "CC"), pvalueCutoff = 0.05) {
  GSEA_results <- lapply(names(ranked_gene_lists), function(treatment) {
    lapply(ontologies, function(ont) {
      gseGO(
        geneList = ranked_gene_lists[[treatment]],
        ont = ont,
        keyType = "SYMBOL",
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        verbose = TRUE,
        OrgDb = org.Mm.eg.db,
        pAdjustMethod = "BH"
      )
    })
  })
  names(GSEA_results) <- names(ranked_gene_lists)
  return(GSEA_results)
}
plot_dotplot <- function(GSEA_results, treatment_index = 1, ontology_index = 1, showCategory = 10) {
  gsea_data <- GSEA_results[[treatment_index]][[ontology_index]]
  dotplot(gsea_data, showCategory = showCategory, split = ".sign") + 
    facet_grid(. ~ .sign) +
    ggtitle(paste0("GSEA Dotplot: Treatment ", names(GSEA_results)[treatment_index],
                   " | Ontology: ", names(GSEA_results[[treatment_index]])[ontology_index]))
}

## all DEGs

results_names <- unique(AEG_full_info_df_6comp$group)
cleaned_ranked_gene_lists <- rank_and_clean_genes(AEG_full_info_df_6comp, results_names)

GSEA_GO_alldeg_level <- perform_GSEA_GO(cleaned_ranked_gene_lists, ontologies = c("BP", "MF", "CC"))

# treatment_index are ordered based on R default for characters
# ontology_index are ordered as BP, MM, CC

plot_dotplot(GSEA_GO_alldeg_level, treatment_index = 1, ontology_index = 1, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level, treatment_index = 2, ontology_index = 1, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level, treatment_index = 3, ontology_index = 1, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level, treatment_index = 4, ontology_index = 1, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level, treatment_index = 5, ontology_index = 1, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level, treatment_index = 6, ontology_index = 1, showCategory = 10)

## unique DEGs

## pairwise intersected DEGs

## high order intersected DEGs

## filter and organize results
## View the result

# Final - only Atrophy - DEGs_visualizations_Enrichment ----------------------------------------------------

## filter only for Atrophy vs Control

AEG_full_info_df
DEG_full_info_df 

AEG_full_info_df_1comp_AtroCtrl <- AEG_full_info_df %>% filter(group == "Atrophy vs Control")
DEG_full_info_df_1comp_AtroCtrl <- DEG_full_info_df %>% filter(group == "Atrophy vs Control")

DEG_df_sum_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl %>%
  group_by(group, direction) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = count, values_fill = 0) %>%
  mutate(total = upregulated + downregulated)


cat(DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol, sep = "\n")

writeLines(DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol, "DEG_Gene_Symbol_AtroCtrl.txt")

## Bar Plots

DEG_sum_AtroCtrl_long <- DEG_df_sum_AtroCtrl %>%
  pivot_longer(cols = c("upregulated", "downregulated", "total"),
               names_to = "category", values_to = "count") %>% 
  mutate(category = factor(category, levels = c("downregulated", "upregulated", "total")))

ggplot(subset(DEG_sum_AtroCtrl_long), aes(x = group, y = count, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), vjust = -0.5, size = 3.5) +
  theme_classic() +
  labs(title = "",
       x = "Conditions",
       y = "DEGs counts",
       fill = "DEGs category") +
  scale_fill_manual(values = c("#FFD662FF", "#00539CFF", "#949398FF")) +
  theme(axis.text.x = element_text(angle = 0, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## Bar plot of top scoring genes

top_10_high_indices <- order(DEG_full_info_df_1comp_AtroCtrl$log2FoldChange, decreasing = TRUE)[1:10]
top_10_low_indices <- order(DEG_full_info_df_1comp_AtroCtrl$log2FoldChange)[1:10]

DEG_AtroCtrl_top10 <- rbind(DEG_full_info_df_1comp_AtroCtrl[top_10_high_indices, ], DEG_full_info_df_1comp_AtroCtrl[top_10_low_indices, ])

ggplot(subset(DEG_AtroCtrl_top10), aes(x=reorder(Gene_Symbol,log2FoldChange), y = log2FoldChange, fill = direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0, color = "black") +
  theme_classic() +
  labs(title = "",
       x = "DEGs",
       y = "Fold change",
       fill = "DEGs category") +
  scale_fill_manual(values = c("#FFD662FF", "#00539CFF")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## heat map of top scoring genes

DEG_AtroCtrl_top10$log2FoldChange_color <- ifelse(DEG_AtroCtrl_top10$log2FoldChange > 0,
                                                  "upregulated", "downregulated")

ggplot(DEG_AtroCtrl_top10, aes(x = reorder(Gene_Symbol, log2FoldChange), y = 1, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(high = "darkgreen", mid = "white", low = "firebrick", midpoint = 0, 
                       name = "Log2 Fold change") +
  labs(title = "", 
       x = "Gene symbol", 
       y = "Fold change") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))

## volcano plots

padj_threshold <- 0.05
logFC_threshold <- 1
generate_nice_volcano_plot <- function(data, treatment) {
  data$Significance <- ifelse(data$log2FoldChange > logFC_threshold & data$padj < padj_threshold, 
                              "Upregulated",
                              ifelse(data$log2FoldChange < -logFC_threshold & data$padj < padj_threshold, 
                                     "Downregulated", "Not Significant"))
  
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Not Significant" = "#949398FF", "Upregulated" = "#00539CFF", "Downregulated" = "#FFD662FF")) +
    theme_minimal() +
    labs(title = paste(treatment), x = NULL, y = NULL) +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          axis.title = element_blank()) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
}
results_names <- unique(AEG_full_info_df_1comp_AtroCtrl$group)
plots <- list()
for (cond in results_names) {
  cond_data <- AEG_full_info_df_1comp_AtroCtrl %>% filter(group == cond)
  cond_plot <- generate_nice_volcano_plot(cond_data, cond)
  plots[[cond]] <- cond_plot
}
combined_volcano_plot <- wrap_plots(plots) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "",
    theme = theme(plot.title = element_text(hjust = 0.5))) &
  theme(
    legend.position = "top",
    plot.margin = margin(5, 5, 5, 5),
    axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
combined_volcano_plot <- combined_volcano_plot &
  labs(x = "Log2 Fold Change", y = "-Log10 P-value adj")
print(combined_volcano_plot)

## Venn and upset diagrams

## Venn

gene_sets_final_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl %>%
  group_by(group) %>%
  summarize(Gene_Symbols = list(unique(Gene_Symbol)), .groups = "drop") %>%
  deframe()

ggvenn(
  gene_sets_final_AtroCtrl, 
  fill_color = viridis(2),
  stroke_size = 0.5, 
  set_name_size = 4, 
  text_size = 3)

## enrichment

## over-representation analysis ORA (GO) "Atrophy vs Control"

## all DEGs 

perform_ORA_GO <- function(gene_set, treatment_name, ontologies = c("BP", "MF", "CC")) {
  results <- lapply(ontologies, function(ont) {
    enrich_result <- enrichGO(
      gene = gene_set,
      OrgDb = org.Mm.eg.db,
      ont = ont, 
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      keyType = "SYMBOL", # depends on your gene IDs
      universe = rownames(rawcount_data) # all genes expressed as background
    )
    as.data.frame(enrich_result) %>%
      mutate(ontology = ont, treatment = treatment_name)
  })
  bind_rows(results)
}

gene_sets_final_AtroCtrl$group <- as.character(gene_sets_final_AtroCtrl$group)
deg_sets_AtroCtrl <- list(gene_sets_final_AtroCtrl$`Atrophy vs Control`)
names(deg_sets_AtroCtrl) <- "Atrophy vs Control"
ORA_GO_alldeg_level_AtroCtrl <- lapply(names(deg_sets_AtroCtrl), function(treatment) {
  perform_ORA_GO(deg_sets_AtroCtrl[[treatment]], treatment)
})
ORA_GO_alldeg_level_AtroCtrl_df <- bind_rows(ORA_GO_alldeg_level_AtroCtrl)
ORA_GO_alldeg_level_AtroCtrl_df$level <- "all_DEGs"

ORA_GO_AtroCtrl_top10 <- ORA_GO_alldeg_level_AtroCtrl_df %>%
  group_by(treatment, ontology) %>%
  top_n(10, wt = -log10(pvalue)) %>%
  ungroup()

ggplot(ORA_GO_AtroCtrl_top10, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)), fill = -log10(pvalue))) + 
  geom_bar(stat = "identity", show.legend = FALSE, width = 0.8) + 
  scale_fill_gradient(low = "blue", high = "firebrick") +
  facet_grid(ontology ~ treatment, scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(title = "", x = "-log10(p-value)", y = "GO Term") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    strip.text = element_text(size = 14, face = "bold", color = "black"), 
    strip.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14), 
    panel.grid.minor = element_blank()
  )

## over-representation analysis ORA (KEGG) "Atrophy vs Control"

## all DEGs 
## over-representation analysis ORA (KEGG) "Atrophy vs Control"

## all DEGs 

## get right ids for KEGG in universe and DEGs input

gene_symbols_rawcount_data <- rownames(rawcount_data)
geneIDs_universe <- bitr(geneID = gene_symbols_rawcount_data, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Mm.eg.db)
gene_symbols <- gene_sets_final_AtroCtrl$`Atrophy vs Control`
mapped_genes <- bitr(geneID = gene_symbols, 
                     fromType = "SYMBOL",
                     toType = "ENTREZID", 
                     OrgDb = org.Mm.eg.db)
if (nrow(mapped_genes) == 0) {
  stop("No genes were successfully mapped to Entrez IDs.")
}
valid_genes <- mapped_genes$ENTREZID
deg_sets_AtroCtrl_geneID <- list(valid_genes)
names(deg_sets_AtroCtrl_geneID) <- "Atrophy vs Control"

## run the analysis

perform_ORA_KEGG <- function(gene_set, treatment_name) {
  enrich_result <- enrichKEGG(
    gene = gene_set,
    organism = 'mmu', # For mouse
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    universe = geneIDs_universe$ENTREZID
  )
  as.data.frame(enrich_result) %>%
    mutate(treatment = treatment_name)
}

ORA_KEGG_alldeg_level_AtroCtrl <- lapply(names(deg_sets_AtroCtrl_geneID), function(treatment) {
  perform_ORA_KEGG(deg_sets_AtroCtrl_geneID[[treatment]], treatment)
})
ORA_KEGG_alldeg_level_AtroCtrl_df <- bind_rows(ORA_KEGG_alldeg_level_AtroCtrl)
ORA_KEGG_alldeg_level_AtroCtrl_df$level <- "all_DEGs"

ORA_KEGG_AtroCtrl_top <- ORA_KEGG_alldeg_level_AtroCtrl_df %>%
  group_by(treatment) %>%
  top_n(5, wt = -log10(pvalue)) %>%
  ungroup()

ggplot(ORA_KEGG_AtroCtrl_top, aes(x = -log10(pvalue), y = reorder(Description, -log10(pvalue)), fill = -log10(pvalue))) + 
  geom_bar(stat = "identity", show.legend = FALSE, width = 0.8) + 
  scale_fill_gradient(low = "blue", high = "firebrick") +
  facet_grid(. ~ treatment, scales = "free_y", space = "free_y") +
  theme_bw() +
  labs(title = "", x = "-log10(p-value)", y = "KEGG Pathway") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
    strip.text = element_text(size = 14, face = "bold", color = "black"), 
    strip.background = element_rect(fill = "white", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), 
    axis.title = element_text(size = 14), 
    panel.grid.minor = element_blank()
  )

## fold change stacked bar plot

FoldChange_DEG_full_info_df_1comp_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl %>%
  dplyr::select(Gene_Symbol, log2FoldChange)

pathway_genes_KEGG_AtroCtrl <- ORA_KEGG_AtroCtrl_top %>%
  dplyr::select(ID, Description, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  left_join(mapped_genes, by = c("geneID" = "ENTREZID")) %>%
  rename(Gene_Symbol = SYMBOL) %>%
  left_join(FoldChange_DEG_full_info_df_1comp_AtroCtrl, by = c("Gene_Symbol" = "Gene_Symbol"))

pathway_genes_KEGG_AtroCtrl_clean <- pathway_genes_KEGG_AtroCtrl %>%
  mutate(Description = str_replace_all(Description, 
                                       c("signaling pathway" = "", 
                                         "- Mus musculus \\(house mouse\\)" = ""))) %>%
  arrange(log2FoldChange)

assign_gene_colors <- function(df) {
  d3_colors <- scale_fill_d3(palette = "category20")$palette(22)
  df <- df %>%
    group_by(Description) %>%
    arrange(Gene_Symbol) %>% 
    mutate(gene_order = row_number()) %>%  
    ungroup()
  df$color <- d3_colors[(df$gene_order - 1) %% 22 + 1]
  return(df)
}

pathway_genes_KEGG_AtroCtrl_clean_colored <- assign_gene_colors(pathway_genes_KEGG_AtroCtrl_clean)

ggplot(data=subset(pathway_genes_KEGG_AtroCtrl_clean), 
       aes(x = Description, y = log2FoldChange, fill = log2FoldChange)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.6) +
  scale_fill_viridis(direction = -1) +
  labs(title = "",
       x = "KEGG Pathway",
       y = "log2 Fold change",
       fill = "Pathway DEGs") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "right") +
  facet_wrap(~Description, scales = "free_x", ncol = 10)

ggplot(data = subset(pathway_genes_KEGG_AtroCtrl_clean_colored), 
       aes(x = Description, y = log2FoldChange, fill = Gene_Symbol)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.6) +
  scale_fill_manual(values = setNames(pathway_genes_KEGG_AtroCtrl_clean_colored$color, 
                                      pathway_genes_KEGG_AtroCtrl_clean_colored$Gene_Symbol)) +  
  labs(title = "",
       x = "KEGG Pathway",
       y = "log2 Fold change",
       fill = "Gene Symbol") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        legend.position = "none",
        legend.title = element_text(size = 12),  
        legend.text = element_text(size = 10)) + 
  facet_wrap(~Description, scales = "free_x", ncol = 10)

get_legend <- function(plot, description) {
  g_legend <- ggplotGrob(plot)$grobs[[which(sapply(ggplotGrob(plot)$grobs, function(x) x$name) == "guide-box")]]
  legend_title <- textGrob(description, gp = gpar(fontsize = 14, fontface = "bold"))
  g_legend$grobs[[1]]$children[[1]] <- legend_title
  return(g_legend)
}
unique_descriptions <- unique(pathway_genes_KEGG_AtroCtrl_clean$Description)
plots <- list()
legends <- list()
for (desc in unique_descriptions) {
  plot <- ggplot(data = subset(pathway_genes_KEGG_AtroCtrl_clean_colored, Description == desc), 
                 aes(x = Description, y = log2FoldChange, fill = Gene_Symbol)) +
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", size = 0.6) +
    scale_fill_d3(palette = "category20") +
    labs(title = "", x = "", y = "log2 Fold change", fill = desc) +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 10), legend.position = "right")
  plots[[desc]] <- plot
  legends[[desc]] <- get_legend(plot, desc)
}
plot_grid(plotlist = legends, ncol = 5)

## revisit

geneList <- c("12345" = 1.5, "23456" = -1.2, "34567" = 0.8)  # Replace with your data
deg_sets_AtroCtrl_geneID
mmu04110 <- pathview(
  gene.data = deg_sets_AtroCtrl_geneID,       
  pathway.id = "mmu04110",    
  species = "mmu",            
  limit = list(gene = max(abs(geneList)), cpd = 1))

## gene set enrichment analysis GSEA (GO) "Atrophy vs Control"

## GO pathways

## get the ranked expression matrix and then perform the defined functions
## or use gseGO directly for simple lists

rank_and_clean_genes <- function(genes_full_info_df, results_names) {
  ranked_gene_lists <- lapply(results_names, function(treatment) {
    genes_full_info_df %>%
      filter(group == treatment) %>%
      arrange(desc(log2FoldChange)) %>%
      pull(log2FoldChange, name = Gene_Symbol)  # Named vector: log2FoldChange as values, Gene_Symbol as names
  })
  names(ranked_gene_lists) <- results_names
  cleaned_ranked_gene_lists <- lapply(ranked_gene_lists, function(gene_list) {
    gene_list <- gene_list[!is.na(gene_list)]  # Remove NAs
    sort(gene_list, decreasing = TRUE)  # Sort by descending log2FoldChange
  })
  return(cleaned_ranked_gene_lists)
}
perform_GSEA_GO <- function(ranked_gene_lists, ontologies = c("BP", "MF", "CC"), pvalueCutoff = 0.05) {
  GSEA_results <- lapply(names(ranked_gene_lists), function(treatment) {
    lapply(ontologies, function(ont) {
      gseGO(
        geneList = ranked_gene_lists[[treatment]],
        ont = ont,
        keyType = "SYMBOL",
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = pvalueCutoff,
        verbose = TRUE,
        OrgDb = org.Mm.eg.db,
        pAdjustMethod = "BH"
      )
    })
  })
  names(GSEA_results) <- names(ranked_gene_lists)
  return(GSEA_results)
}
plot_dotplot <- function(GSEA_results, treatment_index = 1, ontology_index = 1, showCategory = 10) {
  gsea_data <- GSEA_results[[treatment_index]][[ontology_index]]
  dotplot(gsea_data, showCategory = showCategory, split = ".sign") + 
    facet_grid(. ~ .sign) +
    ggtitle(paste0("GSEA Dotplot: Treatment ", names(GSEA_results)[treatment_index],
                   " | Ontology: ", names(GSEA_results[[treatment_index]])[ontology_index]))
}

## all DEGs

results_names_AtroCtrl <- "Atrophy vs Control"
cleaned_ranked_gene_lists_AtroCtrl <- rank_and_clean_genes(AEG_full_info_df, results_names_AtroCtrl)

GSEA_GO_alldeg_level_AtroCtrl <- gseGO(geneList=cleaned_ranked_gene_lists_AtroCtrl$`Atrophy vs Control`, 
                                       ont ="ALL", 
                                       keyType = "SYMBOL", 
                                       minGSSize = 3, 
                                       maxGSSize = 800, 
                                       pvalueCutoff = 0.05, 
                                       verbose = TRUE, 
                                       OrgDb = org.Mm.eg.db, 
                                       pAdjustMethod = "BH")

# foldChange_vector <- AEG_full_info_df_1comp_AtroCtrl %>% 
#   pull(log2FoldChange, name = Gene_Symbol)

foldChange <- cleaned_ranked_gene_lists_AtroCtrl[[results_names_AtroCtrl]]
summary(foldChange)

dotplot(GSEA_GO_alldeg_level_AtroCtrl, showCategory=10, split=".sign") + facet_grid(.~.sign)
cnetplot(GSEA_GO_alldeg_level_AtroCtrl,
         node_label= "category",
         categorySize="pvalue",
         foldChange=foldChange,
         showCategory = 10,
         layout="kk") +
  scale_color_gradient2(name = "fold change", high = "firebrick", mid = "gray50", low = "gray50", midpoint = 0)
cnetplot(GSEA_GO_alldeg_level_AtroCtrl, foldChange=cleaned_ranked_gene_lists_AtroCtrl, circular = TRUE, colorEdge = TRUE) 
edo <- pairwise_termsim(GSEA_GO_alldeg_level_AtroCtrl)
emapplot(edo, pie="count", cex_category=1, layout="kk")

# or with a defined function for larger lists

GSEA_GO_alldeg_level_AtroCtrl <- perform_GSEA_GO(cleaned_ranked_gene_lists_AtroCtrl, ontologies = c("BP", "MF", "CC"))

# treatment_index are ordered based on R default for characters
# ontology_index are ordered as BP, MM, CC

plot_dotplot(GSEA_GO_alldeg_level_AtroCtrl, treatment_index = 1, ontology_index = 1, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level_AtroCtrl, treatment_index = 1, ontology_index = 2, showCategory = 10)
plot_dotplot(GSEA_GO_alldeg_level_AtroCtrl, treatment_index = 1, ontology_index = 3, showCategory = 10)

plot_cnetplot <- function(GSEA_results, treatment_index = 1, ontology_index = 1, showCategory = 10, foldChange = NULL) {
  gsea_data <- GSEA_results[[treatment_index]][[ontology_index]]
  cnetplot(
    gsea_data, 
    categorySize = "pvalue", 
    showCategory = showCategory,
    foldChange = foldChange
  ) +
    ggtitle(paste0("GSEA Cnetplot: Treatment ", names(GSEA_results)[treatment_index],
                   " | Ontology: ", names(GSEA_results[[treatment_index]])[ontology_index]))
}
plot_cnetplot(GSEA_GO_alldeg_level_AtroCtrl, treatment_index = 1, ontology_index = 1, showCategory = 10, foldChange = foldChange_vector)

## gene set enrichment analysis GSEA (KEGG) "Atrophy vs Control"

## KEGG pathways

## get right ids for KEGG in universe and DEGs input

mapped_genes <- bitr(
  geneID = unique(AEG_full_info_df$Gene_Symbol), 
  fromType = "SYMBOL", 
  toType = "ENTREZID", 
  OrgDb = org.Mm.eg.db)

AEG_full_info_df_geneID <- AEG_full_info_df %>%
  left_join(mapped_genes, by = c("Gene_Symbol" = "SYMBOL"))

AEG_full_info_df_geneID_clean <- AEG_full_info_df_geneID[!is.na((AEG_full_info_df_geneID$ENTREZID)), ]

AEG_full_info_df_geneID_only <- AEG_full_info_df_geneID_clean %>%
  mutate(Gene_Symbol = ENTREZID) %>%
  dplyr::select(-ENTREZID)

## get the ranked expression matrix and then perform the defined functions
## or use gseGO directly for simple lists

## all DEGs

results_names_AtroCtrl <- "Atrophy vs Control"
cleaned_ranked_gene_lists_AtroCtrl_geneID <- rank_and_clean_genes(AEG_full_info_df_geneID_only, results_names_AtroCtrl)

GSEA_KEGG_alldeg_level_AtroCtrl <- gseKEGG(geneList=cleaned_ranked_gene_lists_AtroCtrl_geneID$`Atrophy vs Control`, 
                                           organism = 'mmu', 
                                           minGSSize = 3, 
                                           maxGSSize = 800, 
                                           verbose = TRUE,
                                           nPermSimple = 10000)

foldChange <- cleaned_ranked_gene_lists_AtroCtrl_geneID[[results_names_AtroCtrl]]
summary(foldChange)

dotplot(GSEA_KEGG_alldeg_level_AtroCtrl, showCategory=10, split=".sign") + facet_grid(.~.sign)
cnetplot(GSEA_KEGG_alldeg_level_AtroCtrl,
         node_label= "category",
         categorySize="pvalue",
         foldChange=foldChange,
         showCategory = 10,
         layout="kk") +
  scale_color_gradient2(name = "fold change", high = "firebrick", mid = "gray50", low = "gray50", midpoint = 0)
cnetplot(GSEA_KEGG_alldeg_level_AtroCtrl, foldChange=cleaned_ranked_gene_lists_AtroCtrl, circular = TRUE, colorEdge = TRUE) 
edo <- pairwise_termsim(GSEA_KEGG_alldeg_level_AtroCtrl)
emapplot(edo, pie="count", cex_category=1, layout="kk")

# Final - only Atrophy - TF analysis --------------------------------------

## TF analysis

## use expression matrix data to get read count

## TPM and FPKM are very similar, but to compare between samples TPM is better

expression_data_AtroCtrl <- expression_data %>% dplyr::select(c("Gene_Symbol","Description","treatment1-1_TPM","treatment1-2_TPM"))
expression_data_AtroCtrl$treatment1_avg_TPM <- (expression_data_AtroCtrl$`treatment1-1_TPM`+expression_data_AtroCtrl$`treatment1-2_TPM`)/2

topindhigh <- 50 ## 10, 25, or 50?

## ChEA3 TF analysis

ChEA3_Integrated_meanRank_DEG_AtroCtrl <- read.table("ChEA3_Integrated_meanRank_DEG_AtroCtrl.tsv", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
ChEA3_Integrated_meanRank_DEG_AtroCtrl <- data.frame(lapply(ChEA3_Integrated_meanRank_DEG_AtroCtrl, function(x) gsub('"', '', x)))
ChEA3_Integrated_meanRank_DEG_AtroCtrl$Score <- as.numeric(ChEA3_Integrated_meanRank_DEG_AtroCtrl$Score)
ChEA3_Integrated_meanRank_DEG_AtroCtrl$Rank <- as.numeric(ChEA3_Integrated_meanRank_DEG_AtroCtrl$Rank)
TF_ChEA3_Integrated_meanRank_DEG_AtroCtrl <- ChEA3_Integrated_meanRank_DEG_AtroCtrl[[3]]
top_indices <- seq(1:topindhigh) ## already ordered by rank
TF_AtroCtrl_ChEA3 <- TF_ChEA3_Integrated_meanRank_DEG_AtroCtrl[top_indices]
TF_AtroCtrl_ChEA3_mod <- sapply(TF_AtroCtrl_ChEA3, function(x) {paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))})

ggplot(ChEA3_Integrated_meanRank_DEG_AtroCtrl, aes(x = Rank, y = (Score))) +
  geom_point(aes(color = Score), size = 1.5) +
  scale_y_log10() +
  scale_color_gradient(high = "darkgreen", low = "firebrick") +
  labs(x = "TF rank", y = "Score", title = "ChEA3") +
  theme_bw()

## BART TF analysis

DEG_Gene_Symbol_AtroCtrl_bart_results <- read.table("/home/shadi/Desktop/transcriptome_alzghoul_SS/TF_bart/bart2-master/bart2output_AtroCtrl/DEG_Gene_Symbol_AtroCtrl_bart_results.txt", 
                                                    header = TRUE)
TF_DEG_Gene_Symbol_AtroCtrl_bart_results <- DEG_Gene_Symbol_AtroCtrl_bart_results[[1]]
top_indices <- seq(1:topindhigh) ## already ordered by rank
TF_AtroCtrl_BART <- TF_DEG_Gene_Symbol_AtroCtrl_bart_results[top_indices]
TF_AtroCtrl_BART_mod <- sapply(TF_AtroCtrl_BART, function(x) {paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))})

ggplot(DEG_Gene_Symbol_AtroCtrl_bart_results, aes(x = re_rank, y = -log(irwin_hall_pvalue))) +
  geom_point(aes(color = pvalue), size = 1.5) +
  scale_y_log10() +
  scale_color_gradient(high = "darkgreen", low = "firebrick") +
  labs(x = "TF rank", y = "P-value", title = "BRAT") +
  theme_bw()

## RcisTarget TF analysis

DEG_GeneSymbol_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol

# download the mm10 database files
# motif ranking from https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/
# motif annotations from https://resources.aertslab.org/cistarget/motif2tf/

motif_ranking_file <- "/home/shadi/Desktop/transcriptome_alzghoul_SS/mm10motif_db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motifRankings <- importRankings(motif_ranking_file)
class(motifRankings)
motifAnnotations <- importAnnotations("/home/shadi/Desktop/transcriptome_alzghoul_SS/mm10motif_db/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")
class(motifAnnotations)

## simple run

RcisTarget_AtroCtrl_motifEnrichmentTable_wGenes <- cisTarget(DEG_GeneSymbol_AtroCtrl,
                                                             motifRankings=motifRankings,
                                                             motifAnnot=motifAnnotations)

## advanced run

# calculate AUC
RcisTarget_AtroCtrl_motifs_AUC <- calcAUC(DEG_GeneSymbol_AtroCtrl, motifRankings)

auc <- getAUC(RcisTarget_AtroCtrl_motifs_AUC)
hist(auc, main="", xlab="AUC histogram", breaks=100, col="#ff000050", border="darkred")
nes3 <- (3*sd(auc)) + mean(auc)
abline(v=nes3, col="red")

# select significant motifs, add TF annotation & format as table
RcisTarget_AtroCtrl_motifEnrichmentTable <- addMotifAnnotation(RcisTarget_AtroCtrl_motifs_AUC, 
                                                               motifAnnot=motifAnnotations)
dim(RcisTarget_AtroCtrl_motifEnrichmentTable)
head(RcisTarget_AtroCtrl_motifEnrichmentTable[,-"TF_lowConf", with=FALSE])

# identify significant genes for each motif
RcisTarget_AtroCtrl_motifEnrichmentTable_wGenes <- addSignificantGenes(RcisTarget_AtroCtrl_motifEnrichmentTable, 
                                                                       geneSets=DEG_GeneSymbol_AtroCtrl,
                                                                       rankings=motifRankings, 
                                                                       nCores=1,
                                                                       method="iCisTarget")
selectedMotifs <- c(sample(RcisTarget_AtroCtrl_motifEnrichmentTable$motif, 2))
par(mfrow=c(2,2))
getSignificantGenes(DEG_GeneSymbol_AtroCtrl, 
                    rankings=motifRankings,
                    signifRankingNames=selectedMotifs,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="aprox")

resultsSubset <- RcisTarget_AtroCtrl_motifEnrichmentTable_wGenes[1:10,]
showLogo(resultsSubset)

# anotatedTfs <- lapply(split(RcisTarget_AtroCtrl_motifEnrichmentTable_wGenes$TF_highConf,
#                             RcisTarget_AtroCtrl_motifEnrichmentTable$geneSet),
#                       function(x) {
#                         genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
#                         genesSplit <- unique(unlist(strsplit(genes, "; ")))
#                         return(genesSplit)
#                       })
# anotatedTfs_mod <- sapply(anotatedTfs, function(x) {paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))})

TF_AtroCtrl_RcisTarget_mod <- {
  allTFs <- RcisTarget_AtroCtrl_motifEnrichmentTable_wGenes$TF_highConf
  genes <- gsub(" \\(.*\\). ", "; ", allTFs, fixed=FALSE) 
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesFormatted <- sapply(genesSplit, function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  })
  genesFormatted
}

gene_sets_TF <- list(
  ChEA3 = TF_AtroCtrl_ChEA3_mod,
  BART = TF_AtroCtrl_BART_mod,
  RcisTarget = TF_AtroCtrl_RcisTarget_mod
  # DEGs = DEG_GeneSymbol_AtroCtrl
)

ggvenn(gene_sets_TF, 
       fill_color = viridis(4),
       stroke_size = 0.5, 
       set_name_size = 4, 
       text_size = 3)

Reduce(intersect, gene_sets_TF)

# network
signifMotifNames <- RcisTarget_AtroCtrl_motifEnrichmentTable$motif[1:3]
incidenceMatrix <- getSignificantGenes(DEG_GeneSymbol_AtroCtrl, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)

## RcisTarget + background TF analysis

DEG_GeneSymbol_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol
AEG_GeneSymbol_AtroCtrl <- AEG_full_info_df_1comp_AtroCtrl$Gene_Symbol ## this is your bg

gplots::venn(list(Background=AEG_GeneSymbol_AtroCtrl, DEG=DEG_GeneSymbol_AtroCtrl))

motif_ranking_file <- "/home/shadi/Desktop/transcriptome_alzghoul_SS/mm10motif_db/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
motifRankings <- importRankings(motif_ranking_file)
class(motifRankings)
motifRankings_bg <- importRankings(motif_ranking_file, columns=AEG_GeneSymbol_AtroCtrl)
bgRanking <- reRank(motifRankings_bg) 

RcisTarget_AtroCtrl_motifEnrichmentTable_bg <- cisTarget(DEG_GeneSymbol_AtroCtrl, bgRanking, 
                                                         aucMaxRank=0.03*getNumColsInDB(bgRanking),
                                                         geneErnMaxRank=getNumColsInDB(bgRanking),
                                                         geneErnMethod = "icistarget",
                                                         motifAnnot=motifAnnotations)
showLogo(RcisTarget_AtroCtrl_motifEnrichmentTable_bg)

RcisTarget_AtroCtrl_motifs_bg_AUC <- calcAUC(DEG_GeneSymbol_AtroCtrl, bgRanking)
auc <- getAUC(RcisTarget_AtroCtrl_motifs_bg_AUC)
hist(auc, main="", xlab="AUC histogram", breaks=100, col="#ff000050", border="darkred")
nes3 <- (3*sd(auc)) + mean(auc)
abline(v=nes3, col="red")

TF_AtroCtrl_RcisTarget_bg_mod <- {
  allTFs <- RcisTarget_AtroCtrl_motifEnrichmentTable_bg$TF_highConf
  genes <- gsub(" \\(.*\\). ", "; ", allTFs, fixed=FALSE) 
  genesSplit <- unique(unlist(strsplit(genes, "; ")))
  genesFormatted <- sapply(genesSplit, function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  })
  genesFormatted
}

gene_sets_TF <- list(
  ChEA3 = TF_AtroCtrl_ChEA3_mod,
  BART = TF_AtroCtrl_BART_mod,
  # RcisTarget = TF_AtroCtrl_RcisTarget_mod,
  RcisTarget_bg = TF_AtroCtrl_RcisTarget_bg_mod
  # DEGs = DEG_GeneSymbol_AtroCtrl
)

ggvenn(gene_sets_TF, 
       fill_color = viridis(3),
       stroke_size = 0.5, 
       set_name_size = 4, 
       text_size = 3)

Reduce(intersect, gene_sets_TF)

## ChEA3 + BRAT + RcisTarget TF and your DEGs

## filter DEGs for fold change

DEG_TF_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl %>%
  filter(Gene_Symbol %in% c(TF_AtroCtrl_ChEA3_mod, TF_AtroCtrl_BART_mod, TF_AtroCtrl_RcisTarget_bg_mod)) %>%
  mutate(TF_pred_method = case_when(
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod & Gene_Symbol %in% TF_AtroCtrl_BART_mod & Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "ChEA3_BART_RcisTargetbg",
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod & Gene_Symbol %in% TF_AtroCtrl_BART_mod ~ "ChEA3_BART",
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod & Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "ChEA3_RcisTargetbg",
    Gene_Symbol %in% TF_AtroCtrl_BART_mod & Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "BART_RcisTargetbg",
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod ~ "ChEA3", 
    Gene_Symbol %in% TF_AtroCtrl_BART_mod ~ "BART",
    Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "RcisTarget_bg",
    TRUE ~ "Other"))

## filter expression matrix for TPM

expression_data_TF_AtroCtrl <- expression_data_AtroCtrl %>% 
  filter(Gene_Symbol %in% c(TF_AtroCtrl_ChEA3_mod, TF_AtroCtrl_BART_mod, TF_AtroCtrl_RcisTarget_bg_mod)) %>%
  mutate(TF_pred_method = case_when(
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod & Gene_Symbol %in% TF_AtroCtrl_BART_mod & Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "ChEA3_BART_RcisTargetbg",
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod & Gene_Symbol %in% TF_AtroCtrl_BART_mod ~ "ChEA3_BART",
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod & Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "ChEA3_RcisTargetbg",
    Gene_Symbol %in% TF_AtroCtrl_BART_mod & Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "BART_RcisTargetbg",
    Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod ~ "ChEA3", 
    Gene_Symbol %in% TF_AtroCtrl_BART_mod ~ "BART",
    Gene_Symbol %in% TF_AtroCtrl_RcisTarget_bg_mod ~ "RcisTarget_bg",
    TRUE ~ "Other"))

sd_avg_TPM <- sd(expression_data_TF_AtroCtrl$treatment1_avg_TPM)

ggplot()+
  # geom_rect(aes(xmin = mean(expression_data_AtroCtrl$treatment1_avg_TPM) - sd(expression_data_AtroCtrl$treatment1_avg_TPM), 
  #               xmax = mean(expression_data_AtroCtrl$treatment1_avg_TPM) + sd(expression_data_AtroCtrl$treatment1_avg_TPM), 
  #               ymin = -Inf, ymax = Inf), 
  #           fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = (mean(expression_data_AtroCtrl$treatment1_avg_TPM)), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(expression_data_TF_AtroCtrl,
                         TF_pred_method == "ChEA3_BART_RcisTargetbg" | 
                           TF_pred_method == "ChEA3_BART" |
                           TF_pred_method == "ChEA3_RcisTargetbg" | 
                           TF_pred_method == "BART_RcisTargetbg"),
             mapping=aes(x=(treatment1_avg_TPM),y=reorder(Gene_Symbol,treatment1_avg_TPM),col=treatment1_avg_TPM),size=3)+
  # geom_errorbar(data=subset(expression_data_TF_AtroCtrl, TF_pred_method == "ChEA3_BART"), mapping=aes(xmin=treatment1_avg_TPM-sd_avg_TPM, xmax=treatment1_avg_TPM+sd_avg_TPM,y=Gene_Symbol,col=treatment1_avg_TPM), width=0.3)+
  # scale_color_manual(values = c("firebrick"))+
  scale_x_continuous(limits = c(0, NA))+
  scale_color_viridis()+
  labs(y="TF",x="Expression level",color="")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))+
  scale_y_discrete(labels = function(x){ifelse(x == "Cebpb", expression(bold("Cebpb")), x)})

# expression_data_TF_AtroCtrl <- expression_data_AtroCtrl %>%
#   filter(Gene_Symbol %in% TF_AtroCtrl_ChEA3_mod)
# DEGexp_TF_AtroCtrl <- cbind(DEG_TF_AtroCtrl, expression_data_TF_AtroCtrl) %>%
#   dplyr::select(-Gene_Symbol)
# DEGexp_TF_AtroCtrl$Gene_Symbol <- rownames(DEGexp_TF_AtroCtrl)
# sd_log2FoldChange <- sd(DEGexp_TF_AtroCtrl$log2FoldChange)

sd_avg_TPM <- sd(DEGexp_TF_AtroCtrl$treatment1_avg_TPM)

ggplot()+
  geom_rect(aes(xmin = mean(DEGexp_TF_AtroCtrl$treatment1_avg_TPM) - sd(DEGexp_TF_AtroCtrl$treatment1_avg_TPM), 
                xmax = mean(DEGexp_TF_AtroCtrl$treatment1_avg_TPM) + sd(DEGexp_TF_AtroCtrl$treatment1_avg_TPM), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(DEGexp_TF_AtroCtrl$treatment1_avg_TPM), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(DEGexp_TF_AtroCtrl, direction == "upregulated"), mapping=aes(x=treatment1_avg_TPM,y=reorder(Gene_Symbol,treatment1_avg_TPM),col=log2FoldChange))+
  geom_errorbar(data=subset(DEGexp_TF_AtroCtrl, direction == "upregulated"), mapping=aes(xmin=treatment1_avg_TPM-sd_avg_TPM, xmax=treatment1_avg_TPM+sd_avg_TPM,y=Gene_Symbol,col=log2FoldChange), width=0.3)+
  # scale_color_manual(values = c("firebrick"))+
  scale_color_viridis()+
  labs(y="TF",x="Atrophy expression level (TPM)",color="Fold change")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# Final - only Atrophy - lncRNA -------------------------------------------

## long non coding RNA - maybe?

DEG_lncRNA_AtroCtrl <- DEG_full_info_df_1comp_AtroCtrl[grepl("^Gm\\d+", DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol) |
                                                         grepl("^LOC\\d+", DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol) |
                                                         grepl("\\d+Rik$", DEG_full_info_df_1comp_AtroCtrl$Gene_Symbol),]
nrow(DEG_full_info_df_1comp_AtroCtrl)
nrow(DEG_lncRNA_AtroCtrl)

gene_sets_lncRNA_AtroCtrl <- DEG_lncRNA_AtroCtrl %>%
  group_by(group) %>%
  summarize(Gene_Symbols = list(unique(Gene_Symbol)), .groups = "drop") %>%
  deframe()
lncRNA_AtroCtrl <- as.vector(gene_sets_lncRNA_AtroCtrl$`Atrophy vs Control`)
expression_data_lncRNA_AtroCtrl <- expression_data_AtroCtrl %>% 
  filter(Gene_Symbol %in% lncRNA_AtroCtrl)
DEGexp_lncRNA_AtroCtrl <- cbind(DEG_lncRNA_AtroCtrl, expression_data_lncRNA_AtroCtrl) %>%
  dplyr::select(-Gene_Symbol)
DEGexp_lncRNA_AtroCtrl$Gene_Symbol <- rownames(DEGexp_lncRNA_AtroCtrl)
sd_avg_TPM_lnc <- sd(DEGexp_lncRNA_AtroCtrl$treatment1_avg_TPM)
sd_log2FoldChange_lnc <- sd(DEGexp_lncRNA_AtroCtrl$log2FoldChange)

ggplot()+
  geom_rect(aes(xmin = mean(DEGexp_TF_AtroCtrl$treatment1_avg_TPM) - sd(DEGexp_TF_AtroCtrl$treatment1_avg_TPM), 
                xmax = mean(DEGexp_TF_AtroCtrl$treatment1_avg_TPM) + sd(DEGexp_TF_AtroCtrl$treatment1_avg_TPM), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(DEGexp_TF_AtroCtrl$treatment1_avg_TPM), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(DEGexp_lncRNA_AtroCtrl), mapping=aes(x=treatment1_avg_TPM,y=reorder(Gene_Symbol,treatment1_avg_TPM),col=log2FoldChange))+
  geom_errorbar(data=subset(DEGexp_lncRNA_AtroCtrl), mapping=aes(xmin=treatment1_avg_TPM-sd_avg_TPM_lnc, xmax=treatment1_avg_TPM+sd_avg_TPM_lnc,y=Gene_Symbol,col=log2FoldChange), width=0.3)+
  scale_color_viridis()+
  labs(y="TF",x="Atrophy expression level (TPM)",color="Fold change")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

ggplot()+
  geom_rect(aes(xmin = mean(DEGexp_TF_AtroCtrl$log2FoldChange) - sd(DEGexp_TF_AtroCtrl$log2FoldChange), 
                xmax = mean(DEGexp_TF_AtroCtrl$log2FoldChange) + sd(DEGexp_TF_AtroCtrl$log2FoldChange), 
                ymin = -Inf, ymax = Inf), 
            fill = "gray", alpha = 0.8) + 
  geom_vline(xintercept = mean(DEGexp_TF_AtroCtrl$log2FoldChange), linetype = "solid", color = "black", size = 1)+
  geom_point(data=subset(DEGexp_lncRNA_AtroCtrl), mapping=aes(x=log2FoldChange,y=reorder(Gene_Symbol,log2FoldChange),col=direction))+
  geom_errorbar(data=subset(DEGexp_lncRNA_AtroCtrl), mapping=aes(xmin=log2FoldChange-sd_log2FoldChange_lnc, xmax=log2FoldChange+sd_log2FoldChange_lnc,y=Gene_Symbol,col=direction), width=0.3)+
  scale_fill_manual(values = c("#FFD662FF", "#00539CFF")) +
  labs(y="TF",x="Atrophy expression level (TPM)",color="Fold change")+
  theme_linedraw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),)+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

gene_sets_lncRNA_AtroCtrl$group <- as.character(gene_sets_lncRNA_AtroCtrl$group)
deg_sets_AtroCtrl_lncRNA <- list(gene_sets_lncRNA_AtroCtrl$`Atrophy vs Control`)
names(deg_sets_AtroCtrl_lncRNA) <- "Atrophy vs Control"
ORA_GO_lncRNA_level_AtroCtrl <- lapply(names(deg_sets_AtroCtrl_lncRNA), function(treatment) {
  perform_ORA_GO(deg_sets_AtroCtrl_lncRNA[[treatment]], treatment)
})
ORA_GO_lncRNA_level_AtroCtrl_df <- bind_rows(ORA_GO_lncRNA_level_AtroCtrl)
ORA_GO_lncRNA_level_AtroCtrl_df$level <- "lncRNA"

# Supp_tables -------------------------------------------------------------

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/transcriptome_alzghoul_SS"

## convert all Supp_tables to proper data.frame class, if needed

Raw_Trim_Map_data <- as.data.frame(Raw_Trim_Map_data)
expression_data <- as.data.frame(expression_data)
AEG_full_info_df_1comp_AtroCtrl <- as.data.frame(AEG_full_info_df_1comp_AtroCtrl)
DEG_full_info_df_1comp_AtroCtrl <- as.data.frame(DEG_full_info_df_1comp_AtroCtrl)
DEG_sum_AtroCtrl_long <- as.data.frame(DEG_sum_AtroCtrl_long)
ORA_GO_alldeg_level_AtroCtrl_df <- as.data.frame(ORA_GO_alldeg_level_AtroCtrl_df)
ORA_GO_AtroCtrl_top10 <- as.data.frame(ORA_GO_AtroCtrl_top10)
ORA_KEGG_AtroCtrl_top <- as.data.frame(ORA_KEGG_AtroCtrl_top)
ORA_KEGG_alldeg_level_AtroCtrl_df <- as.data.frame(ORA_KEGG_alldeg_level_AtroCtrl_df)
pathway_genes_KEGG_AtroCtrl_clean_colored <- as.data.frame(pathway_genes_KEGG_AtroCtrl_clean_colored)
GSEA_GO_alldeg_level_AtroCtrl <- as.data.frame(GSEA_GO_alldeg_level_AtroCtrl@result)
GSEA_KEGG_alldeg_level_AtroCtrl <- as.data.frame(GSEA_KEGG_alldeg_level_AtroCtrl@result)
RcisTarget_AtroCtrl_motifEnrichmentTable_bg <- as.data.frame(RcisTarget_AtroCtrl_motifEnrichmentTable_bg)
ChEA3_Integrated_meanRank_DEG_AtroCtrl <- as.data.frame(ChEA3_Integrated_meanRank_DEG_AtroCtrl)
DEG_Gene_Symbol_AtroCtrl_bart_results <- as.data.frame(DEG_Gene_Symbol_AtroCtrl_bart_results)
expression_data_TF_AtroCtrl <- as.data.frame(expression_data_TF_AtroCtrl)
DEG_TF_AtroCtrl <- as.data.frame(DEG_TF_AtroCtrl)

## create the excel file

Rdataframes_titles <- data.frame(
  
  Sheet = c("Sheet1", "Sheet2", "Sheet3", "Sheet4", "Sheet5", "Sheet6", "Sheet7", "Sheet8", "Sheet9", "Sheet10",
            "Sheet11", "Sheet12", "Sheet13", "Sheet14", "Sheet15", "Sheet16", "Sheet17"
  ),
  
  DataFrameName = c("Raw_Trim_Map_data",
                    "expression_data",
                    "AEG_full_info_df_1comp_AtroCtrl",
                    "DEG_full_info_df_1comp_AtroCtrl",
                    "DEG_sum_AtroCtrl_long",
                    
                    "ORA_GO_alldeg_level_AtroCtrl_df",
                    "ORA_GO_AtroCtrl_top10",
                    "ORA_KEGG_alldeg_level_AtroCtrl_df",
                    "ORA_KEGG_AtroCtrl_top",
                    "pathway_genes_KEGG_AtroCtrl_clean_colored",
                    "GSEA_GO_alldeg_level_AtroCtrl",
                    "GSEA_KEGG_alldeg_level_AtroCtrl",
                    
                    "RcisTarget_AtroCtrl_motifEnrichmentTable_bg",
                    "ChEA3_Integrated_meanRank_DEG_AtroCtrl",
                    "DEG_Gene_Symbol_AtroCtrl_bart_results",
                    "expression_data_TF_AtroCtrl",
                    "DEG_TF_AtroCtrl"
  ),
  
  DataFrameNote = c("Summary statistics info for raw, trimmed, and mapped RNA-Seq data - Macrogen - AtroCtrl_comp",
                    "Expression Profile by genes (mm10) - Macrogen",
                    "All Expressed Genes for atrophy control comparison",
                    "Differentially Expressed Genes for atrophy control comparison",
                    "DEGs counts in atrophy control comparison (AtroCtrl_comp)",
                    
                    "ORA for GO terms results (Biological Process, Molecular Function, Cellular Component) - AtroCtrl_comp",
                    "Top ten scoring ORA for GO terms results - AtroCtrl_comp",
                    "ORA for KEGG pathways results - AtroCtrl_comp",
                    "Top ten scoring ORA for KEGG pathways results - AtroCtrl_comp",
                    "Top five scoring ORA for KEGG pathways with their corresponding genes and fold change - AtroCtrl_comp",
                    "GSEA for GO terms results (Biological Process, Molecular Function, Cellular Component) - AtroCtrl_comp",
                    "GSEA for KEGG pathways results - AtroCtrl_comp",
                    
                    "TF prediction method 1 - The motif enrichment table from RcisTarget background - AtroCtrl_comp",
                    "TF prediction method 2 - The integrated meanRank results from ChEA3 - AtroCtrl_comp",
                    "TF prediction method 3 - The results table from BART - AtroCtrl_comp",
                    "The top-scoring TF from the three methods (RcisTarget, ChEA3, BART) with their expression levels (TPM) - AtroCtrl_comp",
                    "The top-scoring TF from the three methods (RcisTarget, ChEA3, BART) that are part of the DEGs set - AtroCtrl_comp"
  )
)

write_xlsx(list(
  
  ContentTable = Rdataframes_titles,
  
  Sheet1 = Raw_Trim_Map_data,
  Sheet2 = expression_data,
  Sheet3 = AEG_full_info_df_1comp_AtroCtrl,
  Sheet4 = DEG_full_info_df_1comp_AtroCtrl,
  Sheet5 = DEG_sum_AtroCtrl_long,
  
  Sheet6 = ORA_GO_alldeg_level_AtroCtrl_df,
  Sheet7 = ORA_GO_AtroCtrl_top10,
  Sheet8 = ORA_KEGG_AtroCtrl_top,
  Sheet9 = ORA_KEGG_alldeg_level_AtroCtrl_df,
  Sheet10 = pathway_genes_KEGG_AtroCtrl_clean_colored,
  Sheet11 = GSEA_GO_alldeg_level_AtroCtrl,
  Sheet12 = GSEA_KEGG_alldeg_level_AtroCtrl,
  
  Sheet13 = RcisTarget_AtroCtrl_motifEnrichmentTable_bg,
  Sheet14 = ChEA3_Integrated_meanRank_DEG_AtroCtrl,
  Sheet15 = DEG_Gene_Symbol_AtroCtrl_bart_results,
  Sheet16 = expression_data_TF_AtroCtrl,
  Sheet17 = DEG_TF_AtroCtrl),
  
  file.path(sys_dir,"supp_tables/supp_tables_Rdataframes.xlsx"))


