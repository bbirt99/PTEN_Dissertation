####---- master script ----####
#packages#
library(ggplot2)
library(DESeq2)
library(tidyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggrepel)
library(cowplot)
library(dplyr)
library(limma)
library(magick)
library(ggprism)
library(patchwork)
library(missForest)
library(proDA)
library(tibble)
library(biomaRt)
library(clusterProfiler)
library(MsCoreUtils)
library(stringr)
library(MOFA2)
library(basilisk)
library(reshape2)
library(fgsea)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(GSVA)
library(progeny)
library(Mfuzz)
library(Biobase)
library(ReactomePA)
library(ggplotify)
library(gridGraphics)
library(RcisTarget)
library(decoupleR)
library(purrr)
library(OmnipathR)
library(CARNIVAL)
library(gurobi)
library(slam)
library(igraph)
library(ggrastr)
library(ggsignif)
library(dorothea)
library(readr)
library(ggpubr)
library(ppcor)
library(Seurat)
library(harmony)
library(DoubletFinder)
library(dbscan)

#### DESeq2 Object Construction and Results Generation ####

set.seed(123)

counts_unfiltered <- read.csv("RNAseq_1/all_raw_count_preprocess.csv", row.names = 1, check.names = FALSE)
coldata <- read.csv("RNAseq_1/coldata_ext.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = counts_unfiltered,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

dds <- dds[, !colnames(dds) %in% c("PTEN-1", "PTEN-5")]

dds_wald <- DESeq(dds)

res_2w_24h <- results(dds_wald, contrast=c("condition", "2w", "24h"))
res_24h_0h <- results(dds_wald, contrast=c("condition", "24h", "0h"))
res_2w_0h  <- results(dds_wald, contrast=c("condition", "2w", "0h"))

df_2w_24h_whole <- as.data.frame(res_2w_24h)

df_2w_24h <- as.data.frame(res_2w_24h)[, c("log2FoldChange", "padj")]
colnames(df_2w_24h) <- c("LFC_2w_vs_24h", "Padj_2w_vs_24h")

df_24h_0h <- as.data.frame(res_24h_0h)[, c("log2FoldChange", "padj")]
colnames(df_24h_0h) <- c("LFC_24h_vs_0h", "Padj_24h_vs_0h")

df_2w_0h <- as.data.frame(res_2w_0h)[, c("log2FoldChange", "padj")]
colnames(df_2w_0h) <- c("LFC_2w_vs_0h", "Padj_2w_vs_0h")

dds_lrt <- DESeq(dds, test="LRT", reduced=~1)
res_lrt <- results(dds_lrt)

df_lrt <- as.data.frame(res_lrt)[, "padj", drop=FALSE]
colnames(df_lrt) <- "Padj_LRT_Global"

master_df <- cbind(df_24h_0h, df_2w_24h, df_2w_0h, df_lrt)
master_df$Ensembl <- rownames(master_df)

master_df$Symbol <- mapIds(org.Hs.eg.db,
                           keys = as.character(master_df$Ensembl), 
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")

master_df$Entrez <- mapIds(org.Hs.eg.db,
                           keys = as.character(master_df$Ensembl),
                           column = "ENTREZID",
                           keytype = "ENSEMBL",
                           multiVals = "first")

vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)

vst_symbols <- mapIds(org.Hs.eg.db,
                      keys = rownames(vst_matrix),
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

valid_genes <- !is.na(vst_symbols)
vst_matrix_clean <- vst_matrix[valid_genes, ]
vst_symbols_clean <- vst_symbols[valid_genes]

vst_mat_unique <- avereps(vst_matrix_clean, ID = vst_symbols_clean)

chronological_order <- c("PTEN-2", "PTEN-3", "PTEN-4", 
                         "PTEN-6", "PTEN-7", "PTEN-8", 
                         "PTEN-9", "PTEN-10", "PTEN-11", "PTEN-12")

vst_mat_unique <- vst_mat_unique[, chronological_order]

#### LFQ Value Processing and proDA Analysis ####

df <- read.delim("proteomics/proteinGroups.txt", sep = "\t", header = TRUE, 
                 stringsAsFactors = FALSE, check.names = FALSE)

df <- df[df$`Potential contaminant` != "+" & 
           df$Reverse != "+" & 
           df$`Unique peptides` >= 2, ]

lfq_cols <- grep("^LFQ", colnames(df), value = TRUE)
lfq <- df[, lfq_cols]
lfq[lfq == 0] <- NA
rownames(lfq) <- df$`Protein IDs`

group <- factor(rep(c("T_0h", "T_24h", "T_2w"), each = 4),
                levels = c("T_0h", "T_24h", "T_2w"))
groups_list <- split(colnames(lfq), group)

lfq_log2 <- log2(lfq)

leading_ids <- sapply(strsplit(rownames(lfq_log2), ";"), `[`, 1)
clean_uniprot <- gsub("^.*\\|([A-Z0-9-]+)\\|.*$", "\\1", leading_ids)

mapping <- bitr(clean_uniprot, 
                fromType = "UNIPROT", 
                toType = "SYMBOL", 
                OrgDb = org.Hs.eg.db)

colnames(mapping) <- c("uniprotswissprot", "hgnc_symbol")

prot_df <- as.data.frame(lfq_log2)
prot_df$uniprotswissprot <- clean_uniprot
prot_merged <- merge(mapping, prot_df, by = "uniprotswissprot")

prot_mat_symbols <- avereps(prot_merged[, -c(1, 2)], ID = prot_merged$hgnc_symbol)

keep_proda <- apply(prot_mat_symbols, 1, function(x) {
  any(sapply(groups_list, function(cols) sum(!is.na(x[cols])) >= 3))
})

lfq_for_proda <- prot_mat_symbols[keep_proda, , drop = FALSE] 

fit <- proDA(lfq_for_proda, design = ~ 0 + group)

res_24hr_vs_0h <- test_diff(fit, "groupT_24h - groupT_0h") %>%
  arrange(adj_pval)

res_2w_vs_0h <- test_diff(fit, "groupT_2w - groupT_0h") %>%
  arrange(adj_pval)

res_2w_vs_24hr <- test_diff(fit, "groupT_2w - groupT_24h") %>%
  arrange(adj_pval)

prot_master_df <- res_2w_vs_0h %>%
  dplyr::select(name, logFC_2w_0h = diff, padj_2w_0h = adj_pval) %>%
  left_join(
    dplyr::select(res_24hr_vs_0h, name, logFC_24h_0h = diff, padj_24h_0h = adj_pval),
    by = "name"
  ) %>%
  left_join(
    dplyr::select(res_2w_vs_24hr, name, logFC_2w_24h = diff, padj_2w_24h = adj_pval),
    by = "name"
  )

res_lrt <- test_diff(fit, reduced_model = ~ 1)

lrt_subset <- res_lrt %>%
  dplyr::select(name, 
                pval_LRT = pval, 
                padj_LRT = adj_pval)

prot_master_df <- prot_master_df %>%
  left_join(lrt_subset, by = "name") %>%
  arrange(padj_LRT)

proda_mat <- as.matrix(lfq_for_proda)

lfq_mar_imputed <- impute_matrix(proda_mat, method = "knn")

set.seed(123)

global_mean <- mean(proda_mat, na.rm = TRUE)
global_sd   <- sd(proda_mat, na.rm = TRUE)

impute_mean <- global_mean - (1.8 * global_sd)
impute_sd   <- 0.3 * global_sd

noise_matrix <- matrix(rnorm(prod(dim(proda_mat)), mean = impute_mean, sd = impute_sd),
                       nrow = nrow(proda_mat), ncol = ncol(proda_mat))

final_imputed_matrix <- proda_mat

for(i in 1:nrow(final_imputed_matrix)) {
  
  is_mnar <- any(sapply(groups_list, function(cols) sum(is.na(final_imputed_matrix[i, cols])) == length(cols)))
  
  nas <- is.na(final_imputed_matrix[i, ])
  
  if(any(nas)) {
    if(is_mnar) {
      final_imputed_matrix[i, nas] <- noise_matrix[i, nas]
    } else {
      final_imputed_matrix[i, nas] <- lfq_mar_imputed[i, nas]
    }
  }
}

rownames(final_imputed_matrix) <- rownames(proda_mat)
colnames(final_imputed_matrix) <- colnames(proda_mat)

#### PCA Generation ####
shared_colors <- c("0h" = "#FDE725", "24h" = "#21918C", "2w" = "#440154")

my_prism_theme <- theme_prism(base_family = "sans", base_size = 14) +
  theme(
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )

vsd_mat <- as.matrix(vst_mat_unique)

mrna_vars <- apply(vsd_mat, 1, var)

select_mrna <- order(mrna_vars, decreasing = TRUE)[1:500]

mrna_top500 <- vsd_mat[select_mrna, ]

pca_res_mrna <- prcomp(t(mrna_top500), scale. = FALSE)

var_explained_mrna <- pca_res_mrna$sdev^2 / sum(pca_res_mrna$sdev^2)
pc1_var_mrna <- round(100 * var_explained_mrna)[1]
pc2_var_mrna <- round(100 * var_explained_mrna)[2]

plot_data_mrna <- as.data.frame(pca_res_mrna$x)

plot_data_mrna$Timepoint <- factor(c(rep("0h", 3), rep("24h", 3), rep("2w", 4)), 
                                   levels = c("0h", "24h", "2w"))

max_val <- max(abs(c(plot_data_mrna$PC1, plot_data_mrna$PC2)))

lim <- max_val * 1.1

lim <- ceiling(lim)

pca_plot_mrna <- ggplot(plot_data_mrna, aes(x = PC1, y = PC2, color = Timepoint)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.5) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(values = shared_colors) +
  labs(
    title = "Transcriptomics", 
    x = paste0("PC1 (", pc1_var_mrna, "%)"),
    y = paste0("PC2 (", pc2_var_mrna, "%)"),
    color = "Time point"
  ) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  my_prism_theme

prot_vars <- apply(final_imputed_matrix, 1, var, na.rm = TRUE)

select_prot <- order(prot_vars, decreasing = TRUE)[1:500]

prot_top500 <- final_imputed_matrix[select_prot, ]

pca_res_prot <- prcomp(t(prot_top500), scale. = TRUE)

var_explained_prot <- pca_res_prot$sdev^2 / sum(pca_res_prot$sdev^2)
pc1_var_prot <- round(100 * var_explained_prot)[1]
pc2_var_prot <- round(100 * var_explained_prot)[2]

plot_data_prot <- as.data.frame(pca_res_prot$x)

plot_data_prot$Timepoint <- factor(
  c(rep("0h", 4), rep("24h", 4), rep("2w", 4)), 
  levels = c("0h", "24h", "2w")
)

max_val_prot <- max(abs(c(plot_data_prot$PC1, plot_data_prot$PC2)))

lim_prot <- ceiling(max_val_prot * 1.1)

pca_plot_prot <- ggplot(plot_data_prot, aes(x = PC1, y = PC2, color = Timepoint)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.5) +
  geom_point(size = 4, alpha = 0.85) + 
  scale_color_manual(values = shared_colors) +
  scale_x_continuous(breaks = seq(-100, 100, 10)) +
  scale_y_continuous(breaks = seq(-100, 100, 10)) +
  labs(
    title = "Proteomics",
    x = paste0("PC1 (", pc1_var_prot, "%)"),
    y = paste0("PC2 (", pc2_var_prot, "%)"),
    color = "Time point"
  ) +
  coord_fixed(xlim = c(-lim_prot, lim_prot), ylim = c(-lim_prot, lim_prot)) +
  my_prism_theme

combined_pca_figure <- pca_plot_mrna + pca_plot_prot

final_figure <- combined_pca_figure + 
  plot_layout(guides = "collect")
#### MOFA+ Top 3000 Feature Selection, Model Training, and Master Panel Assembly ####

vsd_expr <- as.matrix(vst_mat_unique)
prot_expr <- lfq_for_proda

colnames(vsd_expr) <- c(
  "Sample_0h_Rep1", "Sample_0h_Rep2", "Sample_0h_Rep3",          
  "Sample_24h_Rep1", "Sample_24h_Rep2", "Sample_24h_Rep3",       
  "Sample_2w_Rep1", "Sample_2w_Rep2", "Sample_2w_Rep3", "Sample_2w_Rep4" 
)

colnames(prot_expr) <- c(
  "Sample_0h_Rep1", "Sample_0h_Rep2", "Sample_0h_Rep3", "Sample_0h_Rep4",
  "Sample_24h_Rep1", "Sample_24h_Rep2", "Sample_24h_Rep3", "Sample_24h_Rep4",
  "Sample_2w_Rep1", "Sample_2w_Rep2", "Sample_2w_Rep3", "Sample_2w_Rep4"
)

target_features <- 3000

gene_var <- apply(vsd_expr, 1, var)
top_rna_features <- names(sort(gene_var, decreasing = TRUE))[1:target_features]

protein_var <- apply(prot_expr, 1, var, na.rm = TRUE)
n_prot <- min(target_features, length(protein_var))
top_protein_features <- names(sort(protein_var, decreasing = TRUE))[1:n_prot]

cat("RNA features selected:", length(top_rna_features), "\n")
cat("Protein features selected:", length(top_protein_features), "\n")

vsd_expr_top3k <- vsd_expr[top_rna_features, ]
prot_expr_top3k <- prot_expr[top_protein_features, ]

rna_long_top3k <- vsd_expr_top3k %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  pivot_longer(-feature, names_to = "sample", values_to = "value") %>%
  mutate(view = "RNA")

prot_long_top3k <- prot_expr_top3k %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  pivot_longer(-feature, names_to = "sample", values_to = "value") %>%
  mutate(view = "Protein")

mofa_df_top3k <- bind_rows(rna_long_top3k, prot_long_top3k)
MOFAobject_top3k <- create_mofa(mofa_df_top3k)

data_opts <- get_default_data_options(MOFAobject_top3k)
model_opts <- get_default_model_options(MOFAobject_top3k)
model_opts$num_factors <- 3  
model_opts$likelihoods <- c("RNA" = "gaussian", "Protein" = "gaussian") 

train_opts <- get_default_training_options(MOFAobject_top3k)
train_opts$convergence_mode <- "slow" 
train_opts$seed <- 123

MOFAobject_top3k <- prepare_mofa(
  object = MOFAobject_top3k,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

model_filepath <- "MDA_MB_468_PTEN_Resistance_top3k.hdf5"
MOFAobject_top3k <- run_mofa(MOFAobject_top3k, outfile = model_filepath, use_basilisk = TRUE)

MOFAobject_top3k <- load_model(model_filepath)

sample_names <- unname(unlist(samples_names(MOFAobject_top3k)))
metadata <- data.frame(sample = sample_names) %>%
  mutate(
    Timepoint = case_when(
      grepl("0h", sample) ~ "0h",
      grepl("24h", sample) ~ "24h",
      grepl("2w", sample) ~ "2w",
      TRUE ~ "Unknown"
    ),
    Timepoint = factor(Timepoint, levels = c("0h", "24h", "2w"))
  )
samples_metadata(MOFAobject_top3k) <- metadata

var_matrix <- get_variance_explained(MOFAobject_top3k)$r2_per_factor[]
var_df <- as.data.frame(var_matrix) %>%
  rownames_to_column("Factor") %>%
  pivot_longer(cols = -Factor, names_to = "Omics_Layer", values_to = "Variance_Explained") %>%
  filter(Factor %in% c("Factor1", "Factor2"))

var_plot <- ggplot(var_df, aes(x = factor(Omics_Layer, levels = c("single_group.RNA", "single_group.Protein")), 
                               y = Variance_Explained, fill = Omics_Layer)) +
  geom_bar(stat = "identity", color = "black", linewidth = 1, position = position_dodge()) + 
  facet_wrap(~Factor, strip.position = "bottom") +
  scale_fill_manual(values = c("single_group.RNA" = "#E69F00", "single_group.Protein" = "#228B22")) +
  scale_x_discrete(labels = c("single_group.RNA" = "Transcriptomics", "single_group.Protein" = "Proteomics")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Variance Explained (%)", x = NULL, title = "Variance Explained") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none", 
    axis.line = element_line(color = "black", linewidth = 1.2),
    axis.text.x = element_text(face = "bold", angle = 45, hjust = 1, vjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14)
  )

factor_plot <- plot_factor(MOFAobject_top3k, factors = c(1, 2), color_by = "Timepoint", dot_size = 4) +
  facet_wrap(~factor, strip.position = "bottom", nrow = 1) +
  theme_classic(base_size = 14) +
  scale_fill_manual(values = c("0h" = "#FDE725FF", "24h" = "#21908CFF", "2w" = "#440154FF")) +
  scale_color_manual(values = c("0h" = "black", "24h" = "black", "2w" = "black")) +
  labs(title = "Latent Factor Projection", x = NULL) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0, "pt"),
    axis.line.x.bottom = element_line(color = "black", linewidth = 1),
    axis.line.y = element_line(color = "black", linewidth = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x.bottom = element_text(face = "bold", size = 12, margin = margin(t = 12, b = 5)),
    plot.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(
      title = "Timepoint", 
      override.aes = list(shape = 21, size = 5, fill = c("#FDE725FF", "#21908CFF", "#440154FF"), color = "black")
    ), 
    color = "none"
  )

plot_mrna_weights <- plot_top_weights(MOFAobject_top3k, view = "RNA", factor = 1, nfeatures = 20) 

levels(plot_mrna_weights$data$feature_id) <- stringr::str_remove(levels(plot_mrna_weights$data$feature_id), "_RNA")

plot_mrna_weights <- plot_mrna_weights +
  facet_wrap(~factor, labeller = as_labeller(c("Factor1" = ""))) +
  labs(title = "Top 20 Transcripts (Factor 1)") +
  theme(plot.title = element_text(face = "bold", size = 14))

plot_prot_weights <- plot_top_weights(MOFAobject_top3k, view = "Protein", factor = 1, nfeatures = 20, sign = "all")

levels(plot_prot_weights$data$feature_id) <- stringr::str_remove(levels(plot_prot_weights$data$feature_id), "_Protein")

plot_prot_weights <- plot_prot_weights +
  facet_wrap(~factor, labeller = as_labeller(c("Factor1" = ""))) +
  labs(title = "Top 20 Proteins (Factor 1)") +
  theme(plot.title = element_text(face = "bold", size = 14))

set.seed(123)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(x = m_df$gene_symbol, f = m_df$gs_name)

weights_rna <- get_weights(MOFAobject_top3k, views = "RNA", factors = 1)[[1]][, 1]
weights_rna_sorted <- sort(weights_rna, decreasing = TRUE)
fgsea_rna <- fgsea(pathways = pathways, stats = weights_rna_sorted, minSize = 10, maxSize = 500)
top_rna <- fgsea_rna %>% filter(padj < 0.05) %>% arrange(desc(NES))

weights_prot <- get_weights(MOFAobject_top3k, views = "Protein", factors = 1)[[1]][, 1]
names(weights_prot) <- gsub("_Protein", "", names(weights_prot)) 
weights_prot_sorted <- sort(weights_prot, decreasing = TRUE)
fgsea_prot <- fgsea(pathways = pathways, stats = weights_prot_sorted, minSize = 10, maxSize = 500)
top_prot <- fgsea_prot %>% filter(padj < 0.05) %>% arrange(desc(NES))

gsea_merged <- inner_join(top_rna, top_prot, by = "pathway", suffix = c("_RNA", "_Prot"))
gsea_plot_data <- gsea_merged %>%
  mutate(pathway_clean = str_replace(pathway, "HALLMARK_", "")) %>%
  arrange(NES_RNA) %>% 
  mutate(pathway_clean = factor(pathway_clean, levels = pathway_clean),
         is_overlap = abs(NES_RNA - NES_Prot) < 0.01)

dumbbell_plot <- ggplot(gsea_plot_data, aes(y = pathway_clean)) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_segment(aes(x = NES_RNA, xend = NES_Prot, yend = pathway_clean), color = "gray50", linewidth = 1.5, alpha = 0.5) +
  geom_point(aes(x = NES_Prot, fill = "Proteomics"), size = 6, shape = 21, color = "black", stroke = 1, alpha = 0.85) +
  geom_point(data = subset(gsea_plot_data, !is_overlap), aes(x = NES_RNA, fill = "Transcriptomics"), size = 6, shape = 21, color = "black", stroke = 1, alpha = 0.85) +
  geom_point(data = subset(gsea_plot_data, is_overlap), aes(x = NES_RNA, fill = "Transcriptomics"), size = 3, shape = 21, color = "black", stroke = 1, alpha = 1) +
  scale_fill_manual(values = c("Transcriptomics" = "#E69F00", "Proteomics" = "#228B22")) +
  theme_prism(base_size = 14, base_family = "sans") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Modality", title = "Pathway Concordance (Factor 1)") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold", size = 14))

isolated_dumbbell <- wrap_elements(panel = dumbbell_plot)

master_figure <- (var_plot | factor_plot) / 
  (plot_mrna_weights | plot_prot_weights) / 
  (isolated_dumbbell) +
  plot_layout(heights = c(0.8, 1.2, 1.5)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))
#### Multi-Omics Volcano Plots (2x3 Grid) ####

generate_prism_rna_volcano <- function(data_frame, contrast_prefix, panel_title, x_lim = 6, y_lim = 10) {
  x_col <- paste0("LFC_", contrast_prefix)
  p_col <- paste0("Padj_", contrast_prefix)
  
  df_clean <- data_frame[!is.na(data_frame[[p_col]]), ]
  
  df_clean$status <- "Not Significant"
  df_clean$status[df_clean[[x_col]] > 1 & df_clean[[p_col]] < 0.05] <- "Upregulated"
  df_clean$status[df_clean[[x_col]] < -1 & df_clean[[p_col]] < 0.05] <- "Downregulated"
  
  df_labeled <- df_clean[!is.na(df_clean$Symbol) & df_clean$Symbol != "<NA>" & df_clean$Symbol != "", ]
  top_5 <- head(df_labeled[order(df_labeled[[p_col]]), ], 5)
  
  ggplot(df_clean, aes(x = .data[[x_col]], y = -log10(.data[[p_col]]))) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#7F7F7F", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#7F7F7F", linewidth = 0.5) +
    geom_point(aes(color = status), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Downregulated" = "#1F77B4", "Not Significant" = "#7F7F7F", "Upregulated" = "#D62728")) +
    geom_text_repel(data = top_5, aes(label = Symbol), size = 2.8, max.overlaps = 20, show.legend = FALSE, fontface = "italic") +
    scale_x_continuous(limits = c(-x_lim, x_lim), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_lim), expand = c(0, 0)) +
    labs(title = panel_title, x = expression(Log[2]~Fold~Change), y = expression(-Log[10](Adjusted~P-value))) +
    theme_prism(base_size = 11, border = TRUE) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_lim), expand = FALSE)
}

plot_rna_24h <- generate_prism_rna_volcano(master_df, "24h_vs_0h", "mRNA: 0h vs 24h", x_lim = 6, y_lim = 75)
plot_rna_2w  <- generate_prism_rna_volcano(master_df, "2w_vs_24h", "mRNA: 24h vs 2w", x_lim = 8, y_lim = 175)
plot_rna_net <- generate_prism_rna_volcano(master_df, "2w_vs_0h", "mRNA: 0h vs 2w", x_lim = 8, y_lim = 200)

generate_prism_prot_volcano <- function(data_frame, contrast_suffix, panel_title, x_lim = 4, y_lim = 8, lfc_thresh = 0.58) {
  x_col <- paste0("logFC_", contrast_suffix)
  p_col <- paste0("padj_", contrast_suffix)
  
  df_clean <- data_frame[!is.na(data_frame[[p_col]]), ]
  
  df_clean$status <- "Not Significant"
  df_clean$status[df_clean[[x_col]] > lfc_thresh & df_clean[[p_col]] < 0.05] <- "Upregulated"
  df_clean$status[df_clean[[x_col]] < -lfc_thresh & df_clean[[p_col]] < 0.05] <- "Downregulated"
  
  df_labeled <- df_clean[!is.na(df_clean$name) & df_clean$name != "<NA>" & df_clean$name != "", ]
  top_5 <- head(df_labeled[order(df_labeled[[p_col]]), ], 5)
  
  ggplot(df_clean, aes(x = .data[[x_col]], y = -log10(.data[[p_col]]))) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "#7F7F7F", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#7F7F7F", linewidth = 0.5) +
    geom_point(aes(color = status), alpha = 0.7, size = 1.2) +
    scale_color_manual(values = c("Downregulated" = "#1F77B4", "Not Significant" = "#7F7F7F", "Upregulated" = "#D62728")) +
    geom_text_repel(data = top_5, aes(label = name), size = 2.8, max.overlaps = 20, show.legend = FALSE, fontface = "italic") +
    scale_x_continuous(limits = c(-x_lim, x_lim), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, y_lim), expand = c(0, 0)) +
    labs(title = panel_title, x = expression(Log[2]~Fold~Change), y = expression(-Log[10](Adjusted~P-value))) +
    theme_prism(base_size = 11, border = TRUE) +
    theme(legend.position = "none", plot.title = element_text(size = 11, face = "bold", hjust = 0.5)) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_lim), expand = FALSE)
}

plot_prot_24h <- generate_prism_prot_volcano(prot_master_df, "24h_0h", "Protein: 0h vs 24h", x_lim = 4, y_lim = 4, lfc_thresh = 0.58)
plot_prot_2w  <- generate_prism_prot_volcano(prot_master_df, "2w_24h", "Protein: 24h vs 2w", x_lim = 5, y_lim = 7, lfc_thresh = 0.58)
plot_prot_net <- generate_prism_prot_volcano(prot_master_df, "2w_0h", "Protein: 0h vs 2w", x_lim = 5, y_lim = 7, lfc_thresh = 0.58)

full_multiomics_grid <- (plot_rna_24h + plot_prot_24h +
                           plot_rna_2w  + plot_prot_2w  +
                           plot_rna_net + plot_prot_net) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "A")
#### GSVA Multi-Omics Heatmaps ####

msigdb_h <- msigdbr(species = "Homo sapiens", collection  = "H")
hallmark <- split(msigdb_h$gene_symbol, msigdb_h$gs_name)

rna_param <- gsvaParam(
  exprData = vst_mat_unique,
  geneSets = hallmark,
  kcdf = "Gaussian",
  minSize = 10,
  maxSize = 500,
  tau = 1,
  absRanking = FALSE,
  maxDiff = TRUE
)
gsva_results_mrna <- gsva(rna_param)

prot_param <- gsvaParam(
  exprData = final_imputed_matrix,
  geneSets = hallmark,
  kcdf = "Gaussian",
  minSize = 10,
  maxSize = 500,
  tau = 1,
  absRanking = FALSE,
  maxDiff = TRUE
)
gsva_results_prot <- gsva(prot_param)

common_pathways <- intersect(rownames(gsva_results_mrna), rownames(gsva_results_prot))

timepoints_mrna <- factor(c(rep("0h", 3), rep("24h", 3), rep("2w", 4)), levels = c("0h", "24h", "2w"))
timepoints_prot <- factor(c(rep("0h", 4), rep("24h", 4), rep("2w", 4)), levels = c("0h", "24h", "2w"))

shared_colors <- list(Timepoint = c("0h" = "#FDE725FF", "24h" = "#21908CFF", "2w" = "#440154FF"))

col_anno_mrna = HeatmapAnnotation(
  Timepoint = timepoints_mrna,
  col = shared_colors,
  show_annotation_name = FALSE
)

col_anno_prot = HeatmapAnnotation(
  Timepoint = timepoints_prot,
  col = shared_colors,
  show_annotation_name = FALSE
)

col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

clean_pathway_names <- function(pathways) {
  clean_pw <- gsub("^HALLMARK_", "", pathways)
  clean_pw <- gsub("_", " ", clean_pw)
  clean_pw <- str_to_title(clean_pw)
  
  acronyms <- c("Pi3k", "Akt", "Myc", "E2f", "Jak", "Stat3", "Tnfa", "Nfkb", 
                "Uv", "Tgf", "Kras", "Il6", "Il2", "Stat5", "Dna", "G2m", "Wnt", "V1", "V2", "Up", "Dn")
  
  for (acr in acronyms) {
    clean_pw <- gsub(paste0("\\b", acr, "\\b"), toupper(acr), clean_pw)
  }
  
  clean_pw <- gsub("\\bMtorc1\\b", "mTORC1", clean_pw)
  clean_pw <- gsub("\\bMtor\\b", "mTOR", clean_pw)
  clean_pw <- gsub("\\bP53\\b", "p53", clean_pw)
  
  return(clean_pw)
}

mrna_sub_all <- gsva_results_mrna[common_pathways, ]
prot_sub_all <- gsva_results_prot[common_pathways, ]

clean_pathways_all <- clean_pathway_names(rownames(mrna_sub_all))

ht_mrna_all = Heatmap(
  t(scale(t(mrna_sub_all))), 
  name = "Enrichment Score", 
  column_title = "Transcriptomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  top_annotation = col_anno_mrna,
  cluster_columns = FALSE, 
  show_row_names = TRUE,
  row_names_side = "left",
  row_labels = clean_pathways_all,           
  row_names_gp = gpar(fontsize = 11),      
  width = unit(6, "cm"),                  
  show_column_names = FALSE,
  show_row_dend = FALSE                  
)

ht_prot_all = Heatmap(
  t(scale(t(prot_sub_all))), 
  name = "Enrichment Score",
  column_title = "Proteomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  top_annotation = col_anno_prot,
  cluster_columns = FALSE,
  width = unit(6, "cm"),                
  show_row_names = FALSE,                
  show_heatmap_legend = FALSE,           
  show_column_names = FALSE,  
  heatmap_legend_param = list(direction = "vertical")
)

final_layout_all = ht_mrna_all + ht_prot_all

draw(final_layout_all, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right")

#### PROGENy Pathway Activity Heatmaps ####

progeny_rna <- progeny(as.matrix(vst_mat_unique), scale = FALSE, organism = "Human", top = 100, perm = 1)
progeny_prot <- progeny(final_imputed_matrix, scale = TRUE, organism = "Human", top = 500, perm = 1)

mrna_sub_progeny <- t(progeny_rna)
prot_sub_progeny <- t(progeny_prot)

timepoints_mrna <- factor(c(rep("0h", 3), rep("24h", 3), rep("2w", 4)), levels = c("0h", "24h", "2w"))
timepoints_prot <- factor(c(rep("0h", 4), rep("24h", 4), rep("2w", 4)), levels = c("0h", "24h", "2w"))

shared_colors <- list(Timepoint = c("0h" = "#FDE725FF", "24h" = "#21908CFF", "2w" = "#440154FF"))

col_anno_mrna_progeny <- HeatmapAnnotation(
  Timepoint = timepoints_mrna,
  col = shared_colors,
  show_annotation_name = FALSE
)

col_anno_prot_progeny <- HeatmapAnnotation(
  Timepoint = timepoints_prot,
  col = shared_colors,
  show_annotation_name = FALSE
)

col_fun_progeny <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mrna_progeny <- Heatmap(
  mrna_sub_progeny,
  name = "Pathway Z-Score",
  column_title = "Transcriptomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun_progeny,
  top_annotation = col_anno_mrna_progeny,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11),
  width = unit(6, "cm"),
  height = unit(7, "cm"),
  show_column_names = FALSE,
  show_row_dend = FALSE
)

ht_prot_progeny <- Heatmap(
  prot_sub_progeny,
  name = "Pathway Z-Score",
  column_title = "Proteomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun_progeny,
  top_annotation = col_anno_prot_progeny,
  cluster_columns = FALSE,
  width = unit(6, "cm"),
  height = unit(7, "cm"),
  show_row_names = FALSE,
  show_heatmap_legend = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(direction = "vertical")
)

final_layout_progeny <- ht_mrna_progeny + ht_prot_progeny

draw(final_layout_progeny,
     merge_legends = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right"
)
#### Transcriptomic Mfuzz Clustering and ORA Analysis ####

set.seed(123)

expr <- vst_mat_unique
sig_rna <- master_df

lrt_mask <- !is.na(sig_rna$Padj_LRT_Global) & sig_rna$Padj_LRT_Global < 0.05
lrt_genes <- unique(sig_rna$Symbol[lrt_mask])
lrt_genes <- lrt_genes[!is.na(lrt_genes) & lrt_genes != ""]

valid_genes <- intersect(rownames(expr), lrt_genes)
expr_sig <- expr[valid_genes, ]

time_matrix <- matrix(
  c(
    rowMeans(expr_sig[, 1:3], na.rm = TRUE),
    rowMeans(expr_sig[, 4:6], na.rm = TRUE),
    rowMeans(expr_sig[, 7:10], na.rm = TRUE)
  ),
  ncol = 3,
  dimnames = list(rownames(expr_sig), c("T_0h", "T_24h", "T_336h"))
)

eset <- new("ExpressionSet", exprs = as.matrix(time_matrix))
eset_std <- standardise(eset)

cl_rna_m2 <- mfuzz(eset_std, c = 6, m = 2.00)

bg_symbols <- rownames(expr_sig)
bg_mapped <- bitr(bg_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez <- unique(bg_mapped$ENTREZID)

run_cluster_ora <- function(mfuzz_obj, cluster_num, membership_threshold = 0.5) {
  cluster_mem <- mfuzz_obj$membership[, cluster_num]
  target_symbols <- names(cluster_mem)[cluster_mem >= membership_threshold]
  
  mapped <- bitr(
    target_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  gene_entrez <- unique(mapped$ENTREZID)
  
  ora_res <- enrichPathway(
    gene         = gene_entrez,
    universe     = universe_entrez,
    organism     = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable     = TRUE
  )
  
  return(ora_res)
}

ora_list <- list()
cluster_labels <- c(
  "Cluster 1",
  "Cluster 2",
  "Cluster 3",
  "Cluster 4",
  "Cluster 5",
  "Cluster 6"
)

for (i in 1:6) {
  res <- run_cluster_ora(cl_rna_m2, cluster_num = i, membership_threshold = 0.5)
  df <- as.data.frame(res)
  if (nrow(df) > 0) {
    df_top <- df %>%
      slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
      mutate(
        Cluster = cluster_labels[i],
        GeneRatioNumeric = sapply(GeneRatio, function(x) {
          parts <- as.numeric(strsplit(x, "/")[[1]])
          parts[1] / parts[2]
        }),
        Description_wrapped = str_wrap(Description, width = 38)
      )
    ora_list[[i]] <- df_top
  }
}

combined_mfuzz_ora_all <- bind_rows(ora_list)

expr_matrix <- exprs(eset_std)
mem_matrix  <- cl_rna_m2$membership

target_clusters <- c(2, 4, 5)
cluster_names   <- c("2" = "Cluster 2", "4" = "Cluster 4", "5" = "Cluster 5")

traj_list <- list()
for (cl_idx in target_clusters) {
  core_genes <- names(which(mem_matrix[, cl_idx] >= 0.5))
  if (length(core_genes) > 0) {
    df_cl <- as.data.frame(expr_matrix[core_genes, , drop = FALSE])
    df_cl$Gene       <- rownames(df_cl)
    df_cl$Membership <- mem_matrix[core_genes, cl_idx]
    df_cl$Cluster    <- cluster_names[as.character(cl_idx)]
    traj_list[[as.character(cl_idx)]] <- df_cl
  }
}

traj_df <- bind_rows(traj_list) %>%
  pivot_longer(cols = c("T_0h", "T_24h", "T_336h"), names_to = "Time", values_to = "Expression") %>%
  mutate(
    Time = factor(Time, levels = c("T_0h", "T_24h", "T_336h"), labels = c("0h", "24h", "2w")),
    Cluster = factor(Cluster, levels = c("Cluster 2", "Cluster 4", "Cluster 5"))
  )

centroid_df <- traj_df %>%
  group_by(Cluster, Time) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

panel_a_fullwidth <- ggplot() +
  geom_line(
    data = traj_df, 
    aes(x = Time, y = Expression, group = Gene, color = Membership), 
    alpha = 0.35, linewidth = 0.5
  ) +
  geom_line(
    data = centroid_df, 
    aes(x = Time, y = Expression, group = 1), 
    color = "black", linewidth = 1.15
  ) +
  scale_color_gradientn(
    colors = c("#FF0000", "#FF7F00", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF", "#8B00FF"),
    name = "Membership"
  ) +
  facet_wrap(~ Cluster, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "Std. Expr.") +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8.5),
    axis.text.x = element_text(size = 9, face = "bold", margin = margin(t = 2)),
    axis.title.y = element_text(size = 9.5, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

target_order <- c("Cluster 2", "Cluster 4", "Cluster 5")

ora_top5_generatio <- combined_mfuzz_ora_all %>%
  mutate(
    Cluster_Short = str_replace(Cluster, " \\(.*\\)", ""),
    Description_wrapped = str_wrap(Description, width = 45)
  ) %>%
  filter(Cluster_Short %in% target_order) %>%
  group_by(Cluster_Short) %>%
  arrange(desc(GeneRatioNumeric), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    Cluster_Short = factor(Cluster_Short, levels = target_order)
  )

panel_b_unified <- ggplot(ora_top5_generatio, aes(x = GeneRatioNumeric, y = reorder(Description_wrapped, GeneRatioNumeric))) +
  geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.9) +
  scale_color_viridis_c(name = expression(-log[10](p[adj])), option = "magma", end = 0.85) +
  scale_size_continuous(name = "Gene Count", range = range(2.5, 5.5)) +
  scale_x_continuous(n.breaks = 4) +
  facet_wrap(~ Cluster_Short, scales = "free", ncol = 1) +
  theme_bw(base_size = 10) +
  labs(x = "Gene Ratio", y = NULL) +
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.9),
    axis.text.x = element_text(size = 8.5, margin = margin(t = 2)),
    axis.title.x = element_text(size = 9.5, face = "bold", margin = margin(t = 4)),
    strip.text = element_text(face = "bold", size = 10.5),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

shared_legend <- cowplot::get_legend(panel_b_unified)
panel_b_clean <- panel_b_unified + theme(legend.position = "none")

bottom_row <- cowplot::plot_grid(
  panel_b_clean, 
  shared_legend, 
  nrow = 1, 
  rel_widths = c(1, 0.22)
)

final_master_figure <- cowplot::plot_grid(
  panel_a_fullwidth, 
  bottom_row, 
  ncol = 1, 
  rel_heights = c(0.30, 1),
  labels = c("A", "B"),
  label_size = 14
)
#### Proteomic Mfuzz Clustering and ORA Analysis ####

set.seed(123)

expr_prot <- final_imputed_matrix
sig_prot  <- prot_master_df

lrt_mask_prot  <- !is.na(sig_prot$padj_LRT) & sig_prot$padj_LRT < 0.05
lrt_genes_prot <- unique(sig_prot$name[lrt_mask_prot])
lrt_genes_prot <- lrt_genes_prot[!is.na(lrt_genes_prot) & lrt_genes_prot != ""]

valid_genes_prot <- intersect(rownames(expr_prot), lrt_genes_prot)
expr_sig_prot    <- expr_prot[valid_genes_prot, ]
expr_sig_prot <- expr_sig_prot[complete.cases(expr_sig_prot[, 1:12]), ]

time_matrix_prot <- matrix(
  c(
    rowMeans(expr_sig_prot[, 1:4], na.rm = TRUE),
    rowMeans(expr_sig_prot[, 5:8], na.rm = TRUE),
    rowMeans(expr_sig_prot[, 9:12], na.rm = TRUE)
  ),
  ncol = 3,
  dimnames = list(rownames(expr_sig_prot), c("T_0h", "T_24h", "T_336h"))
)

eset_prot <- new("ExpressionSet", exprs = as.matrix(time_matrix_prot))
eset_std_prot <- standardise(eset_prot)

cl_prot_m2 <- mfuzz(eset_std_prot, c = 6, m = 2.00)

bg_symbols_prot <- rownames(expr_sig_prot)
bg_mapped_prot  <- bitr(bg_symbols_prot, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez_prot <- unique(bg_mapped_prot$ENTREZID)

run_cluster_ora_prot <- function(mfuzz_obj, cluster_num, membership_threshold = 0.5) {
  cluster_mem <- mfuzz_obj$membership[, cluster_num]
  target_symbols <- names(cluster_mem)[cluster_mem >= membership_threshold]
  
  mapped <- bitr(
    target_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  gene_entrez <- unique(mapped$ENTREZID)
  
  ora_res <- enrichPathway(
    gene          = gene_entrez,
    universe      = universe_entrez_prot,
    organism      = "human",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  
  return(ora_res)
}

ora_list_prot <- list()
cluster_labels_prot <- c(
  "Cluster 1",
  "Cluster 2",
  "Cluster 3",
  "Cluster 4",
  "Cluster 5",
  "Cluster 6"
)

for (i in 1:6) {
  res <- run_cluster_ora_prot(cl_prot_m2, cluster_num = i, membership_threshold = 0.5)
  df <- as.data.frame(res)
  if (nrow(df) > 0) {
    df_top <- df %>%
      slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>%
      mutate(
        Cluster = cluster_labels_prot[i],
        GeneRatioNumeric = sapply(GeneRatio, function(x) {
          parts <- as.numeric(strsplit(x, "/")[[1]])
          parts[1] / parts[2]
        }),
        Description_wrapped = str_wrap(Description, width = 38)
      )
    ora_list_prot[[i]] <- df_top
  }
}

combined_prot_ora_all <- bind_rows(ora_list_prot)

expr_matrix_prot <- exprs(eset_std_prot)
mem_matrix_prot  <- cl_prot_m2$membership

target_clusters_prot <- c(3, 4, 5)
cluster_names_prot   <- c("3" = "Cluster 3", "4" = "Cluster 4", "5" = "Cluster 5")

traj_list_prot <- list()
for (cl_idx in target_clusters_prot) {
  core_genes <- names(which(mem_matrix_prot[, cl_idx] >= 0.5))
  if (length(core_genes) > 0) {
    df_cl <- as.data.frame(expr_matrix_prot[core_genes, , drop = FALSE])
    df_cl$Gene       <- rownames(df_cl)
    df_cl$Membership <- mem_matrix_prot[core_genes, cl_idx]
    df_cl$Cluster    <- cluster_names_prot[as.character(cl_idx)]
    traj_list_prot[[as.character(cl_idx)]] <- df_cl
  }
}

traj_df_prot <- bind_rows(traj_list_prot) %>%
  pivot_longer(cols = c("T_0h", "T_24h", "T_336h"), names_to = "Time", values_to = "Expression") %>%
  mutate(
    Time = factor(Time, levels = c("T_0h", "T_24h", "T_336h"), labels = c("0h", "24h", "2w")),
    Cluster = factor(Cluster, levels = c("Cluster 3", "Cluster 4", "Cluster 5"))
  )

centroid_df_prot <- traj_df_prot %>%
  group_by(Cluster, Time) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

panel_a_fullwidth_prot <- ggplot() +
  geom_line(
    data = traj_df_prot, 
    aes(x = Time, y = Expression, group = Gene, color = Membership), 
    alpha = 0.35, linewidth = 0.5
  ) +
  geom_line(
    data = centroid_df_prot, 
    aes(x = Time, y = Expression, group = 1), 
    color = "black", linewidth = 1.15
  ) +
  scale_color_gradientn(
    colors = c("#FF0000", "#FF7F00", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF", "#8B00FF"),
    name = "Membership"
  ) +
  facet_wrap(~ Cluster, ncol = 3, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(x = NULL, y = "Std. Expr.") +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    axis.text.y = element_text(size = 8.5),
    axis.text.x = element_text(size = 9, face = "bold", margin = margin(t = 2)),
    axis.title.y = element_text(size = 9.5, face = "bold"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

target_order_prot <- c("Cluster 3", "Cluster 4", "Cluster 5")

ora_top5_prot <- combined_prot_ora_all %>%
  mutate(
    Cluster_Short = str_replace(Cluster, " \\(.*\\)", ""),
    Description_wrapped = str_wrap(Description, width = 45)
  ) %>%
  filter(Cluster_Short %in% target_order_prot) %>%
  group_by(Cluster_Short) %>%
  arrange(desc(GeneRatioNumeric), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    Cluster_Short = factor(Cluster_Short, levels = target_order_prot)
  )

panel_b_unified_prot <- ggplot(ora_top5_prot, aes(x = GeneRatioNumeric, y = reorder(Description_wrapped, GeneRatioNumeric))) +
  geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.9) +
  scale_color_viridis_c(name = expression(-log[10](p[adj])), option = "magma", end = 0.85) +
  scale_size_continuous(name = "Protein Count", range = range(2.5, 5.5)) +
  scale_x_continuous(n.breaks = 4) +
  facet_wrap(~ Cluster_Short, scales = "free", ncol = 1) +
  theme_bw(base_size = 10) +
  labs(x = "Gene Ratio", y = NULL) +
  theme(
    axis.text.y = element_text(size = 9, lineheight = 0.9),
    axis.text.x = element_text(size = 8.5, margin = margin(t = 2)),
    axis.title.x = element_text(size = 9.5, face = "bold", margin = margin(t = 4)),
    strip.text = element_text(face = "bold", size = 10.5),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
  )

shared_legend_prot <- cowplot::get_legend(panel_b_unified_prot)
panel_b_clean_prot <- panel_b_unified_prot + theme(legend.position = "none")

bottom_row_prot <- cowplot::plot_grid(
  panel_b_clean_prot, 
  shared_legend_prot, 
  nrow = 1, 
  rel_widths = c(1, 0.22)
)

final_master_figure_prot <- cowplot::plot_grid(
  panel_a_fullwidth_prot, 
  bottom_row_prot, 
  ncol = 1, 
  rel_heights = c(0.30, 1),
  labels = c("A", "B"),
  label_size = 14
)
#### Omnipath Protein-Protein Interaction Network (PKN) Filtering ####

omnipath_net_ke <- omnipath_interactions(
  datasets = c('omnipath', 'kinaseextra', 'signor'),
  directed = TRUE,
  signed = TRUE,
  genesymbols = TRUE)

omnipath_filtered <- omnipath_net_ke %>%
  filter(curation_effort >= 2) %>%
  filter(consensus_direction == TRUE) %>%
  filter(consensus_stimulation | consensus_inhibition) %>%
  filter(source_genesymbol != target_genesymbol) %>%
  mutate(interaction = case_when(
    consensus_stimulation & !consensus_inhibition ~ 1,
    consensus_inhibition & !consensus_stimulation ~ -1,
    TRUE ~ NA_integer_
  )) %>%
  filter(!is.na(interaction)) %>%
  select(source = source_genesymbol, 
         interaction, 
         target = target_genesymbol) %>%
  distinct()

omnipath_filtered_nodes <- unique(c(omnipath_filtered$source, omnipath_filtered$target))

ptk2b_edges <- omnipath_filtered %>%
  filter(source == "PTK2B" | target == "PTK2B")
#### CollecTRI and PROGENy Activity Inference ####

df_2w_0h_collectri <- as.data.frame(res_2w_0h) %>%
  rownames_to_column(var = "ENSEMBL")

df_2w_0h_collectri$Symbol <- mapIds(
  org.Hs.eg.db,
  keys = df_2w_0h_collectri$ENSEMBL,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

clean_stat_df <- df_2w_0h_collectri %>%
  filter(!is.na(Symbol) & !is.na(stat)) %>%
  group_by(Symbol) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  column_to_rownames(var = "Symbol")

stat_matrix <- clean_stat_df[, 'stat', drop = FALSE]

omnipath_filtered_nodes <- unique(c(omnipath_filtered$source, omnipath_filtered$target))

collectri_net <- get_collectri(organism = 'human', split_complexes = FALSE)

tf_activities <- run_ulm(
  mat     = stat_matrix,
  network = collectri_net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH"))

top_10_collectri <- tf_activities %>%
  filter(p_value_adj < 0.05 & score > 0) %>%
  arrange(desc(score)) %>%
  slice_head(n = 10)

top_50_collectri <- tf_activities %>%
  filter(p_value_adj < 0.05) %>%
  filter(source %in% omnipath_filtered_nodes) %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 50)

plot_data_collectri <- tf_activities %>%
  filter(p_value_adj < 0.05 & score > 0) %>%
  filter(source %in% omnipath_filtered_nodes) %>%
  mutate(
    neg_log_fdr = -log10(p_value_adj),
    source = reorder(source, score)
  ) %>%
  arrange(desc(score)) %>%
  slice_head(n = 10)

collectri_prism_plot <- ggplot(plot_data_collectri, aes(x = score, y = source)) +
  geom_segment(aes(x = 0, xend = score, y = source, yend = source), 
               color = "grey60", linewidth = 0.8) +
  geom_point(aes(size = neg_log_fdr, color = score), alpha = 1) +
  scale_color_viridis_c(option = "plasma", direction = 1, end = 0.9) +
  scale_size_continuous(range = c(4, 9)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 12)) +
  labs(
    x = "CollecTRI Activity Score (ULM t-value)",
    y = "",
    color = "Activity Score",
    size = "Significance\n(-log10 FDR)"
  ) +
  theme_prism(base_size = 11, base_family = "Arial", border = FALSE) +
  theme(
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right"
  )

carnival_tf_input <- top_50_collectri$score
names(carnival_tf_input) <- top_50_collectri$source

top_50_collectri_nomyc <- tf_activities %>%
  filter(p_value_adj < 0.05) %>%
  filter(source %in% omnipath_filtered_nodes) %>%
  filter(source != "MYC") %>%
  arrange(desc(abs(score))) %>%
  slice_head(n = 50)

carnival_tf_input_nomyc <- top_50_collectri_nomyc$score
names(carnival_tf_input_nomyc) <- top_50_collectri_nomyc$source

progeny_net <- get_progeny(organism = 'human', top = 100)

pathway_activities <- run_ulm(
  mat     = stat_matrix,
  network = progeny_net,
  .source = "source",
  .target = "target",
  .mor    = "weight",
  minsize = 5
)

pathway_map <- list(
  EGFR       = c("EGFR", "ERBB2"),
  Hypoxia    = "HIF1A",
  "JAK-STAT" = c("JAK1", "JAK2", "JAK3"),
  MAPK       = c("BRAF", "ARAF", "RAF1"),
  NFkB       = "NFKB1",
  PI3K       = c("PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"),
  TGFb       = c("TGFBR1", "TGFBR2", "BMPR1A", "BMPR1B", "BMPR2"),
  TNFa       = c("TNFRSF1A", "TNFRSF1B"),
  Trail      = c("CASP8", "CASP10"),
  VEGF       = c("FLT1", "FLT3", "KDR", "PDGFRA", "PDGFRB"),
  p53        = "TP53",
  Androgen   = "AR",
  Estrogen   = c("ESR1", "ESR2"),
  WNT        = "DVL1"
)

expanded_weights <- pathway_activities %>%
  filter(source %in% names(pathway_map)) %>%
  mutate(nodes = map(source, ~ pathway_map[[.x]])) %>%
  unnest(nodes) %>%
  filter(nodes %in% omnipath_filtered_nodes)

progeny_vector <- expanded_weights$score
names(progeny_vector) <- expanded_weights$nodes
progeny_vector <- progeny_vector[!duplicated(names(progeny_vector))]

key_nodes <- c("PTEN", "MYC", "RELA", "NFKB1", "NFKB2", "REL", "RELB", "KAT7", "JUN", "IRF3", "STAT3")

interaction_summary <- omnipath_filtered %>%
  filter(source %in% key_nodes | target %in% key_nodes) %>%
  group_by(Node = ifelse(source %in% key_nodes & target %in% key_nodes, "Internal", 
                         ifelse(source %in% key_nodes, source, target))) %>%
  summarise(
    Total_Edges = n(),
    As_Source = sum(source %in% key_nodes),
    As_Target = sum(target %in% key_nodes)
  )

#### RcisTarget Motif Enrichment (2w vs 0h) ####

up_genes_0h <- master_df %>% 
  filter(LFC_2w_vs_0h > 1 & Padj_2w_vs_0h < 0.05) %>% 
  filter(!is.na(Symbol)) %>% 
  pull(Symbol) %>%
  unique()

down_genes_0h <- master_df %>% 
  filter(LFC_2w_vs_0h < -1 & Padj_2w_vs_0h < 0.05) %>% 
  filter(!is.na(Symbol)) %>% 
  pull(Symbol) %>%
  unique()

geneLists_0h <- list(
  Chronic_0h_Up = up_genes_0h,
  Chronic_0h_Down = down_genes_0h
)

data(motifAnnotations_hgnc_v9, package="RcisTarget")
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9

motifRankings <- importRankings("~/Downloads/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

motif_enrich_0h <- cisTarget(
  geneLists_0h, 
  motifRankings,
  motifAnnot = motifAnnotations_hgnc
)

results_0h_df <- as.data.frame(motif_enrich_0h)

artifact_regex <- "POLR2|TAF[0-9]|SIN3|TBP|CD59|HNRNP|EXOSC|ZNF239|ZNF583|AVEN"

plot_data_rcis_dedup <- results_0h_df %>%
  filter(geneSet == "Chronic_0h_Up") %>%
  mutate(
    Clean_TF = str_remove_all(TF_highConf, "\\s*\\(directAnnotation\\)\\s*\\.?"),
    Clean_TF = str_remove_all(Clean_TF, "\\s*\\(inferredBy_Orthology\\)\\s*\\.?"),
    Clean_TF = str_trim(Clean_TF),
    Primary_TF = str_extract(Clean_TF, "^[A-Z0-9]+")
  ) %>%
  filter(!is.na(Primary_TF) & Primary_TF != "" & Primary_TF != "-") %>%
  filter(!str_detect(Clean_TF, artifact_regex)) %>%
  group_by(Primary_TF) %>%
  slice_max(order_by = NES, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10) %>%
  mutate(
    Display_Label = str_replace_all(Clean_TF, "; ", " / "),
    Display_Label = str_wrap(Display_Label, width = 28),
    Display_Label = reorder(Display_Label, NES)
  )

rcistarget_prism_plot <- ggplot(plot_data_rcis_dedup, aes(x = NES, y = Display_Label)) +
  geom_segment(
    aes(x = 0, xend = NES, y = Display_Label, yend = Display_Label), 
    color = "grey60", 
    linewidth = 0.8
  ) +
  geom_point(
    aes(size = nEnrGenes, color = NES), 
    alpha = 1.0
  ) +
  scale_color_viridis_c(option = "viridis", direction = 1, end = 0.9) +
  scale_size_continuous(range = c(4, 9)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 9.5)) +
  labs(
    x = "Maximum Normalized Enrichment Score (NES)",
    y = "",
    color = "NES Score",
    size = "Target Genes\n(nEnrGenes)"
  ) +
  theme_prism(base_size = 11, base_family = "Arial", border = FALSE) +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right",
    axis.text.y = element_text(face = "bold", size = 9, lineheight = 0.85)
  )

top50_table_0h <- results_0h_df %>%
  filter(geneSet == "Chronic_0h_Up") %>%
  mutate(
    Clean_TF = str_remove_all(TF_highConf, "\\s*\\(directAnnotation\\)\\s*\\.?"),
    Clean_TF = str_remove_all(Clean_TF, "\\s*\\(inferredBy_Orthology\\)\\s*\\.?"),
    Clean_TF = str_trim(Clean_TF)
  ) %>%
  filter(Clean_TF != "" & !is.na(Clean_TF) & Clean_TF != "-") %>%
  filter(!str_detect(Clean_TF, artifact_regex)) %>%
  group_by(Clean_TF) %>%
  slice_max(order_by = NES, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(NES)) %>%
  mutate(
    Clean_TF = str_replace_all(Clean_TF, "; ", "/"),
    Rank = row_number()
  ) %>%
  select(Rank, Motif_TF = Clean_TF, NES, nEnrGenes, AUC) %>%
  head(50)

# Composite Panel Assembly (CollecTRI + RcisTarget)
panel_a_collectri <- collectri_prism_plot + 
  theme(plot.title = element_blank()) + 
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))

panel_b_rcistarget <- rcistarget_prism_plot + 
  theme(plot.title = element_blank()) + 
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))

composite_collectri_rcis_figure <- panel_a_collectri + panel_b_rcistarget + 
  plot_layout(widths = c(1, 1.2)) +
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag = element_text(face = "bold", size = 16)
  )
#### CARNIVAL Causal Network Reconstruction ####

Sys.setenv(GUROBI_HOME = "/iridisfs/ixsoftware/gurobi/11.0.3/install/gurobi1103/linux64")

if (!file.exists("/iridisfs/ixsoftware/gurobi/11.0.3/install/gurobi1103/linux64/lib/libgurobi110.so")) {
  stop("[FATAL] Gurobi shared library not found at specified path.")
}

base_dir <- "/iridisfs/scratch/bab1u25/Carnival_2w_0h"

pkn             <- omnipath_filtered
T2w_T0h_tfs     <- carnival_tf_input_nomyc
T2w_T0h_progeny <- progeny_vector

my_opts <- defaultLpSolveCarnivalOptions()

my_opts$solver        <- "gurobi"
my_opts$solverPath    <- "/iridisfs/ixsoftware/gurobi/11.0.3/install/gurobi1103/linux64/bin/gurobi_cl"
my_opts$workdir       <- file.path(base_dir, "Carnival_tmp")
my_opts$outputFolder  <- file.path(base_dir, "Carnival_output")
my_opts$betaWeight    <- 0.3
my_opts$timelimit     <- 7200
my_opts$threads       <- 22
my_opts$mipGap        <- 0.05
my_opts$poolrelGap    <- 0.0001
my_opts$limitPop      <- 500
my_opts$poolCap       <- 100
my_opts$cleanTmpFiles <- 1
my_opts$keepLPFiles   <- 0

pten_pos <- c(PTEN = 1)

carnival_baseline_res <- runVanillaCarnival(
  measurements          = T2w_T0h_tfs,
  perturbations         = pten_pos,
  priorKnowledgeNetwork = pkn,
  weights               = T2w_T0h_progeny,
  carnivalOptions       = my_opts
)
#### CARNIVAL Kinase Cascade Extraction and Network Filtering ####

res <- carnival_baseline_res
sif  <- res$weightedSIF      
attr <- res$nodesAttributes

attr_dict <- setNames(attr$AvgAct, attr$Node)

kinases <- OmnipathR::import_omnipath_annotations(resources = "kinase.com")
kinase_list <- unique(kinases$genesymbol)

kin_nodes <- attr[attr$Node %in% kinase_list, ]
kin_nodes$AbsAvgAct <- abs(kin_nodes$AvgAct)

deg <- table(c(sif$Node1, sif$Node2))
kin_nodes$Degree <- as.integer(deg[kin_nodes$Node])
kin_nodes$Degree[is.na(kin_nodes$Degree)] <- 0L

kin_nodes <- kin_nodes[order(-kin_nodes$AbsAvgAct, -kin_nodes$Degree), ]

top_kinases <- kin_nodes[kin_nodes$AbsAvgAct == 100, ]

hc <- sif[sif$Weight == 100, ]

g <- graph_from_data_frame(hc[, c("Node1", "Node2", "Sign")], directed = TRUE)

input_node <- attr$Node[attr$NodeType == "P"]     
target_hub <- "MAPK14"

sp <- shortest_paths(g, from = input_node, to = target_hub, mode = "out")
core_path <- names(sp$vpath[[1]])

core <- unique(c(core_path, "CSNK2A1", "RELA", "SNAI1"))

one_hop <- unique(unlist(lapply(core, function(n) {
  c(names(neighbors(g, n, mode = "in")),
    names(neighbors(g, n, mode = "out")))
})))
one_hop <- setdiff(one_hop, core)

additions <- c("IKBKB", "STAT3", "TGFBR2", "TP53", "NFKB1", "RAF1", "SMAD4",
               "E2F1", "TNF", "SP1", "PRKD1", "PLK1", "IRF5", "HDAC3", "STK11")

selected <- union(core, additions)
sub_g <- induced_subgraph(g, selected)

edges_df <- as_data_frame(sub_g, what = "edges")
edges_df$Interaction <- ifelse(edges_df$Sign == 1, "activates", "inhibits")
names(edges_df)[names(edges_df) == "from"] <- "Source"
names(edges_df)[names(edges_df) == "to"]   <- "Target"
edges_df <- edges_df[, c("Source", "Interaction", "Target", "Sign")]

node_type_lookup <- setNames(attr$NodeType, attr$Node)

node_df <- data.frame(
  Node     = selected,
  AvgAct   = attr_dict[selected],
  Role     = ifelse(selected %in% core, "Core cascade", "Extension"),
  IsKinase = ifelse(selected %in% kinase_list, "Y", "N"),
  NodeType = ifelse(node_type_lookup[selected] == "P", "Input",
                    ifelse(node_type_lookup[selected] == "M", "Terminal (sink)",
                           "Intermediate"))
)
#### PTK2B Evidence and Downstream Target Panels ####

rna_master <- master_df
target_gene <- "PTK2B"

generate_masked_kinome_volcano <- function(data_frame, contrast_prefix, kinase_list, target_gene = "PTK2B", panel_title, x_lim = 3, y_lim = 25, lfc_thresh = 0.5) {
  
  x_col <- paste0("LFC_", contrast_prefix)
  p_col <- paste0("Padj_", contrast_prefix)
  
  df_clean <- data_frame[!is.na(data_frame[[p_col]]), ]
  
  df_clean$Gene_Class <- "Background"
  df_clean$Gene_Class[df_clean$Symbol %in% kinase_list] <- "Kinome"
  df_clean$Gene_Class[df_clean$Symbol == target_gene] <- "Target"
  
  df_clean$Gene_Class <- factor(df_clean$Gene_Class, levels = c("Background", "Kinome", "Target"))
  df_clean <- df_clean[order(df_clean$Gene_Class), ]
  
  target_df <- df_clean[df_clean$Gene_Class == "Target", ]
  
  p <- ggplot(df_clean, aes(x = .data[[x_col]], y = -log10(.data[[p_col]]))) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "#A9A9A9", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#A9A9A9", linewidth = 0.5) +
    geom_point(aes(color = Gene_Class, size = Gene_Class, alpha = Gene_Class)) +
    scale_color_manual(values = c("Background" = "#D3D3D3", 
                                  "Kinome" = "#1F77B4", 
                                  "Target" = "#FFD700")) +
    scale_size_manual(values = c("Background" = 0.8, 
                                 "Kinome" = 1.8,  
                                 "Target" = 4.0)) + 
    scale_alpha_manual(values = c("Background" = 0.15, 
                                  "Kinome" = 0.85, 
                                  "Target" = 1.0)) +
    geom_text_repel(data = target_df, 
                    aes(label = Symbol), 
                    size = 5, 
                    fontface = "bold.italic",
                    color = "black",
                    box.padding = 1.5,
                    point.padding = 0.5,
                    segment.color = "black",
                    segment.size = 0.6,
                    min.segment.length = 0) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_lim), expand = FALSE) +
    labs(
      x = expression(bold(Log[2]~Fold~Change)), 
      y = expression(bold(-Log[10](Adjusted~P-value)))) +
    theme_prism(base_size = 14, border = TRUE) +  
    theme(
      legend.position = "none",
      axis.line = element_line(linewidth = 0.5), 
      axis.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold")
    )
  
  return(p)
}

plot_rna_10a <- generate_masked_kinome_volcano(
  data_frame = rna_master,
  contrast_prefix = "2w_vs_0h",   
  kinase_list = kinase_list, 
  target_gene = "PTK2B",
  panel_title = "Transcriptomic Kinase Landscape (0h vs 2w)", 
  x_lim = 3,   
  y_lim = 25,  
  lfc_thresh = 0.5 
)

top_kinases_rna <- rna_master %>%
  filter(Symbol %in% kinase_list) %>%
  arrange(desc(LFC_2w_vs_0h)) %>%
  select(Symbol, LFC_2w_vs_0h, Padj_2w_vs_0h) %>%
  head(10)

vsd_expr <- vst_mat_unique

ptk2b_rna <- data.frame(
  Sample = colnames(vsd_expr),
  Expression = as.numeric(vsd_expr["PTK2B", ])
)

ptk2b_df <- ptk2b_rna %>%
  mutate(
    Timepoint = case_when(
      Sample %in% c("PTEN-2", "PTEN-3", "PTEN-4") ~ "0h",
      Sample %in% c("PTEN-6", "PTEN-7", "PTEN-8") ~ "24h",
      Sample %in% c("PTEN-9", "PTEN-10", "PTEN-11", "PTEN-12") ~ "2w",
      TRUE ~ "Unknown" 
    ),
    Timepoint = factor(Timepoint, levels = c("0h", "24h", "2w"))
  )

max_val <- max(ptk2b_df$Expression, na.rm = TRUE)
step <- 0.08 

ptk2b_plot <- ggplot(ptk2b_df, aes(x = Timepoint, y = Expression, fill = Timepoint)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "black", linewidth = 0.8, width = 0.5) +
  geom_jitter(width = 0.1, size = 4, shape = 21, color = "black", stroke = 1) +
  geom_signif(comparisons = list(c("0h", "24h")), 
              annotations = "ns", 
              y_position = max_val + (step * 1), 
              tip_length = 0.02, textsize = 4, vjust = 0) +
  geom_signif(comparisons = list(c("24h", "2w")), 
              annotations = "*", 
              y_position = max_val + (step * 2), 
              tip_length = 0.02, textsize = 5, vjust = 0.5) +
  geom_signif(comparisons = list(c("0h", "2w")), 
              annotations = "**", 
              y_position = max_val + (step * 3.5), 
              tip_length = 0.02, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = c("0h" = "#FDE725FF", 
                               "24h" = "#21908CFF", 
                               "2w" = "#440154FF")) +
  theme_prism(base_size = 14, border = FALSE) + 
  labs(
    x = NULL, 
    y = "PTK2B mRNA Expression (VST)" 
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.02)), 
    breaks = seq(7.6, 8.6, by = 0.2) 
  ) +
  coord_cartesian(
    ylim = c(min(ptk2b_df$Expression) - 0.1, max_val + (step * 3.5) + 0.05) 
  ) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5), 
    axis.title.y = element_text(face = "bold", margin = margin(r = 15)),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12) 
  )

pyk2_activators <- c(
  "ITGB1", "ITGB3", "ITGAV", "PXN", "SRC", "FYN",
  "CALM1", "CALM2", "CALM3", "PRKCA", "PRKCB", "PRKCG", "CCR2", "CXCR1", "CXCR2",
  "IL1R1", "TNFRSF1A", "TGFBR1", "IRAK2", "IRAK3"
)

valid_genes <- pyk2_activators[pyk2_activators %in% rownames(vsd_expr)]
activator_matrix <- vsd_expr[valid_genes, ]

mean_0h  <- rowMeans(activator_matrix[, 1:3])
mean_24h <- rowMeans(activator_matrix[, 4:6])
mean_2w  <- rowMeans(activator_matrix[, 7:10])

activator_summary <- data.frame(
  Gene = rownames(activator_matrix),
  Mean_0h = mean_0h,
  Mean_24h = mean_24h,
  Mean_2w = mean_2w
)

integrin_genes <- grep("^ITG", rownames(vsd_expr), value = TRUE)
integrin_matrix <- vsd_expr[integrin_genes, ]

mean_0h_int  <- rowMeans(integrin_matrix[, 1:3])
mean_24h_int <- rowMeans(integrin_matrix[, 4:6])
mean_2w_int  <- rowMeans(integrin_matrix[, 7:10])

integrin_summary <- data.frame(
  Gene = rownames(integrin_matrix),
  Mean_0h = mean_0h_int,
  Mean_24h = mean_24h_int,
  Mean_2w = mean_2w_int
)

integrin_summary$Shift_0h_to_2w <- integrin_summary$Mean_2w - integrin_summary$Mean_0h
integrin_summary <- integrin_summary[order(-integrin_summary$Shift_0h_to_2w), ]

integrin_stats <- rna_master %>%
  filter(grepl("^ITG", Symbol)) %>%
  select(Symbol, LFC_2w_vs_0h, Padj_2w_vs_0h) %>%
  arrange(desc(LFC_2w_vs_0h))

pyk2_activators_triads <- c(
  "PXN", "SRC", "FYN",                      
  "CALM1", "CALM2", "CALM3",                
  "PRKCA", "PRKCB", "PRKCG",                
  "CCR2", "CXCR1", "CXCR2",                 
  "IL1R1", "TNFRSF1A", "TGFBR1",            
  "IRAK2", "IRAK3"                          
)

target_stats <- rna_master %>%
  filter(Symbol %in% pyk2_activators_triads) %>%
  select(Symbol, LFC_2w_vs_0h, Padj_2w_vs_0h) %>%
  arrange(desc(LFC_2w_vs_0h))

akt_non_canon <- c("PRKDC", "ATM", "TBK1", "IKBKE", "ILK", "PTK6", "TNK2", "SRC", "CSNK2A1", "CSNK2A2", "PAK1", "PAK2", "CDK2")

akt_stats <- rna_master %>%
  filter(Symbol %in% akt_non_canon) %>%
  dplyr::select(Symbol, LFC_2w_vs_0h, Padj_2w_vs_0h) %>%
  arrange(desc(LFC_2w_vs_0h))

prot <- final_imputed_matrix

panel_10c_genes <- c(
  "ITGAM", "ITGB2",           
  "ITGA5", "ITGAV", "ITGB1",  
  "PXN", "SRC",               
  "GRB2", "SOS1",             
  "IRAK2", "IRAK3", "TGFBR1", 
  "CALM2"                     
)

mrna_sub_10c <- vsd_expr[panel_10c_genes, ]

prot_sub_10c <- matrix(NA, nrow = length(panel_10c_genes), ncol = ncol(prot), 
                       dimnames = list(panel_10c_genes, colnames(prot)))

for (gene in panel_10c_genes) {
  if (gene %in% rownames(prot)) {
    prot_sub_10c[gene, ] <- as.numeric(prot[gene, ])
  }
}

safe_scale <- function(x) {
  if(all(is.na(x))) return(x)
  return((x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
}

mrna_z_10c <- t(apply(mrna_sub_10c, 1, safe_scale))
prot_z_10c <- t(apply(prot_sub_10c, 1, safe_scale))

timepoints_mrna <- factor(c(rep("0h", 3), rep("24h", 3), rep("2w", 4)), levels = c("0h", "24h", "2w"))
timepoints_prot <- factor(c(rep("0h", 4), rep("24h", 4), rep("2w", 4)), levels = c("0h", "24h", "2w"))

shared_colors <- list(Timepoint = c("0h" = "#FDE725FF", "24h" = "#21908CFF", "2w" = "#440154FF"))

col_anno_mrna_10c = HeatmapAnnotation(
  Timepoint = timepoints_mrna,
  col = shared_colors,
  show_annotation_name = FALSE
)

col_anno_prot_10c = HeatmapAnnotation(
  Timepoint = timepoints_prot,
  col = shared_colors,
  show_annotation_name = FALSE
)

col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))

ht_mrna_10c = Heatmap(
  mrna_z_10c, 
  name = "Z-Score", 
  column_title = "Transcriptomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  top_annotation = col_anno_mrna_10c,
  cluster_columns = FALSE, 
  cluster_rows = FALSE,      
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11, fontface = "bold.italic"), 
  width = unit(6, "cm"), 
  height = unit(7, "cm"), 
  show_column_names = FALSE,
  show_row_dend = FALSE
)

ht_prot_10c = Heatmap(
  prot_z_10c, 
  name = "Z-Score",
  column_title = "Proteomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  na_col = "grey80", 
  top_annotation = col_anno_prot_10c,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  width = unit(6, "cm"), 
  height = unit(7, "cm"), 
  show_row_names = FALSE,                
  show_heatmap_legend = FALSE,           
  show_column_names = FALSE,  
  heatmap_legend_param = list(direction = "vertical")
)

final_layout_10c = ht_mrna_10c + ht_prot_10c

panel_10d_genes <- c(
  "IKBKB", "IKBKE", "TBK1",   
  "RELA", "NFKB1", "STAT3",   
  "TNF", "IL6",       
  "CCL2", "TGFB1"             
)

mrna_sub_10d <- vsd_expr[panel_10d_genes, ]

prot_sub_10d <- matrix(NA, nrow = length(panel_10d_genes), ncol = ncol(prot), 
                       dimnames = list(panel_10d_genes, colnames(prot)))

for (gene in panel_10d_genes) {
  if (gene %in% rownames(prot)) {
    prot_sub_10d[gene, ] <- as.numeric(prot[gene, ])
  }
}

mrna_z_10d <- t(apply(mrna_sub_10d, 1, safe_scale))
prot_z_10d <- t(apply(prot_sub_10d, 1, safe_scale))

col_anno_mrna_10d = HeatmapAnnotation(Timepoint = timepoints_mrna, col = shared_colors, show_annotation_name = FALSE)
col_anno_prot_10d = HeatmapAnnotation(Timepoint = timepoints_prot, col = shared_colors, show_annotation_name = FALSE)

ht_mrna_10d = Heatmap(
  mrna_z_10d, 
  name = "Z-Score", 
  column_title = "Transcriptomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  top_annotation = col_anno_mrna_10d,
  cluster_columns = FALSE, 
  cluster_rows = FALSE,      
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11, fontface = "bold.italic"), 
  width = unit(6, "cm"), 
  height = unit(6, "cm"),  
  show_column_names = FALSE
)

ht_prot_10d = Heatmap(
  prot_z_10d, 
  name = "Z-Score",
  column_title = "Proteomics",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  na_col = "grey80",         
  top_annotation = col_anno_prot_10d,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  width = unit(6, "cm"), 
  height = unit(6, "cm"), 
  show_row_names = FALSE,                
  show_heatmap_legend = FALSE,           
  show_column_names = FALSE
)

final_layout_10d = ht_mrna_10d + ht_prot_10d
#### siPTK2B Knockdown and Downstream Perturbation Analysis ####

baseline <- read.csv("RNAseq_si/baseline(in).csv", check.names = FALSE) %>%
  slice(-1) %>%
  mutate(across(-Gene, as.numeric))

siptk2b <- read.csv("RNAseq_si/data(in).csv", check.names = FALSE) %>%
  slice(-1) %>%
  mutate(across(-Gene, as.numeric))

merged_data <- inner_join(baseline, siptk2b, by = "Gene", suffix = c("_Ctrl", "_siPTK2B")) %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene, .keep_all = TRUE)

merged_data <- as.data.frame(merged_data)
rownames(merged_data) <- merged_data$Gene
count_matrix <- merged_data[, -1]

resistant_counts <- count_matrix[, c("rep4_Ctrl", "rep5_Ctrl", "rep6_Ctrl", 
                                     "rep4_siPTK2B", "rep5_siPTK2B", "rep6_siPTK2B")]

resistant_metadata <- data.frame(
  row.names = colnames(resistant_counts),
  Condition = factor(c(rep("Control", 3), rep("siPTK2B", 3)), levels = c("Control", "siPTK2B"))
)

dds_resistant <- DESeqDataSetFromMatrix(countData = round(resistant_counts),
                                        colData = resistant_metadata,
                                        design = ~ Condition)

keep_res <- rowSums(counts(dds_resistant) >= 10) >= 3
dds_resistant <- dds_resistant[keep_res,]
dds_resistant <- DESeq(dds_resistant)
res_resistant <- results(dds_resistant, contrast = c("Condition", "siPTK2B", "Control"))

vsd_resistant <- vst(dds_resistant, blind = FALSE)
vst_mat_resistant <- assay(vsd_resistant)

plotPCA(vsd_resistant, intgroup="Condition")

ptk2b_kd_rna <- data.frame(
  Sample = colnames(vst_mat_resistant),
  Expression = as.numeric(vst_mat_resistant["PTK2B", ])
)

ptk2b_kd_df <- ptk2b_kd_rna %>%
  mutate(
    Condition = case_when(
      Sample %in% c("rep4_Ctrl", "rep5_Ctrl", "rep6_Ctrl") ~ "siCTRL",
      Sample %in% c("rep4_siPTK2B", "rep5_siPTK2B", "rep6_siPTK2B") ~ "siPTK2B",
      TRUE ~ "Drop" 
    )
  ) %>%
  filter(Condition != "Drop") %>%
  mutate(Condition = factor(Condition, levels = c("siCTRL", "siPTK2B")))

max_val_sirna <- max(ptk2b_kd_df$Expression, na.rm = TRUE)
step_sirna <- 0.08 

ptk2b_kd_plot <- ggplot(ptk2b_kd_df, aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "black", linewidth = 0.8, width = 0.5) +
  geom_jitter(width = 0.1, size = 4, shape = 21, color = "black", stroke = 1) +
  geom_signif(comparisons = list(c("siCTRL", "siPTK2B")), 
              annotations = "***", 
              y_position = max_val_sirna + (step_sirna * 1.5), 
              tip_length = 0.02, textsize = 6, vjust = 0.5) +
  scale_fill_manual(values = c("siCTRL" = "#440154FF", 
                               "siPTK2B" = "pink")) +
  theme_prism(base_size = 14, border = FALSE) +
  labs(
    x = NULL, 
    y = "PTK2B mRNA Expression (VST)" 
  ) +
  coord_cartesian(
    ylim = c(min(ptk2b_kd_df$Expression) - 0.1, max_val_sirna + (step_sirna * 3) + 0.05) 
  ) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5), 
    axis.title.y = element_text(face = "bold", margin = margin(r = 15)),
    axis.text.x = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 12)
  )

res_df_sirna <- as.data.frame(res_resistant)
res_df_sirna$symbol <- rownames(res_df_sirna)
res_df_sirna <- res_df_sirna[!is.na(res_df_sirna$padj), ]

res_df_sirna$Significance <- "Not Significant"
res_df_sirna$Significance[res_df_sirna$log2FoldChange > 1 & res_df_sirna$padj < 0.05] <- "Upregulated"
res_df_sirna$Significance[res_df_sirna$log2FoldChange < -1 & res_df_sirna$padj < 0.05] <- "Downregulated"
res_df_sirna$Significance <- factor(res_df_sirna$Significance, levels = c("Downregulated", "Not Significant", "Upregulated"))

top_genes_sirna <- res_df_sirna %>%
  filter(Significance != "Not Significant") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  slice_head(n = 10)

volcano_plot_sirna <- ggplot(res_df_sirna, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Downregulated" = "navy", "Not Significant" = "grey85", "Upregulated" = "firebrick3")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(data = top_genes_sirna, aes(label = symbol),
                  size = 3.5,
                  fontface = "italic",
                  box.padding = 0.8,
                  point.padding = 0.5,
                  force = 10,
                  max.overlaps = 15,
                  min.segment.length = 0,
                  segment.color = "grey40") +
  theme_prism(base_size = 14, border = TRUE) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10(Adjusted P-value)") +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12, face = "bold"))

upstream_genes_si <- c(
  "ITGAM", "ITGB2",           
  "ITGA5", "ITGAV", "ITGB1",  
  "PXN", "SRC",               
  "GRB2", "SOS1",             
  "IRAK2", "IRAK3", "TGFBR1", 
  "CALM2"                     
)

downstream_genes_si <- c(
  "IKBKB", "IKBKE", "TBK1",   
  "RELA", "NFKB1", "STAT3",   
  "TNF", "IL6",               
  "CCL2", "TGFB1"             
)

all_target_genes_si <- c(upstream_genes_si, downstream_genes_si)
target_samples_si <- c("rep4_Ctrl", "rep5_Ctrl", "rep6_Ctrl", 
                       "rep4_siPTK2B", "rep5_siPTK2B", "rep6_siPTK2B")

mrna_sub_si <- vst_mat_resistant[all_target_genes_si, target_samples_si]
mrna_z_si <- t(apply(mrna_sub_si, 1, safe_scale))

conditions_sirna <- factor(c(rep("siCTRL", 3), rep("siPTK2B", 3)), levels = c("siCTRL", "siPTK2B"))
condition_colors_si <- list(Condition = c("siCTRL" = "#440154FF", "siPTK2B" = "pink"))

col_anno_sirna_hm = HeatmapAnnotation(
  Condition = conditions_sirna,
  col = condition_colors_si,
  show_annotation_name = FALSE,
  simple_anno_size = unit(0.5, "cm")
)

row_split_factor_si <- factor(
  c(rep("Upstream Activators", length(upstream_genes_si)), 
    rep("Downstream Effectors", length(downstream_genes_si))),
  levels = c("Upstream Activators", "Downstream Effectors") 
)

ht_knockdown = Heatmap(
  mrna_z_si, 
  name = "Z-Score", 
  col = col_fun,
  top_annotation = col_anno_sirna_hm,
  row_split = row_split_factor_si,
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_gap = unit(6, "mm"), 
  row_title_rot = 0, 
  cluster_columns = FALSE, 
  cluster_rows = FALSE,      
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 11, fontface = "bold.italic"), 
  width = unit(6, "cm"), 
  height = unit(12, "cm"), 
  show_column_names = FALSE,
  heatmap_legend_param = list(direction = "vertical")
)

msigdb_h <- msigdbr(species = "Homo sapiens", collection  = "H")
hallmark <- split(msigdb_h$gene_symbol, msigdb_h$gs_name)

sirna_param <- gsvaParam(
  exprData = vst_mat_resistant,
  geneSets = hallmark,
  kcdf = "Gaussian",
  minSize = 10,
  maxSize = 500,
  tau = 1,
  absRanking = FALSE,
  maxDiff = TRUE
)

resistant_gsva_results_sirna <- gsva(sirna_param)

target_pathways <- c(
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE"
)

filtered_pathways_sirna <- intersect(target_pathways, rownames(resistant_gsva_results_sirna))
sirna_sub <- resistant_gsva_results_sirna[filtered_pathways_sirna, ]

timepoints_sirna_gsva <- factor(
  c(rep("Baseline", 3), rep("siRNA", 3)),
  levels = c("Baseline", "siRNA")
)

sirna_colors_gsva <- list(Timepoint = c("Baseline" = "#440154FF", "siRNA" = "pink"))

col_anno_sirna_gsva = HeatmapAnnotation(
  Timepoint = timepoints_sirna_gsva,
  col = sirna_colors_gsva,
  show_annotation_name = FALSE
)

ht_sirna_gsva = Heatmap(
  t(scale(t(sirna_sub))), 
  name = "Enrichment Score", 
  column_title = "Resistant (2W)",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  col = col_fun,
  top_annotation = col_anno_sirna_gsva,
  cluster_columns = FALSE, 
  show_row_names = TRUE,
  row_names_side = "left",
  row_labels = clean_pathway_names(rownames(sirna_sub)),            
  row_names_gp = gpar(fontsize = 11),      
  width = unit(6, "cm"),
  height = unit(4.5, "cm"),
  show_column_names = FALSE,
  show_row_dend = FALSE                   
)

collectri_net <- get_collectri(organism = 'human', split_complexes = FALSE)

tf_sirna_tstat <- run_ulm(mat = res_resistant[, 'stat', drop = FALSE], 
                          network = collectri_net, 
                          .source = "source", 
                          .target = "target", 
                          .mor = "mor", 
                          minsize = 5) %>%
  mutate(p_value_adj = p.adjust(p_value, method = "BH"))

sig_tf_sirna_tstat <- tf_sirna_tstat %>%
  filter(p_value_adj < 0.05) %>%
  filter(source %in% omnipath_filtered_nodes) %>%
  arrange(desc(abs(score))) %>%  
  slice_head(n = 50)

carnival_sirna_input <- sig_tf_sirna_tstat$score
names(carnival_sirna_input) <- sig_tf_sirna_tstat$source

t_stats_vector <- res_resistant$stat
names(t_stats_vector) <- rownames(res_resistant)
t_stats_matrix <- as.matrix(t_stats_vector)

my_progeny_weights <- progeny(t_stats_matrix, 
                              scale = FALSE, 
                              organism = "Human", 
                              top = 100)

raw_scores <- as.numeric(my_progeny_weights[1, ])
names(raw_scores) <- colnames(my_progeny_weights)

scaled_scores <- sapply(raw_scores, function(x) {
  if (x > 0) {
    return(x / max(raw_scores[raw_scores > 0]))
  } else if (x < 0) {
    return(-x / min(raw_scores[raw_scores < 0]))
  } else {
    return(0)
  }
})

pathway_df <- data.frame(
  source = names(scaled_scores),
  score = as.numeric(scaled_scores),
  stringsAsFactors = FALSE
)

expanded_weights_sirna <- pathway_df %>%
  filter(source %in% names(pathway_map)) %>%
  mutate(nodes = map(source, ~ pathway_map[[.x]])) %>%
  unnest(nodes) %>%
  filter(nodes %in% omnipath_filtered_nodes)

progeny_vector_sirna <- expanded_weights_sirna$score
names(progeny_vector_sirna) <- expanded_weights_sirna$nodes
progeny_vector_sirna <- progeny_vector_sirna[!duplicated(names(progeny_vector_sirna))]
#### TCGA Basal-Only GSVA & PTK2B Correlation Pipeline ####

expr <- read.table("cbio_tcga_all_expr.txt", header = TRUE, sep = "\t", check.names = FALSE)
clinical_data <- read.delim("data_clinical_patient.txt", header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)

cleaned_data <- expr %>%
  filter(!is.na(Hugo_Symbol) & Hugo_Symbol != "")

if (any(duplicated(cleaned_data$Hugo_Symbol))) {
  expr_means <- rowMeans(cleaned_data[, -c(1, 2)], na.rm = TRUE)
  cleaned_data <- cleaned_data %>%
    mutate(Row_Mean = expr_means) %>%
    group_by(Hugo_Symbol) %>%
    slice_max(order_by = Row_Mean, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(-Row_Mean)
}

cleaned_data <- cleaned_data %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  column_to_rownames(var = "Hugo_Symbol")

colnames(cleaned_data) <- sub("-01$", "", colnames(cleaned_data))

log_transformed_data <- log2(cleaned_data + 1)

basal_ids <- clinical_data[clinical_data$SUBTYPE == "BRCA_Basal", "PATIENT_ID"]
basal_ids <- intersect(basal_ids, colnames(log_transformed_data))
basal_expr <- as.matrix(log_transformed_data[, basal_ids])

msigdb_h <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark <- split(msigdb_h$gene_symbol, msigdb_h$gs_name)

set.seed(123)
basal_param <- gsvaParam(
  exprData = basal_expr,
  geneSets = hallmark,
  kcdf = "Gaussian", 
  minSize = 10,
  maxSize = 500
)

gsvsa_basal_matrix <- as.matrix(gsva(basal_param))
gsva_basal_z <- t(scale(t(gsvsa_basal_matrix))) 

clean_pathways_tcga <- gsub("^HALLMARK_", "", rownames(gsva_basal_z))
clean_pathways_tcga <- gsub("_", " ", clean_pathways_tcga)
clean_pathways_tcga <- str_to_title(clean_pathways_tcga)

acronyms <- c("Pi3k", "Akt", "Myc", "E2f", "Jak", "Stat3", "Tnfa", "Nfkb", 
              "Uv", "Tgf", "Kras", "Il6", "Il2", "Stat5", "Dna", "G2m", "Wnt", "V1", "V2", "Up", "Dn")

for (acr in acronyms) {
  clean_pathways_tcga <- gsub(paste0("\\b", acr, "\\b"), toupper(acr), clean_pathways_tcga)
}
clean_pathways_tcga <- gsub("\\bMtorc1\\b", "mTORC1", clean_pathways_tcga)
clean_pathways_tcga <- gsub("\\bMtor\\b", "mTOR", clean_pathways_tcga)
clean_pathways_tcga <- gsub("\\bP53\\b", "p53", clean_pathways_tcga)

rownames(gsva_basal_z) <- clean_pathways_tcga

ptk2b_basal_vec <- basal_expr["PTK2B", ]
basal_ptk2b_expr <- data.frame(
  Sample_ID = names(ptk2b_basal_vec),
  PTK2B = as.numeric(ptk2b_basal_vec)
)

ptk2b_basal_ord <- basal_ptk2b_expr[order(basal_ptk2b_expr$PTK2B), ]

gsva_basal_z_ordered <- gsva_basal_z[, ptk2b_basal_ord$Sample_ID]

cor_results <- apply(gsva_basal_z_ordered, 1, function(pathway_scores) {
  test <- cor.test(pathway_scores, ptk2b_basal_ord$PTK2B, method = "spearman", exact = FALSE)
  c(R = test$estimate[["rho"]], p = test$p.value)
})

cor_df <- as.data.frame(t(cor_results))
colnames(cor_df) <- c("R", "p_value")

cor_df$adj_p_value <- p.adjust(cor_df$p_value, method = "fdr")

cor_df$stars <- cut(cor_df$adj_p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"), right = FALSE)
cor_df$anno_text <- paste0(round(cor_df$R, 2), " ", cor_df$stars)

row_order_idx <- order(cor_df$R, decreasing = TRUE)
gsva_basal_z_sorted <- gsva_basal_z_ordered[row_order_idx, ] 
cor_df_sorted <- cor_df[row_order_idx, ]

presentation_pathways <- c(
  "IL6 JAK STAT3 Signaling",
  "Inflammatory Response",
  "TNFA Signaling Via NFKB",
  "Interferon Gamma Response",
  "PI3K AKT mTOR Signaling",
  "mTORC1 Signaling",
  "Epithelial Mesenchymal Transition",
  "Hypoxia",
  "MYC Targets V1",
  "E2F Targets"
)

kept_pathways <- intersect(rownames(gsva_basal_z_sorted), presentation_pathways)
gsva_basal_z_final <- gsva_basal_z_sorted[kept_pathways, ]
cor_df_final <- cor_df_sorted[kept_pathways, ]

heatmap_col_tcga <- colorRamp2(seq(-2.5, 2.5, length = 201), colorRampPalette(c("blue", "white", "red"))(201))

top_anno_tcga <- HeatmapAnnotation(
  PTK2B = anno_barplot(
    ptk2b_basal_ord$PTK2B, 
    bar_width = 1, 
    gp = gpar(col = NA, fill = "grey70"),
    border = FALSE,
    height = unit(2, "cm"),
    baseline = min(ptk2b_basal_ord$PTK2B) 
  ),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  annotation_name_rot = 90  
)

right_anno_tcga <- rowAnnotation(
  Stats = anno_text(
    cor_df_final$anno_text, 
    gp = gpar(fontsize = 12, fontface = "italic"), 
    just = "left"
  )
)

final_heatmap_tcga <- Heatmap(
  gsva_basal_z_final,
  name = "GSVA\n(Z-score)",
  col = heatmap_col_tcga,
  top_annotation = top_anno_tcga,
  right_annotation = right_anno_tcga,
  cluster_columns = FALSE, 
  cluster_rows = FALSE,    
  height = unit(7, "cm"), 
  show_column_names = FALSE, 
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"), 
  border = TRUE
)

draw(final_heatmap_tcga, heatmap_legend_side = "right", annotation_legend_side = "right")
#### CIBERSORTx Immune Deconvolution Pipeline ####

tcga_expr <- read.delim("TCGA+cibersort/data_mrna_seq_v2_rsem.txt")

p_clin <- read_tsv("TCGA+cibersort/data_clinical_patient.txt", comment = "#")
s_clin <- read_tsv("TCGA+cibersort/data_clinical_sample.txt", comment = "#")

colnames(tcga_expr) <- gsub("\\.", "-", gsub("\\.01$", "", colnames(tcga_expr)))

tcga_expr <- tcga_expr[!(tcga_expr$Hugo_Symbol %in% c("", " ", NA)), ]

numeric_cols <- sapply(tcga_expr, is.numeric)

row_means <- rowMeans(tcga_expr[, numeric_cols], na.rm = TRUE)

tcga_ordered <- tcga_expr[order(row_means, decreasing = TRUE), ]
tcga_clean <- tcga_ordered[!duplicated(tcga_ordered$Hugo_Symbol), ]

rownames(tcga_clean) <- tcga_clean$Hugo_Symbol 
tcga_clean$Entrez_Gene_Id <- NULL
tcga_clean$Hugo_Symbol<- NULL

basal_ids <- p_clin$PATIENT_ID[which(p_clin$SUBTYPE == "BRCA_Basal")]

basal_clin <- p_clin[p_clin$PATIENT_ID %in% basal_ids,]

ciber_res <- read.csv("TCGA+cibersort/pancan_BRCA_res.csv")

tnbc_ciber_master <- merge(ciber_res, basal_clin, by.x = "Mixture", by.y = "PATIENT_ID" )

target_genes <- c("PTK2B", "PTEN", "SRC", "AKT1")

genes_subset <- as.data.frame(t(tcga_clean[target_genes, basal_ids]))
colnames(genes_subset) <- target_genes
genes_subset$Mixture <- rownames(genes_subset)

tnbc_ciber_master <- merge(tnbc_ciber_master, genes_subset, by = "Mixture")

plot_data <- tnbc_ciber_master |>
  dplyr::select(Mixture, PTK2B, Macrophages.M1, T.cells.regulatory..Tregs., Absolute.score..sig.score.) |>
  mutate(PTK2B_log2 = log2(PTK2B + 1)) |>
  pivot_longer(
    cols = c(Macrophages.M1, T.cells.regulatory..Tregs., Absolute.score..sig.score.), 
    names_to = "Cell_Type",
    values_to = "Fraction"
  ) |>
  mutate(Cell_Type = factor(Cell_Type, 
                            levels = c("Macrophages.M1", "T.cells.regulatory..Tregs.", "Absolute.score..sig.score."),
                            labels = c("M1 Macrophages", "Regulatory T Cells", "All Immune Cells")))

correlation_plot <- ggplot(plot_data, aes(x = PTK2B_log2, y = Fraction, color = Cell_Type, fill = Cell_Type)) +
  geom_point(alpha = 0.5, size = 2) + 
  geom_smooth(method = "lm", alpha = 0.2) + 
  facet_wrap(~ Cell_Type, scales = "free_y") + 
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 4, color = "black") + 
  scale_color_manual(values = c("#E64B35", "#4DBBD5", "#77DD77")) + 
  scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#77DD77")) +
  theme_prism(base_size = 14) +
  labs(
    x = "Normalised PTK2B Expression",
    y = "CIBERSORT Immune Fraction"
  ) +
  theme(
    legend.position = "none", 
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

common_samples <- intersect(colnames(gsva_basal_z_ordered), ptk2b_basal_ord$Sample_ID)
common_samples <- intersect(common_samples, ciber_res$Mixture)

ptk2b_vec <- ptk2b_basal_ord$PTK2B[match(common_samples, ptk2b_basal_ord$Sample_ID)]
immune_vec <- ciber_res$Absolute.score..sig.score.[match(common_samples, ciber_res$Mixture)]
gsva_matched <- gsva_basal_z_ordered[, common_samples]

pcor_all_results <- apply(gsva_matched, 1, function(pathway_scores) {
  tryCatch({
    test <- pcor.test(ptk2b_vec, pathway_scores, immune_vec, method = "spearman")
    c(partial_R = test$estimate, p_value = test$p.value)
  }, error = function(e) {
    c(partial_R = NA, p_value = NA)
  })
})

pcor_all_df <- as.data.frame(t(pcor_all_results))
pcor_all_df$adj_p_value <- p.adjust(pcor_all_df$p_value, method = "fdr")
pcor_all_df$stars <- cut(pcor_all_df$adj_p_value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", "ns"), right = FALSE)
pcor_all_df <- pcor_all_df[order(pcor_all_df$partial_R, decreasing = TRUE), ]

#### scRNA-seq Pre-Processing (GSE176078) ####

base_dir <- "GSE_176078_new"
sample_folders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

seurat_list <- list()

for (folder in sample_folders) {
  
  patient_id <- basename(folder)
  
  mtx_file      <- file.path(folder, "count_matrix_sparse.mtx")
  barcodes_file <- file.path(folder, "count_matrix_barcodes.tsv")
  features_file <- file.path(folder, "count_matrix_genes.tsv") 
  meta_file     <- file.path(folder, "metadata.csv")
  
  counts <- ReadMtx(mtx = mtx_file, 
                    cells = barcodes_file, 
                    features = features_file, 
                    feature.column = 1)
  
  obj <- CreateSeuratObject(counts = counts, 
                            project = patient_id, 
                            min.cells = 3, 
                            min.features = 200)
  
  if (file.exists(meta_file)) {
    sample_meta <- read.csv(meta_file, row.names = 1) 
    obj <- AddMetaData(obj, metadata = sample_meta)
  }
  
  obj$Patient_ID <- patient_id
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  obj <- subset(obj, subset = nCount_RNA > 250)
  
  seurat_list[[patient_id]] <- obj
}

base_obj <- seurat_list[[1]]
rest_of_objs <- seurat_list[-1]

tnbc_merged <- merge(x = base_obj, 
                     y = rest_of_objs, 
                     add.cell.ids = names(seurat_list), 
                     project = "GSE176078_Atlas")
#### scRNA-seq Analysis & PROGENy Footprinting ####

seurat_raw <- tnbc_merged
seurat_filtered <- subset(seurat_raw, subset = percent.mt < 20)
seurat_filtered <- JoinLayers(seurat_filtered)

run_df_modern <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  
  sweep.res <- paramSweep(obj, PCs = 1:30, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  doublet_rate <- (ncol(obj) / 1000) * 0.008
  nExp <- round(doublet_rate * ncol(obj))
  
  obj <- doubletFinder(obj, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp, sct = FALSE)
  colnames(obj@meta.data)[grep("DF.classifications", colnames(obj@meta.data))] <- "Doublet_Status"
  return(obj)
}

seurat_list_df <- SplitObject(seurat_filtered, split.by = "Patient_ID")
seurat_list_df <- lapply(seurat_list_df, run_df_modern)

seurat_merged_df <- merge(seurat_list_df[[1]], y = seurat_list_df[-1])
seurat_merged_df <- JoinLayers(seurat_merged_df)
seurat_singlets <- subset(seurat_merged_df, subset = Doublet_Status == "Singlet")

seurat_final <- seurat_singlets

cancer_clusters <- c("Cancer Basal SC", "Cancer Cycling", "Cancer LumA SC", "Cancer LumB SC")

tumor_obj <- subset(seurat_final, 
                    subset = celltype_minor %in% cancer_clusters & 
                      Patient_ID != "CID44991" & 
                      celltype_minor != "Cancer Her2 SC")

tumor_obj$Patient_ID <- droplevels(as.factor(tumor_obj$Patient_ID))
tumor_obj$celltype_minor <- droplevels(as.factor(tumor_obj$celltype_minor))

tumor_obj <- NormalizeData(tumor_obj, verbose = FALSE)
tumor_obj <- FindVariableFeatures(tumor_obj, nfeatures = 2000, verbose = FALSE)
tumor_obj <- ScaleData(tumor_obj, verbose = FALSE)
tumor_obj <- RunPCA(tumor_obj, npcs = 30, verbose = FALSE)
tumor_obj <- RunHarmony(tumor_obj, group.by.vars = "Patient_ID", dims.use = 1:30)
tumor_obj <- RunUMAP(tumor_obj, reduction = "harmony", dims = 1:30)

umap_coords <- Embeddings(tumor_obj, reduction = "umap")
dbscan_res <- dbscan(umap_coords, eps = 0.55, minPts = 15)
tumor_obj$dbscan_clusters <- as.character(dbscan_res$cluster)

basal_final <- subset(tumor_obj, subset = dbscan_clusters != "0")
basal_final$dbscan_clusters <- droplevels(as.factor(basal_final$dbscan_clusters))
Idents(basal_final) <- "dbscan_clusters"

expr_matrix_sparse <- GetAssayData(basal_final, assay = "RNA", layer = "data")
expr_matrix_dense <- as.matrix(expr_matrix_sparse)

progeny_matrix <- progeny(expr_matrix_dense, scale=FALSE, organism="Human", top=500, perm=1)

basal_final[["progeny"]] <- CreateAssayObject(data = t(progeny_matrix))
basal_final <- ScaleData(basal_final, assay = "progeny")
progeny_scores <- t(GetAssayData(basal_final, assay = "progeny", layer = "scale.data"))

basal_final$PROGENy_PI3K <- progeny_scores[, "PI3K"]
basal_final$PROGENy_NFkB <- progeny_scores[, "NFkB"]
basal_final$PROGENy_JAKSTAT <- progeny_scores[, "JAK-STAT"]

ptk2b_expr <- FetchData(basal_final, vars = "PTK2B", layer = "data")
basal_final$PTK2B_Status <- ifelse(ptk2b_expr$PTK2B > 0, "PTK2B_Positive", "PTK2B_Negative")
basal_final$PTK2B_Status <- factor(basal_final$PTK2B_Status, levels = c("PTK2B_Negative", "PTK2B_Positive"))

sc_plot_data <- FetchData(basal_final, vars = c("PTK2B_Status", "PROGENy_PI3K", "PROGENy_NFkB", "PROGENy_JAKSTAT"))

tnbc_secretome <- list(c("IL6", "CXCL8", "CCL2", "VEGFA", "TNF", "IL1A", "IL1B", "CSF1"))
pyk2_empirical_secretome <- list(c("CCL2", "CXCL1", "CXCL8", "TNF", "S100A8", "S100A9"))

basal_final <- AddModuleScore(
  object = basal_final,
  features = tnbc_secretome,
  name = "Secretome_Score"
)

basal_final <- AddModuleScore(
  object = basal_final,
  features = pyk2_empirical_secretome,
  name = "PYK2_secretome"
)

secretome_df <- FetchData(basal_final, vars = "Secretome_Score1")
colnames(secretome_df) <- "Secretome_Score"
sc_plot_data$Secretome_Score <- secretome_df$Secretome_Score

pyk2_secretome_df <- FetchData(basal_final, vars = "PYK2_secretome1")
colnames(pyk2_secretome_df) <- "PYK2_secretome"
sc_plot_data$PYK2_secretome <- pyk2_secretome_df$PYK2_secretome

sc_comparisons <- list(c("PTK2B_Negative", "PTK2B_Positive"))
x_axis_mask <- scale_x_discrete(labels = c("PTK2B_Negative" = "PTK2B-Negative", "PTK2B_Positive" = "PTK2B-Positive"))

pA_sc <- DimPlot(seurat_final, reduction = "umap", group.by = "celltype_major", pt.size = 0.1) +
  theme_void() + 
  ggtitle(NULL) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)) +
  guides(color = guide_legend(override.aes = list(size = 3)))

pB_sc <- FeaturePlot(basal_final, features = "PTK2B", pt.size = 0.5, order = TRUE) +
  scale_color_viridis(option = "plasma") +
  theme_void() +
  ggtitle("PTK2B") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right")

base_violin_theme <- theme_classic() + 
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_pi3k_sc <- ggplot(sc_plot_data, aes(x = PTK2B_Status, y = PROGENy_PI3K, fill = PTK2B_Status)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("PTK2B_Negative" = "grey85", "PTK2B_Positive" = "darkred")) +
  base_violin_theme + x_axis_mask +
  labs(title = "PI3K", x = "", y = "PROGENy Score") +
  stat_compare_means(comparisons = sc_comparisons, method = "wilcox.test", label = "p.signif")

p_nfkb_sc <- ggplot(sc_plot_data, aes(x = PTK2B_Status, y = PROGENy_NFkB, fill = PTK2B_Status)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("PTK2B_Negative" = "grey85", "PTK2B_Positive" = "darkred")) +
  base_violin_theme + x_axis_mask +
  labs(title = "NF-kB", x = "", y = "") +
  stat_compare_means(comparisons = sc_comparisons, method = "wilcox.test", label = "p.signif")

p_jakstat_sc <- ggplot(sc_plot_data, aes(x = PTK2B_Status, y = PROGENy_JAKSTAT, fill = PTK2B_Status)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("PTK2B_Negative" = "grey85", "PTK2B_Positive" = "darkred")) +
  base_violin_theme + x_axis_mask +
  labs(title = "JAK-STAT", x = "", y = "") +
  stat_compare_means(comparisons = sc_comparisons, method = "wilcox.test", label = "p.signif")

pD_sc <- ggplot(sc_plot_data, aes(x = PTK2B_Status, y = PYK2_secretome, fill = PTK2B_Status)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("PTK2B_Negative" = "grey85", "PTK2B_Positive" = "darkorange")) +
  theme_classic() + x_axis_mask +
  labs(title = "siPTK2B-Derived Secretome Signature", x = "", y = "Module Score") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(comparisons = sc_comparisons, method = "wilcox.test", label = "p.signif")

layout_design_sc <- "
  AABB
  CDEF
"

final_figure_sc <- pA_sc + pB_sc + p_pi3k_sc + p_nfkb_sc + p_jakstat_sc + pD_sc + 
  plot_layout(design = layout_design_sc, heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold", size = 18, family = "sans"))
