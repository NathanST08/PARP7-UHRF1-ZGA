#### Figure 1 #### 

library(MetaboAnalystR) 
mSet <- InitDataObjects("conc", "msetora", FALSE, 150) 
cmpd.vec <- c("Glycerophosphocholine", "Citicoline", "Betaine", "Norvaline 2", 
              "beta-N-Methylamino-L-alanine", "thiamine", "Acetylcarnitine", 
              "NAD", "Nicotinamide ribotide", "Deoxyadenosine monophosphate", 
              "Cytidine 5-diphosphate", "Creatine", "Putrescine 2", 
              "Phosphocreatine", "Taurine", "1-Methylhistamine", 
              "Thiamine pyrophosphate", "2-Aminooctanoic acid", 
              "methylsuccinic acid", "L-Cystathionine", "S-Cysteinosuccinic acid") 

mSet <- Setup.MapData(mSet, cmpd.vec) 
mSet <- CrossReferencing(mSet, "name") 
mSet <- CreateMappingResultTable(mSet) 
mSet <- SetMetabolomeFilter(mSet, F) 
mSet <- SetCurrentMsetLib(mSet, "kegg_pathway", 2) 
mSet <- CalculateHyperScore(mSet) 

mSet <- PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA) 
mSet <- PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA) 
print("Pathway analysis completed successfully. Plots have been saved.") 

library(ggplot2) 
library(dplyr) 
library(forcats) 

data_text <- "Pathway\ttotal\texpected\thits\tRaw_p\tHolm_p\tFDR
Thiamine metabolism\t7\t0.0553\t2\t0.00118\t0.0953\t0.0666
Glycine, serine and threonine metabolism\t32\t0.253\t3\t0.00164\t0.132\t0.0666
Nicotinate and nicotinamide metabolism\t15\t0.118\t2\t0.00568\t0.448\t0.153
One carbon pool by folate\t26\t0.205\t2\t0.0167\t1\t0.339
Arginine and proline metabolism\t35\t0.276\t2\t0.0295\t1\t0.419
Glycerophospholipid metabolism\t36\t0.284\t2\t0.0311\t1\t0.419
Taurine and hypotaurine metabolism\t8\t0.0632\t1\t0.0616\t1\t0.713
Citrate cycle (TCA cycle)\t20\t0.158\t1\t0.148\t1\t1
Ether lipid metabolism\t20\t0.158\t1\t0.148\t1\t1
Pyruvate metabolism\t22\t0.174\t1\t0.161\t1\t1
Propanoate metabolism\t22\t0.174\t1\t0.161\t1\t1
Glycolysis / Gluconeogenesis\t24\t0.19\t1\t0.175\t1\t1
Lipoic acid metabolism\t28\t0.221\t1\t0.201\t1\t1
Lysine degradation\t29\t0.229\t1\t0.207\t1\t1
Cysteine and methionine metabolism\t33\t0.261\t1\t0.232\t1\t1
Valine, leucine and isoleucine degradation\t40\t0.316\t1\t0.275\t1\t1
Primary bile acid biosynthesis\t46\t0.363\t1\t0.31\t1\t1
Purine metabolism\t70\t0.553\t1\t0.433\t1\t1" 

df <- read.table(text = data_text, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 

target_pathways <- c( 
  "Glycerophospholipid metabolism", 
  "Taurine and hypotaurine metabolism", 
  "Nicotinate and nicotinamide metabolism",  
  "Ether lipid metabolism", 
  "Pyruvate metabolism",                    
  "Primary bile acid biosynthesis", 
  "Purine metabolism" 
) 

df_filtered <- df %>% 
  filter(Pathway %in% target_pathways) %>% 
  mutate(RichFactor = hits / total) %>% 
  arrange(RichFactor) %>% 
  mutate(Pathway = factor(Pathway, levels = Pathway)) 

y_axis_levels <- levels(df_filtered$Pathway) 
y_colors <- ifelse(y_axis_levels == "Nicotinate and nicotinamide metabolism", "red", "black") 

p <- ggplot(df_filtered, aes(x = RichFactor, y = Pathway)) +
  geom_point(aes(size = hits, color = FDR)) +
  scale_color_gradient(low = "#e41a1c", high = "#fee08b", name = "P-adjusted\n(FDR)") +
  scale_size_continuous(range = c(4, 9), breaks = c(1, 2, 3), name = "Gene Count\n(Hits)") +
  labs(x = "Rich Factor (Hits / Total Genes)", y = "") +
  theme_bw() +
  theme( 
    axis.text.y = element_text(size = 12, color = y_colors, face = "bold"), 
    axis.text.x = element_text(size = 11, color = "black"), 
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)), 
    legend.title = element_text(size = 11, face = "bold"), 
    legend.text = element_text(size = 10), 
    panel.grid.major.y = element_line(color = "grey90"), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1) 
  ) 

print(p) 
ggsave("Specific_Pathways_Bubble.pdf", plot = p, width = 7.5, height = 5.5, dpi = 300) 
print("Plot successfully generated and saved.") 

#### Figure 3 #### 

library(MetaboAnalystR) 
mSet <- InitDataObjects("conc", "msetora", FALSE, 150) 
cmpd.vec <- c("Lysine", "Glutamic Acid", "Glutamine", "Nicotinic acid", 
              "O-Phosphotyrosine", "Cytidine monophosphate(CMP)", "Xanthosine", 
              "Histamine", "Cytosine", "Imidazole", "Tyramine", "Methyladenosine", 
              "N-Acetylputrescine", "succinyl-CoA", "Agmatine", 
              "Nicotinamide Adenine Dinucleotide Hydrate(NAD)", 
              "Deoxycytidine 5-monophosphate", "Hypotaurine", "Adenine", 
              "Hydroxykynurenine", "Retinal", "Picolinic acid", "Alanine", 
              "5,6-dihydrothymine") 

mSet <- Setup.MapData(mSet, cmpd.vec) 
mSet <- CrossReferencing(mSet, "name") 
mSet <- CreateMappingResultTable(mSet) 
mSet <- SetMetabolomeFilter(mSet, F) 
mSet <- SetCurrentMsetLib(mSet, "kegg_pathway", 2) 
mSet <- CalculateHyperScore(mSet) 

mSet <- PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA) 
mSet <- PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA) 

mSet <- CalculateHyperScore(mSet) 
mSet <- PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA) 
mSet <- PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA) 
print("Pathway analysis completed successfully. Plots have been saved to the working directory.") 

library(ggplot2) 
library(dplyr) 
library(forcats) 

data_text <- "Pathway\ttotal\texpected\thits\tRaw_p\tHolm_p\tFDR
Urea Cycle\t28\t0.419\t4\t0.000541\t0.053\t0.0348
Glucose-Alanine Cycle\t13\t0.195\t3\t0.000711\t0.069\t0.0348
Nicotinate and Nicotinamide Metabolism\t35\t0.524\t4\t0.0013\t0.125\t0.0424
Purine Metabolism\t73\t1.09\t5\t0.00302\t0.287\t0.0741
Glutamate Metabolism\t48\t0.719\t4\t0.00429\t0.403\t0.0841
Warburg Effect\t57\t0.853\t4\t0.00803\t0.747\t0.093
Lysine Degradation\t30\t0.449\t3\t0.00865\t0.796\t0.093
Malate-Aspartate Shuttle\t10\t0.15\t2\t0.00879\t0.8\t0.093
Glycine and Serine Metabolism\t59\t0.883\t4\t0.00909\t0.818\t0.093
Ammonia Recycling\t31\t0.464\t3\t0.00949\t0.845\t0.093" 

df <- read.table(text = data_text, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 

df_processed <- df %>% 
  mutate(RichFactor = hits / total) %>% 
  arrange(RichFactor) %>% 
  mutate(Pathway = factor(Pathway, levels = Pathway)) 

y_axis_levels <- levels(df_processed$Pathway) 
y_colors <- ifelse(y_axis_levels == "Nicotinate and Nicotinamide Metabolism", "red", "black") 

p <- ggplot(df_processed, aes(x = RichFactor, y = Pathway)) +
  geom_point(aes(size = hits, color = FDR)) +
  scale_color_gradient(low = "#e41a1c", high = "#fee08b", name = "P-adjusted\n(FDR)") +
  scale_size_continuous(range = c(4, 9), breaks = c(2, 3, 4, 5), name = "Gene Count\n(Hits)") +
  labs(x = "Rich Factor (Hits / Total Genes)", y = "") +
  theme_bw() +
  theme( 
    axis.text.y = element_text(size = 12, color = y_colors, face = "bold"), 
    axis.text.x = element_text(size = 11, color = "black"), 
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)), 
    legend.title = element_text(size = 11, face = "bold"), 
    legend.text = element_text(size = 10), 
    panel.grid.major.y = element_line(color = "grey90"), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", size = 1) 
  ) 

print(p) 
ggsave("Updated_Metabolism_Bubble.pdf", plot = p, width = 7.5, height = 5.5, dpi = 300) 
print("Plot successfully generated.") 

#### Figure 4 #### 

options(repos = structure(c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) 
library(readxl) 
library(dplyr) 
library(openxlsx) 
library(tidyr) 
library(ggplot2) 
library(pheatmap) 
library(RColorBrewer) 

file_path <- "Down+Minor-FPKM-12.xlsx" 
data <- read_excel(file_path) 

cat("=== Basic Data Information ===\n") 
cat("Data dimensions:", dim(data), "\n") 
cat("Column names:\n") 
print(colnames(data)) 

if (!"gene_name" %in% colnames(data)) { 
  cat("Warning: 'gene_name' column not found. Available columns are:\n") 
  print(colnames(data)) 
  cat("Will use the first column as gene names.\n") 
  gene_col <- colnames(data)[1] 
} else { 
  gene_col <- "gene_name" 
} 

cat("Used gene name column:", gene_col, "\n") 

expression_data <- data[, !colnames(data) %in% gene_col]  
rownames(expression_data) <- data[[gene_col]]  

cat("Gene row names:\n") 
print(rownames(expression_data)) 
cat("Sample column names:\n") 
print(colnames(expression_data)) 

sample_names <- colnames(expression_data) 
ctl_samples <- sample_names[grepl("CTL", sample_names, ignore.case = TRUE)] 
rbn_samples <- sample_names[grepl("RBN", sample_names, ignore.case = TRUE)] 

cat("CTL samples:", ctl_samples, "\n") 
cat("RBN samples:", rbn_samples, "\n") 
cat("Number of genes:", nrow(expression_data), "\n") 

expression_data_log <- log2(expression_data + 1) 
expression_data_scaled <- t(scale(t(expression_data_log))) 

cat("Summary of scaled data:\n") 
print(summary(as.vector(expression_data_scaled))) 

annotation_col <- data.frame( 
  Group = c(rep("CTL", length(ctl_samples)), rep("RBN", length(rbn_samples))), 
  row.names = c(ctl_samples, rbn_samples) 
) 

annotation_colors <- list( 
  Group = c(CTL = "#1f77b4", RBN = "#ff7f0e") 
) 

color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100) 

pheatmap( 
  expression_data_scaled, 
  main = "Down-regulated Minor ZGA Genes (12 genes)\nZ-score normalized log2(FPKM+1)", 
  color = color_palette, 
  scale = "none", 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  clustering_distance_rows = "euclidean", 
  clustering_method = "complete", 
  show_rownames = TRUE, 
  show_colnames = TRUE, 
  annotation_col = annotation_col, 
  annotation_colors = annotation_colors, 
  fontsize_row = 12, 
  fontsize_col = 10, 
  fontsize = 11, 
  angle_col = 45, 
  border_color = "gray60", 
  treeheight_row = 20, 
  treeheight_col = 20, 
  cellwidth = 30, 
  cellheight = 20, 
  display_numbers = FALSE, 
  legend = TRUE, 
  filename = "Down_minor_12_genes_heatmap.pdf", 
  width = 8, 
  height = 6
) 

column_order <- c(ctl_samples, rbn_samples) 
expression_data_ordered <- expression_data_scaled[, column_order, drop = FALSE] 

pheatmap( 
  expression_data_ordered, 
  main = "Down-regulated Minor ZGA Genes Expression\nZ-score normalized log2(FPKM+1)", 
  color = color_palette, 
  scale = "none", 
  cluster_rows = TRUE, 
  cluster_cols = FALSE,  
  show_rownames = TRUE, 
  show_colnames = TRUE, 
  annotation_col = annotation_col[column_order, , drop = FALSE], 
  annotation_colors = annotation_colors, 
  fontsize_row = 12, 
  fontsize_col = 10, 
  angle_col = 45, 
  border_color = "gray60", 
  gaps_col = length(ctl_samples),  
  filename = "Down_minor_12_genes_heatmap_ordered.pdf", 
  width = 8, 
  height = 6
) 

library(reshape2) 

heatmap_data <- expression_data_ordered %>% 
  as.data.frame() %>% 
  mutate(Gene = rownames(.)) %>% 
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Zscore") %>% 
  mutate( 
    Group = ifelse(grepl("RBN", Sample), "RBN", "CTL"), 
    Sample = factor(Sample, levels = column_order)  
  ) 

p <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = Zscore)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2( 
    low = "#2166AC", 
    mid = "#F7F7F7", 
    high = "#B2182B", 
    midpoint = 0, 
    name = "Z-score\nlog2(FPKM+1)", 
    limits = c(-2, 2)  
  ) +
  labs( 
    title = "Down-regulated Minor ZGA Genes (12 genes)", 
    subtitle = "Z-score normalized log2(FPKM+1) values", 
    x = "", 
    y = "" 
  ) +
  theme_minimal() +
  theme( 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(size = 10, hjust = 0.5), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"), 
    axis.text.y = element_text(size = 11, face = "bold"), 
    axis.title = element_text(size = 12), 
    legend.title = element_text(size = 9), 
    legend.text = element_text(size = 8), 
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "gray50", fill = NA, linewidth = 0.5) 
  ) +
  geom_vline(xintercept = length(ctl_samples) + 0.5, linetype = "dashed", 
             color = "black", linewidth = 0.8) 

ggsave("Down_minor_12_genes_heatmap_ggplot.pdf", p, width = 10, height = 6) 

pheatmap( 
  expression_data_ordered, 
  main = "Down-regulated Minor ZGA Genes\nwith Expression Values", 
  color = color_palette, 
  scale = "none", 
  cluster_rows = TRUE, 
  cluster_cols = FALSE, 
  show_rownames = TRUE, 
  show_colnames = TRUE, 
  annotation_col = annotation_col[column_order, , drop = FALSE], 
  annotation_colors = annotation_colors, 
  fontsize_row = 11, 
  fontsize_col = 10, 
  fontsize = 10, 
  angle_col = 45, 
  border_color = "gray60", 
  display_numbers = TRUE, 
  number_format = "%.1f", 
  number_color = "black", 
  fontsize_number = 8, 
  gaps_col = length(ctl_samples), 
  filename = "Down_minor_12_genes_heatmap_with_values.pdf", 
  width = 9, 
  height = 6
) 

write.csv(expression_data_scaled, "Down_minor_12_genes_normalized_data.csv") 
write.csv(expression_data_log, "Down_minor_12_genes_log2_data.csv") 

cat("\n=== Data Summary ===\n") 
cat("Gene list:\n") 
print(rownames(expression_data)) 

mean_expression <- data.frame( 
  Gene = rownames(expression_data), 
  CTL_mean = rowMeans(expression_data[, ctl_samples, drop = FALSE], na.rm = TRUE), 
  RBN_mean = rowMeans(expression_data[, rbn_samples, drop = FALSE], na.rm = TRUE) 
) %>% 
  mutate( 
    Log2FC = log2((RBN_mean + 1) / (CTL_mean + 1))  
  ) 

print(mean_expression) 
write.csv(mean_expression, "Down_minor_12_genes_expression_summary.csv", row.names = FALSE) 

#### Figure 5 #### 

library(ChIPseeker) 
library(TxDb.Mmusculus.UCSC.mm10.knownGene) 
library(org.Mm.eg.db) 
library(clusterProfiler)  

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
target_genes <- read.table("ZGA_gene_list_1171.txt", header = FALSE)$V1  
target_genes_entrez <- bitr(target_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                            OrgDb = org.Mm.eg.db)$ENTREZID  

peaks <- readPeakFile("CTL-1_peaks.narrowPeak", header = FALSE) 
peak_anno <- annotatePeak( 
  peaks, 
  tssRegion = c(-3000, 3000), 
  TxDb = txdb, 
  annoDb = "org.Mm.eg.db" 
) 

anno_df <- as.data.frame(peak_anno)  
target_anno <- anno_df[anno_df$geneId %in% target_genes_entrez, ]  

anno_summary <- as.data.frame(table(gsub(" \\(.*", "", target_anno$annotation)))  
colnames(anno_summary) <- c("Category", "Count") 
anno_summary$Percent <- round(anno_summary$Count / sum(anno_summary$Count) * 100, 1) 
legend_labels <- paste0(anno_summary$Category, " (", anno_summary$Percent, "%)")  

library(ggplot2) 

pie <- ggplot(anno_summary, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_manual( 
    values = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3"), 
    labels = legend_labels  
  ) +
  ggtitle("CUT&Tag Genomic Annotation (Target Genes)") +
  theme( 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "right",  
    legend.text = element_text(size = 10),  
    legend.spacing.y = unit(0.2, "cm")  
  ) 

pdf("CTL-1_major-ZGA_genes_annotation_pie.pdf", width = 10, height = 6)  
print(pie) 
dev.off() 

peaks <- readPeakFile("RBN-1_peaks.narrowPeak", header = FALSE) 
peak_anno <- annotatePeak( 
  peaks, 
  tssRegion = c(-3000, 3000), 
  TxDb = txdb, 
  annoDb = "org.Mm.eg.db" 
) 

anno_df <- as.data.frame(peak_anno)  
target_anno <- anno_df[anno_df$geneId %in% target_genes_entrez, ]  

anno_summary <- as.data.frame(table(gsub(" \\(.*", "", target_anno$annotation)))  
colnames(anno_summary) <- c("Category", "Count") 
anno_summary$Percent <- round(anno_summary$Count / sum(anno_summary$Count) * 100, 1) 
legend_labels <- paste0(anno_summary$Category, " (", anno_summary$Percent, "%)")  

pie <- ggplot(anno_summary, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_manual( 
    values = c("#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3"), 
    labels = legend_labels  
  ) +
  ggtitle("CUT&Tag Genomic Annotation (Target Genes)") +
  theme( 
    plot.title = element_text(hjust = 0.5), 
    legend.position = "right",  
    legend.text = element_text(size = 10),  
    legend.spacing.y = unit(0.2, "cm")  
  ) 

pdf("RBN-1_major-ZGA_genes_annotation_pie.pdf", width = 10, height = 6)  
print(pie) 
dev.off() 

#### Figure 6 #### 

library(ggplot2) 
library(dplyr) 

file_path <- file.choose() 

if (grepl("\\.csv$", file_path, ignore.case = TRUE)) { 
  df <- read.csv(file_path, stringsAsFactors = FALSE) 
} else { 
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, quote = "") 
} 

if("Category" %in% colnames(df)) { 
  df$Category <- ifelse(df$Category == "Reactome Gene Sets", "Reactome", df$Category) 
} 

p_col <- grep("LogP|P.*val", colnames(df), value = TRUE, ignore.case = TRUE)[1] 
if (is.na(p_col)) { 
  stop("P-value or LogP column not found") 
} 

p_vals <- as.numeric(df[[p_col]]) 

if (all(p_vals <= 0, na.rm = TRUE)) { 
  raw_p <- 10^(p_vals) 
} else { 
  raw_p <- p_vals
} 

df$adj_p <- p.adjust(raw_p, method = "BH") 
df$minus_log_q <- -log10(df$adj_p) 

if ("InTerm_InList" %in% colnames(df)) { 
  df$GeneCount <- as.numeric(sub("/.*", "", df$InTerm_InList)) 
} else if ("Symbols" %in% colnames(df)) { 
  df$GeneCount <- sapply(strsplit(as.character(df$Symbols), ","), length) 
} else { 
  df$GeneCount <- 50
} 

df <- df %>% 
  arrange(Category, minus_log_q) %>% 
  mutate(Description = factor(Description, levels = unique(Description))) 

p <- ggplot(df, aes(x = minus_log_q, y = Description)) +
  geom_point(aes(size = GeneCount, color = minus_log_q)) +
  scale_color_gradient( 
    low = "#377eb8", 
    high = "#e41a1c", 
    name = expression("-log"[10]*"(q-value)") 
  ) +
  scale_size_continuous( 
    name = "Gene Count", 
    breaks = c(50, 100, 150), 
    range = c(3, 10) 
  ) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme( 
    strip.background = element_rect(fill = "#d9ead3", color = "black"), 
    strip.text.y = element_text(angle = 270, size = 11, face = "bold"), 
    axis.title.x = element_blank(), 
    axis.text.y = element_text(size = 11, color = "black"), 
    axis.text.x = element_text(size = 11, color = "black"), 
    panel.grid.major = element_line(color = "grey90"), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA) 
  ) 

print(p) 
ggsave(filename = "Adjusted_P_Pathway.pdf", plot = p, width = 10, height = 6) 
ggsave(filename = "Adjusted_P_Pathway.tiff", plot = p, width = 10, height = 6, dpi = 600, compression = "lzw")

#### Figure 7 #### 

library(readxl) 
library(pheatmap) 
library(RColorBrewer) 
library(viridis) 
library(ggplot2) 
library(grid) 
library(ComplexHeatmap) 
library(circlize) 

gene_data <- read_excel("extracted_gene list-15-作图.xlsx", sheet = "Sheet1") 

expression_matrix <- as.matrix(gene_data[, -1]) 
rownames(expression_matrix) <- gene_data$gene_name
colnames(expression_matrix) <- c("CTL3", "CTL4", "RBN3", "RBN4", "UHRF1-1", "UHRF1-2") 

expression_matrix_scaled <- t(scale(t(expression_matrix))) 

sample_groups <- data.frame( 
  Group = factor(rep(c("CTL", "RBN", "UHRF1"), each = 2), 
                 levels = c("CTL", "RBN", "UHRF1")), 
  row.names = colnames(expression_matrix) 
) 

main_palette <- viridis(100) 

group_colors <- list( 
  Group = c( 
    CTL = "#1B9E77",  
    RBN = "#D95F02",  
    UHRF1 = "#7570B3" 
  ) 
) 

col_annotation <- HeatmapAnnotation( 
  Group = sample_groups$Group, 
  col = list(Group = c("CTL" = "#1B9E77", "RBN" = "#D95F02", "UHRF1" = "#7570B3")), 
  annotation_name_side = "left", 
  annotation_legend_param = list( 
    title = "Treatment Group", 
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 9) 
  ), 
  gap = unit(1, "mm"), 
  simple_anno_size = unit(0.4, "cm") 
) 

row_means <- rowMeans(expression_matrix_scaled, na.rm = TRUE) 
row_annotation <- rowAnnotation( 
  "Mean Expression" = row_means, 
  col = list("Mean Expression" = colorRamp2( 
    c(min(row_means), 0, max(row_means)), 
    c("blue", "white", "red") 
  )), 
  annotation_name_side = "top", 
  annotation_legend_param = list( 
    title = "Mean Z-score", 
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 9) 
  ), 
  simple_anno_size = unit(0.3, "cm") 
) 

ht <- Heatmap( 
  expression_matrix_scaled, 
  name = "Expression Z-score", 
  col = colorRamp2(seq(-2, 2, length.out = 100), viridis(100)), 
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 9, fontface = "italic"), 
  column_names_side = "top", 
  column_names_gp = gpar(fontsize = 10, fontface = "bold"), 
  column_names_rot = 45, 
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  rect_gp = gpar(col = "white", lwd = 0.5), 
  cell_fun = function(j, i, x, y, width, height, fill) { 
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(fill = fill, col = "white")) 
  }, 
  width = unit(6, "cm"), 
  height = unit(10, "cm"), 
  heatmap_legend_param = list( 
    title = "Z-score", 
    title_gp = gpar(fontsize = 10, fontface = "bold"), 
    labels_gp = gpar(fontsize = 9), 
    legend_height = unit(3, "cm"), 
    grid_width = unit(0.4, "cm"), 
    at = c(-2, -1, 0, 1, 2) 
  ), 
  top_annotation = col_annotation, 
  right_annotation = row_annotation
) 

pdf("Heatmap_SCI_Quality.pdf", width = 8, height = 10, paper = "a4r") 
draw(ht, 
     padding = unit(c(2, 2, 2, 2), "cm"), 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right", 
     merge_legend = TRUE) 
dev.off() 

png("Heatmap_SCI_Quality.png", width = 2400, height = 3000, res = 300) 
draw(ht, 
     padding = unit(c(2, 2, 2, 2), "cm"), 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right", 
     merge_legend = TRUE) 
dev.off() 

pdf("Heatmap_pheatmap_optimized.pdf", width = 8, height = 10) 
pheatmap( 
  expression_matrix_scaled, 
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), 
  border_color = "gray80", 
  border = TRUE, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  annotation_col = sample_groups, 
  annotation_colors = group_colors, 
  annotation_names_col = TRUE, 
  annotation_legend = TRUE, 
  annotation_names_row = TRUE, 
  fontsize = 10, 
  fontsize_row = 9, 
  fontsize_col = 10, 
  fontface = "bold", 
  angle_col = 45, 
  display_numbers = FALSE, 
  number_color = "black", 
  number_format = "%.1f", 
  fontsize_number = 6, 
  cellwidth = 20, 
  cellheight = 12, 
  gaps_row = NULL, 
  gaps_col = c(2, 4), 
  legend = TRUE, 
  legend_breaks = c(-2, -1, 0, 1, 2), 
  legend_labels = c("-2.0", "-1.0", "0.0", "1.0", "2.0"), 
  main = "Gene Expression Heatmap", 
  fontsize_main = 12, 
  treeheight_row = 30, 
  treeheight_col = 30, 
  cutree_rows = NA, 
  cutree_cols = NA, 
  margins = c(8, 8), 
  silent = FALSE
) 
dev.off() 

library(reshape2) 

heatmap_data <- melt(expression_matrix_scaled) 
colnames(heatmap_data) <- c("Gene", "Sample", "Zscore") 

heatmap_data$Group <- rep(sample_groups$Group, each = nrow(expression_matrix_scaled)) 

ggplot_heatmap <- ggplot(heatmap_data, aes(x = Sample, y = Gene, fill = Zscore)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2( 
    low = "#2166AC", 
    mid = "white", 
    high = "#B2182B", 
    midpoint = 0, 
    name = "Z-score", 
    limits = c(-2, 2) 
  ) +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_minimal(base_size = 11) +
  theme( 
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), 
    axis.text.y = element_text(face = "italic"), 
    axis.title = element_blank(), 
    strip.text = element_text(face = "bold", size = 10), 
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "gray30", fill = NA, size = 0.5), 
    legend.position = "right", 
    legend.key.height = unit(2, "cm"), 
    plot.margin = unit(c(1, 1, 1, 1), "cm") 
  ) +
  labs( 
    title = "Gene Expression Heatmap", 
    subtitle = "CTL vs RBN vs UHRF1 Treatment Groups" 
  ) 

ggsave("Heatmap_ggplot2_version.pdf", ggplot_heatmap, width = 10, height = 8, dpi = 300) 
ggsave("Heatmap_ggplot2_version.png", ggplot_heatmap, width = 10, height = 8, dpi = 300) 
