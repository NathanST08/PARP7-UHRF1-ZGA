#### Figure 1  ####

mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
cmpd.vec<-c("Glycerophosphocholine","Citicoline","Betaine","Norvaline 2","beta-N-Methylamino-L-alanine","thiamine","Acetylcarnitine","NAD","Nicotinamide ribotide","Deoxyadenosine monophosphate","Cytidine 5-diphosphate","Creatine","Putrescine 2","Phosphocreatine","Taurine","1-Methylhistamine","Thiamine pyrophosphate","2-Aminooctanoic acid","methylsuccinic acid")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA)

library(ggplot2)
library(dplyr)

kegg_data <- data.frame(
  Pathway = c("Thiamine metabolism", 
              "Nicotinate and nicotinamide metabolism", 
              "Glycine, serine and threonine metabolism", 
              "Arginine and proline metabolism", 
              "Glycerophospholipid metabolism", 
              "Citrate cycle (TCA cycle)", 
              "Pyruvate metabolism", 
              "Glycolysis / Gluconeogenesis"),
  Total = c(7, 15, 32, 35, 36, 20, 22, 24),
  Hits = c(2, 2, 2, 2, 2, 1, 1, 1),
  Raw_P = c(0.000982, 0.00476, 0.0210, 0.0249, 0.0263, 0.136, 0.149, 0.161)
)

kegg_data$Rich_Factor <- kegg_data$Hits / kegg_data$Total

kegg_data <- kegg_data %>%
  arrange(desc(Raw_P)) %>% 
  mutate(Pathway = factor(Pathway, levels = Pathway))

p <- ggplot(kegg_data, aes(x = Rich_Factor, y = Pathway)) +

  geom_point(aes(size = Hits, color = Raw_P)) +

  scale_color_gradient(low = "#DC0000", high = "#4DBBD5", name = "Raw P-value") +

  scale_size_continuous(range = c(5, 10), name = "Metabolite Hits", 
                        breaks = c(1, 2)) + 

  theme_bw() +

  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),

    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),

    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor = element_blank(),

    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +

  labs(
    x = "Rich Factor",
    y = ""
  )

print(p)
ggsave("Metabolism_Enrichment_1E_Golden.pdf", plot = p, width = 8.5, height = 5.5, dpi = 300)

#### Figure 3  ####

mSet<-InitDataObjects("conc", "msetora", FALSE, 150)
cmpd.vec<-c("L-Lysine","Glutamate","D-Glucosamine-6-phosphate","NAD")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_2_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_2_", "png", 150, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_3_", "net", "png", 150, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_3_", "png", 150, width=NA)

library(ggplot2)
library(dplyr)

  Pathway = c("Nitrogen metabolism", 
              "Biotin metabolism", 
              "Arginine biosynthesis", 
              "Butanoate metabolism", 
              "Nicotinate and nicotinamide metabolism", 
              "Alanine, aspartate and glutamate metabolism", 
              "Glutathione metabolism", 
              "Glyoxylate and dicarboxylate metabolism", 
              "Arginine and proline metabolism")
  Total = c(6, 10, 14, 15, 15, 28, 28, 32, 35)
  Hits = c(1, 1, 1, 1, 1, 1, 1, 1, 1)
  Raw_P = c(0.0118, 0.0196, 0.0274, 0.0294, 0.0294, 0.0543, 0.0543, 0.0619, 0.0676)

kegg_data$Rich_Factor <- kegg_data$Hits / kegg_data$Total

kegg_data <- kegg_data %>%
  arrange(desc(Raw_P)) %>% 
  mutate(Pathway = factor(Pathway, levels = Pathway))

p <- ggplot(kegg_data, aes(x = Rich_Factor, y = Pathway)) +
  
  geom_point(aes(size = Hits, color = Raw_P)) +

  scale_color_gradient(low = "#DC0000", high = "#4DBBD5", name = "Raw P-value") +

  scale_size_continuous(limits = c(1, 2), range = c(6, 10), 
                        breaks = c(1), name = "Metabolite Hits") + 

  theme_bw() +

  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_line(color = "grey85", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  ) +

  labs(
    x = "Rich Factor",
    y = ""
  )

print(p)
ggsave("NAD_AminoAcid_Remodeling.pdf", plot = p, width = 9, height = 5.5, dpi = 300)

#### Figure 6  ####

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


### Figure 4 ####

options(repos = structure(c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))) 

library(readxl) 

library(dplyr) 

library(openxlsx) 

library(readxl) 

library(dplyr) 

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

cat("\n=== Analysis Completion Report ===\n") 

cat("✓ File:", file_path, "\n") 

cat("✓ Used gene name column:", gene_col, "\n") 

cat("✓ Number of analyzed genes:", nrow(expression_data), "\n") 

cat("✓ CTL group samples:", paste(ctl_samples, collapse = ", "), "\n") 

cat("✓ RBN group samples:", paste(rbn_samples, collapse = ", "), "\n") 

cat("✓ Generated heatmap files:\n") 

cat("   1. Up_minor_7_genes_heatmap.pdf - Standard heatmap\n") 

cat("   2. Up_minor_7_genes_heatmap_ordered.pdf - Ordered heatmap\n") 

cat("   3. Up_minor_7_genes_heatmap_ggplot.pdf - ggplot2 heatmap\n") 

cat("   4. Up_minor_7_genes_heatmap_with_values.pdf - Heatmap with values\n") 

cat("✓ Data files:\n") 

cat("   1. Up_minor_7_genes_normalized_data.csv - Normalized data\n") 

cat("   2. Up_minor_7_genes_log2_data.csv - log2 transformed data\n") 

cat("   3. Up_minor_7_genes_expression_summary.csv - Expression summary\n") 

cat("\nGene list:\n") 

for (i in 1:nrow(expression_data)) { 
  
  cat(i, ". ", rownames(expression_data)[i], "\n", sep = "") 
  
} 

cat("\n=== Debugging Information ===\n") 

cat("expression_data dimensions:", dim(expression_data), "\n") 

cat("expression_data_scaled dimensions:", dim(expression_data_scaled), "\n") 

cat("Row names example:", head(rownames(expression_data)), "\n") 

cat("Column names example:", head(colnames(expression_data)), "\n")


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

pdf("RBN-1_major-ZGA_genes_annotation_pie.pdf", width = 10, height = 6)  
print(pie) 

dev.off()
