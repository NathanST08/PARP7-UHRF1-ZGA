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
              "Arginine and proline metabolism"),
  Total = c(6, 10, 14, 15, 15, 28, 28, 32, 35),
  Hits = c(1, 1, 1, 1, 1, 1, 1, 1, 1),
  Raw_P = c(0.0118, 0.0196, 0.0274, 0.0294, 0.0294, 0.0543, 0.0543, 0.0619, 0.0676)
)

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
