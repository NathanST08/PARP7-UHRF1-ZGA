#### Figure 1  ####
#### Workflow for obtaining the metabolite expression matrix

rm(list=ls())

library(openxlsx)
library(ropls)

Experimental_Group=read.table('2C.txt',header = T,sep = '\t')
Control_Group=read.table('1C.txt',header = T,sep = '\t')
Blank_Group=read.table('Ctrl.txt',header = T,sep = '\t')

SN_calculation_method='SN'
ND=48

Cell_eg<-Experimental_Group
eg_name<-t(Cell_eg[1,-1])
eg_name<-data.frame(table(eg_name))
eg_name<-as.character(eg_name[1,1])

Cell_cg<-Control_Group
cg_name<-t(Cell_cg[1,-1])
cg_name<-data.frame(table(cg_name))
cg_name<-as.character(cg_name[1,1])

group.name<-paste(eg_name,'vs',cg_name)
dir.create(group.name, recursive = TRUE)

Cell_DATA_data<-merge(x=Cell_eg[-1,],y=Cell_cg[-1,],by.x = "X",by.y = "X",all.x = T,all.y = T)
Cell_DATA_group<-merge(x=Cell_eg[1,],y=Cell_cg[1,],by.x = "X",by.y = "X",all.x = T,all.y = T)
Cell_DATA<-rbind(Cell_DATA_group,Cell_DATA_data)
rownames(Cell_DATA) <-Cell_DATA[,1]
rm('Cell_DATA_group','Cell_DATA_data',eg_name,cg_name)

Cell_eg_name=variable.names(Cell_eg)[-1]
Cell_eg_name=c('mz',Cell_eg_name)
Cell_eg_num=as.numeric(length(Cell_eg_name))-1

Cell_cg_name=variable.names(Cell_cg)[-1]
Cell_cg_name=c('mz',Cell_cg_name)
Cell_cg_num=as.numeric(length(Cell_cg_name))-1

Cell_num=Cell_eg_num+Cell_cg_num
Cell_DATA[is.na(Cell_DATA)]<-0

Ctrl<-Blank_Group
DATA_data<-merge(x=Cell_DATA[-1,],y=Ctrl[-1,],by.x = "X",by.y = "X",all.x = T,all.y = T)
DATA_group<-merge(x=Cell_DATA[1,],y=Ctrl[1,],by.x = "X",by.y = "X",all.x = T,all.y = T)
DATA<-rbind(DATA_group,DATA_data)
rm('DATA_group','DATA_data')

Ctrl_name=variable.names(Ctrl)[-1]
Ctrl_name=c('mz',Ctrl_name)
Ctrl_num=as.numeric(length(Ctrl_name))-1
DATA[is.na(DATA)]<-0

AP01_AL01_Cell=DATA[,c(1,2:(1+Cell_num))]
AP01_AL01_Ctrl=DATA[,c(1,(2+Cell_num):(1+Ctrl_num+Cell_num))]
rm('Cell_eg_name','Cell_cg_name','Cell_num','Ctrl_name','Ctrl_num','DATA','Cell_DATA')

AP02_GP01_Cell<-AP01_AL01_Cell
rownames(AP02_GP01_Cell) <-AP02_GP01_Cell[,1]
data=AP02_GP01_Cell[-1,-1]
x=apply(data,c(1,2),as.numeric)
mz_cell=AP02_GP01_Cell[,1]
num_cell=AP02_GP01_Cell[1,-1]

e=x-ND
f=as.data.frame(e)
g=rbind(num_cell,f)
AP03_DN01_Cell=cbind(mz_cell,g)
xlsx.name=paste0(group.name,'/',group.name," AP03_DN01_Cell.xlsx")
write.xlsx(AP03_DN01_Cell,xlsx.name)
rm(list = c('data','x','e','f','g','mz_cell','num_cell'))

AP02_GP01_Ctrl<-AP01_AL01_Ctrl
rownames(AP02_GP01_Ctrl) <-AP02_GP01_Ctrl[,1]
data=AP02_GP01_Ctrl[-1,-1]
x=apply(data,c(1,2),as.numeric)
mz_ctrl=AP02_GP01_Ctrl[,1]
num_ctrl=AP02_GP01_Ctrl[1,-1]

e=x-ND
f=as.data.frame(e)
g=rbind(num_ctrl,f)
AP03_DN01_Ctrl=cbind(mz_ctrl,g)
xlsx.name=paste0(group.name,'/',group.name," AP03_DN01_Ctrl.xlsx")
write.xlsx(AP03_DN01_Ctrl,xlsx.name)
rm(list = c('data','x','e','f','g','mz_ctrl','num_ctrl'))

mz_Cell=AP03_DN01_Cell[,1]
num_Cell=AP03_DN01_Cell[1,-1]

N=AP03_DN01_Ctrl[-1,-1]
S=AP03_DN01_Cell[-1,-1]

n=apply(N,c(1,2),as.numeric)
mean_result=apply(n,1,mean)
mean_result<-ifelse(mean_result!=0,mean_result,1)
s=apply(S,c(1,2),as.numeric)

SN=s/mean_result

f=as.data.frame(SN)
g=rbind(num_Cell,f)
AP05_SN01_Cell_SN=cbind(mz_Cell,g)
AP05_SN01_Cell_SN=AP05_SN01_Cell_SN[!rownames(AP05_SN01_Cell_SN) %in% c("Cl-PHE_(200.0 -> 154.0)"), ]
xlsx.name=paste0(group.name,'/',group.name,' AP05_SN01_Cell_SN.xlsx')
write.xlsx(AP05_SN01_Cell_SN,xlsx.name)
rm(list=c('N','n','mean_result','S','s','f','g','mz_Cell','num_Cell','SN'))

rpols_data<-AP05_SN01_Cell_SN
data<-rpols_data[-1,-1]
data.t<-data.frame(t(data))
group<-t(rpols_data[1,-1])
e=apply(data.t,c(1,2),as.numeric)
rm(list = c('data.t','data'))

e.plsda <- opls(e, group, predI = 5, orthoI = 0)
AP07_DR02_Cell_PLSDA_RAW <- e.plsda

e.oplsda <- opls(e, group, predI=1, orthoI = NA)
AP07_DR03_Cell_OPLSDA_RAW <- e.oplsda
oplsda_modle1 <- AP07_DR03_Cell_OPLSDA_RAW@modelDF

rm(list=c('e','group','e.plsda','e.oplsda'))

data=AP05_SN01_Cell_SN[-1,-1]
EX_num=Cell_eg_num
CL_num=Cell_cg_num

Pvalue<-c(rep(0,nrow(data)))
for(i in 1:nrow(data)){
  data.test=t.test(as.numeric(data[i,1:EX_num]),as.numeric(data[i,(EX_num+1):(CL_num+EX_num)]))
  Pvalue[i]<-data.test$p.value
}

log2_FC<-c(rep(0,nrow(data)))
for(i in 1:nrow(data)){
  log2_FC[i]<-log2((mean(as.numeric(data[i,1:EX_num])))/(mean(as.numeric(data[i,(EX_num+1):(CL_num+EX_num)]))))
}

demo_p_FC<-as.data.frame(rownames(data))
colnames(demo_p_FC)<-('mz')
demo_p_FC$Pvalue=Pvalue
demo_p_FC$log2_FC=log2_FC

if (nrow(oplsda_modle1) == 0 && ncol(oplsda_modle1) == 0) {
  data_VIP=as.data.frame(AP07_DR02_Cell_PLSDA_RAW@vipVn)
} else {
  data_VIP=as.data.frame(AP07_DR03_Cell_OPLSDA_RAW@vipVn)
}
colnames(data_VIP)<-c('VIP')
rownames(data_VIP)<-c(demo_p_FC[,1])
data_VIP=cbind(rownames(data_VIP),data_VIP)
data_VIP_rowname=data_VIP$`rownames(data_VIP)`
mz=data_VIP_rowname=gsub('[X]', '', data_VIP_rowname)
data_VIP=cbind(mz,data_VIP)
data_VIP=data_VIP[,-2]

demo_p_FC_VIP<-merge(x=demo_p_FC,y=data_VIP,by.x = "mz",by.y = "mz",all.x = T,all.y = T)
demo_p_FC_VIP[is.na(demo_p_FC_VIP)]<-0

cut_off_pvalue = 0.05
cut_off_logFC = 1
demo_p_FC_VIP$Change = ifelse(demo_p_FC_VIP$Pvalue < cut_off_pvalue & abs(demo_p_FC_VIP$log2_FC) >= cut_off_logFC,
                              ifelse(demo_p_FC_VIP$log2_FC> cut_off_logFC ,'Up','Down'),'Stable')

AP08_SD01_all_metabolites<-demo_p_FC_VIP
AP08_SD01_differential_metabolites<-AP08_SD01_all_metabolites[AP08_SD01_all_metabolites$Pvalue<0.05,]
AP08_SD01_differential_metabolites<-AP08_SD01_differential_metabolites[abs(AP08_SD01_differential_metabolites$log2_FC)>1,]
AP08_SD01_differential_metabolites<-AP08_SD01_differential_metabolites[AP08_SD01_differential_metabolites$VIP>1,]

xlsx.name=paste0(group.name,'/',group.name," AP08_SD01_differential_metabolites.xlsx")
write.xlsx(AP08_SD01_differential_metabolites, file = xlsx.name)

rm(list = c('data','i','Pvalue','log2_FC','demo_p_FC','data_VIP',
            'mz','data.test','demo_p_FC_VIP','data_VIP_rowname','cut_off_pvalue','cut_off_logFC',
            'AP08_SD01_all_metabolites','AP08_SD01_differential_metabolites'))

print("Analysis completed. The specified tables have been successfully exported to the corresponding folders.")


#### Screen of 17 significantly downregulated metabolites (Padj < 0.05) in 2-cell vs. 1-cell.
#### KEGG pathway enrichment analysis for distinct metabolic clusters was performed using the MetaboAnalyst web server (version 6.0.0, https://www.metaboanalyst.ca), which operates via the underlying MetaboAnalystR package (version 4.0)”
mSet<-InitDataObjects("conc", "msetora", FALSE, 150) 
cmpd.vec<-c("Glycerophosphocholine","Citicoline","Betaine","Norvaline 2","beta-N-Methylamino-L-alanine","thiamine","Acetylcarnitine","NAD","Nicotinamide ribotide","Deoxyadenosine monophosphate","Cytidine 5-diphosphate","Creatine","Putrescine 2","Phosphocreatine","Taurine","1-Methylhistamine","Thiamine pyrophosphate") 
mSet<-Setup.MapData(mSet, cmpd.vec); 
mSet<-CrossReferencing(mSet, "name"); 
mSet<-CreateMappingResultTable(mSet) 
mSet<-PerformDetailMatch(mSet, "Norvaline 2"); 
mSet<-GetCandidateList(mSet); 
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
library(forcats) 

data_text <- "Pathway\ttotal\texpected\thits\tRaw p\tHolm p\tFDR
Thiamine metabolism\t7\t0.0507\t2\t0.000982\t0.0796\t0.0796
Nicotinate and nicotinamide metabolism\t15\t0.109\t2\t0.00476\t0.381\t0.193
Glycine, serine and threonine metabolism\t32\t0.232\t2\t0.021\t1\t0.426
Arginine and proline metabolism\t35\t0.253\t2\t0.0249\t1\t0.426
Glycerophospholipid metabolism\t36\t0.261\t2\t0.0263\t1\t0.426
Taurine and hypotaurine metabolism\t8\t0.0579\t1\t0.0566\t1\t0.764
Citrate cycle (TCA cycle)\t20\t0.145\t1\t0.136\t1\t1
Ether lipid metabolism\t20\t0.145\t1\t0.136\t1\t1
Pyruvate metabolism\t22\t0.159\t1\t0.149\t1\t1
Propanoate metabolism\t22\t0.159\t1\t0.149\t1\t1
Glycolysis / Gluconeogenesis\t24\t0.174\t1\t0.161\t1\t1
One carbon pool by folate\t26\t0.188\t1\t0.173\t1\t1
Lipoic acid metabolism\t28\t0.203\t1\t0.186\t1\t1
Lysine degradation\t29\t0.21\t1\t0.192\t1\t1
Valine, leucine and isoleucine degradation\t40\t0.29\t1\t0.255\t1\t1
Primary bile acid biosynthesis\t46\t0.333\t1\t0.288\t1\t1
Purine metabolism\t70\t0.507\t1\t0.406\t1\t1" 

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
  labs( 
    x = "Rich Factor (Hits / Total Genes)", 
    y = "" 
  ) +
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

#### Figure 3  ####
#### A comparison between the PARPi and CTL groups at the late 2-cell stage yielded 25 significantly upregulated metabolites. Subsequently, these metabolites were subjected to KEGG pathway enrichment analysis using the MetaboAnalyst platform (https://www.metaboanalyst.ca).
#### KEGG pathway enrichment analysis for distinct metabolic clusters was performed using the MetaboAnalyst web server (version 6.0.0, https://www.metaboanalyst.ca), which operates via the underlying MetaboAnalystR package (version 4.0)”
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
mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway", 2) 
mSet <- CalculateHyperScore(mSet) 
mSet <- PlotORA(mSet, "ora_0_", "net", "png", 150, width=NA) 
mSet <- PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 150, width=NA) 
mSet <- CalculateHyperScore(mSet) 
mSet <- PlotORA(mSet, "ora_1_", "net", "png", 150, width=NA) 
mSet <- PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 150, width=NA) 

library(ggplot2) 
library(dplyr) 
library(forcats) 

data_text <- "Pathway  total  expected  hits  Raw_p  Holm_p  FDR
Urea Cycle  28  0.419  4  0.000541  0.053  0.0348
Glucose-Alanine Cycle  13  0.195  3  0.000711  0.069  0.0348
Nicotinate and Nicotinamide Metabolism  35  0.524  4  0.0013  0.125  0.0424
Purine Metabolism  73  1.09  5  0.00302  0.287  0.0741
Glutamate Metabolism  48  0.719  4  0.00429  0.403  0.0841
Warburg Effect  57  0.853  4  0.00803  0.747  0.093
Lysine Degradation  30  0.449  3  0.00865  0.796  0.093
Malate-Aspartate Shuttle  10  0.15  2  0.00879  0.8  0.093
Glycine and Serine Metabolism  59  0.883  4  0.00909  0.818  0.093
Ammonia Recycling  31  0.464  3  0.00949  0.845  0.093" 

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
  labs( 
    x = "Rich Factor (Hits / Total Genes)", 
    y = "" 
  ) +
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


#### Figure 4 ####
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
