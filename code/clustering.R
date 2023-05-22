install.packages("Seurat")
install.packages("tidyr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("gridExtra")
install.packages("rjson")
install.packages("harmony")
install.packages("xlsx")

library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(rjson)
library(harmony)
library(xlsx)


set.seed(1234)

#Download patient samples and merge 
path <- "F:/IB/Proect_itiog/data/sampels"
dirs <- list.dirs(path, recursive = F, full.names = F)
for(x in dirs){
  name <-  x
  cts <- Read10X(data.dir = paste0(path, "/", x))
  assign(name, CreateSeuratObject(counts=cts))
}

merge_seurat <- merge(MD01_010_normal_1, y = c(MD01_010_tumor_1,
                               MD043_006_normal_1, MD043_006_tumor_1,
                               MD043_008_normal_1, MD043_008_tumor_1, 
                               NY016_007_normal_1, NY016_007_tumor_1), 
                      add.cell.ids = ls()[c(3:8, 10:11)])

head(merge_seurat)
tail(merge_seurat)

#calculate mitochondial and ribosomal transcript percentages
merge_seurat[["percent.mt"]] <- PercentageFeatureSet(merge_seurat, pattern = "^MT[-\\.]")
merge_seurat[["percent.rb"]] <- PercentageFeatureSet(merge_seurat, pattern = "^RP[SL]")

#create a sample column
merge_seurat$sample <- rownames(merge_seurat@meta.data)

#split sample column
merge_seurat@meta.data <- separate(merge_seurat@meta.data, col= "sample", 
                                   into = c('Patient', 'ID', 'Type', 'num', 'Barcode'),
                                   sep = '_')
merge_seurat@meta.data <- unite(merge_seurat@meta.data, Patient, c(Patient, ID))

meta <- merge_seurat@meta.data
dim(meta)
unique(meta$Patient)
unique(meta$Type)
unique(meta$Barcode)


#Quality control
#Fig.1. Violin plot for the metrics: nFeature, nCount, percent.mt, percent.rb per sample
VlnPlot(merge_seurat, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
        split.by = 'Patient', group.by = 'Patient', log = TRUE, ncol = 4, pt.size=0)
ggsave("F:/IB/Proect_itiog/grafics/itog/vilonplot_per_samples.png")

#Fig.2. Violin plot for the metrics: nFeature, nCount, percent.mt, percent.rb per sample
VlnPlot(merge_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
        log = TRUE, ncol = 4, pt.size=0)
ggsave("F:/IB/Proect_itiog/grafics/itog/vilonplot.png")

#correlation
#Fig.3. Correlation plot nCount and percent.mt
FeatureScatter(merge_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave("F:/IB/Proect_itiog/grafics/itog/nCount_mt.png")

#Fig.4. Correlation plot nCount and  nFeature 
FeatureScatter(merge_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave("F:/IB/Proect_itiog/grafics/itog/nCount_nFeature.png")

FeatureScatter(merge_seurat, feature1 = "percent.rb", feature2 = "percent.mt")
FeatureScatter(merge_seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt")

#Filtering
filtered_seurat <- subset(merge_seurat, subset = nFeature_RNA > 250 & 
                            nFeature_RNA < 2500 & percent.mt < 10)

all.genes <- rownames(filtered_seurat)
length(all.genes) #33538

#Filtered out mitochondrial and ribosomal genes from analysis
var_regex_all = '^MT-|^RP[SL]'
all.genes = grep(var_regex_all, all.genes, invert=T, value=T)
length(all.genes) #33421

filtered_seurat <- subset(filtered_seurat, features = all.genes)

VlnPlot(filtered_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
        log = TRUE, ncol = 4, pt.size=0)
FeatureScatter(filtered_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(filtered_seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt")

#number of cells before filtration
merge_seurat@meta.data %>% group_by(Patient, Type) %>% summarise(n =  n())

#number of cells after filtration
filtered_seurat@meta.data %>% group_by(Patient, Type) %>% summarise(n =  n())


#PCA
filtered_seurat <- NormalizeData(object = filtered_seurat)

filtered_seurat <- FindVariableFeatures(object = filtered_seurat,
                                        selection.method = "vst", nfeatures = 3000)
#Fig.5. Most variable genes 
top10 <- head(VariableFeatures(filtered_seurat), 10)
plot1 <- VariableFeaturePlot(filtered_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,
                     xnudge = 0, ynudge = 0)
plot2
ggsave("F:/IB/Proect_itiog/grafics/itog/variable_genes.png")

#remowe TCR genes, type I Interferon, immunoglobulin genes
IF <- fromJSON(file = 'F:/IB/Proect_itiog/data/ST_TYPE_I_INTERFERON_PATHWAY.v2023.1.Hs.json')
IF$ST_TYPE_I_INTERFERON_PATHWAY$geneSymbols

most_varibl_genes <- VariableFeatures(filtered_seurat)
length(most_varibl_genes) #3000

var_regex = '^TR[ABDG][VJC]|^IG[HJKLM]|^IFNAR1|^IFNB1|^IFNAR1|^IFNB1|^IRF9|^JAK1|^PTPRU|^REG1A|^STAT1|^STAT2|^TYK2'

most_varibl_genes = grep(var_regex, most_varibl_genes, invert=T, value=T)
length(most_varibl_genes) # 2750

filtered_seurat <- ScaleData(object = filtered_seurat, vars.to.regress = 'nCount_RNA')
filtered_seurat <- RunPCA(object = filtered_seurat, features = most_varibl_genes)

ElbowPlot(filtered_seurat)

#Fig.6. PCA by sample type
DimPlot(filtered_seurat, reduction = "pca", group.by = 'Type')
ggsave("F:/IB/Proect_itiog/grafics/itog/PCA_type.png")

#UMAP
filtered_seurat <- RunUMAP(object = filtered_seurat, reduction = "pca", dims = 1:20) 

#Fig.7. Clustering by patient before removing batch effects
DimPlot(filtered_seurat, reduction = "umap", group.by = 'Patient',  repel = T)
ggsave("F:/IB/Proect_itiog/grafics/itog/umap_patient_before_harmony.png")

#harmony
filtered_seurat_harmony <- RunHarmony(filtered_seurat , group.by.vars = "Patient",
                                      dims.use = 1:20)
filtered_seurat_harmony@reductions
filtered_seurat_harmony.emb <- Embeddings(filtered_seurat_harmony, "harmony")
filtered_seurat_harmony.emb[1:10,1:10]

filtered_seurat_harmony <- RunUMAP(filtered_seurat_harmony, reduction = "harmony",
                                   dims = 1:20)
filtered_seurat_harmony <- FindNeighbors(filtered_seurat_harmony, 
                                         reduction = "harmony", dims = 1:20)
filtered_seurat_harmony <- FindClusters(object = filtered_seurat_harmony, 
                                        resolution = 0.4) #13 clusters

table(filtered_seurat_harmony@meta.data$seurat_clusters)

#Fig.8. Clustering by patient after removing batch effects
DimPlot(filtered_seurat_harmony, reduction = "umap", group.by = 'Patient',  repel = T)
ggsave("F:/IB/Proect_itiog/grafics/itog/umap_patient_after_harmony.png")

#subset of marker genes
#CD39 = ENTPD1, HOBIT=ZNF683, CD103=ITGAE, PD1 = PDCD1, TIM3 = HAVCR2 
FeaturePlot(filtered_seurat_harmony, 
            features = c("CD8A", "GZMK","TCF7", "PDCD1","TIGIT",
                         "CD4",  "ITGAE", "CXCL13", "CTLA4", "ENTPD1",
                         "FOXP3","ZNF683","SLC4A10", "HAVCR2", "LAG3"),
            ncol = 5) & NoLegend() & NoAxes() 

#remove 11-13 clusters
sub_obj <- subset(filtered_seurat_harmony, idents = c(0:10))
table(sub_obj@meta.data$seurat_clusters)

DimPlot(sub_obj, reduction = "umap", repel = T, label = T)
ggsave("F:/IB/Proect_itiog/grafics/itog/umap_harmony_11_clusters.png")

#Fig.9. Subset of marker genes
#CD39 = ENTPD1, HOBIT=ZNF683, CD103=ITGAE, PD1 = PDCD1, TIM3 = HAVCR2
FeaturePlot(sub_obj, features = c("CD8A", "GZMK","TCF7", "PDCD1","TIGIT",
                                  "CD4",  "ITGAE", "CXCL13", "CTLA4", "ENTPD1",
                                  "FOXP3","ZNF683","SLC4A10", "HAVCR2", "LAG3"),
            ncol = 5) & NoLegend() & NoAxes()
ggsave("F:/IB/Proect_itiog/grafics/itog/umap_marker_genes.png")


#Differential expression and marker selection
ct_markers <- c("EOMES","GZMK","CRTAM", "NKG7", "GNLY", "S1PR5",
                "LINC02446", "ZNF683","ITGAE", "STMN1", "TUBB", "MKI67",
                "HLA-DRA", "HLA-DQA1", "HLA-DQB1", "CCR6", "IL4I1", 
                "SLC4A10", "SELL", "CCR7", "TCF7", "IL7R", 
                "GPR183", "CD40LG", "PLIN2", "CXCR6", "ALOX5AP",
                "MT1X", "MT1E", "S100A11", "CXCL13", "FAAH2", "NR3C1",
                "MAF", "PTPN13", "KLRB1", "FOXP3", "CCR8", "IL2RA")

DoHeatmap(sub_obj, features = ct_markers,
          cells = sample(1:ncol(filtered_seurat), 5000), angle = 70) + NoLegend()
ggsave("F:/IB/Proect_itiog/grafics/itog/heatmap_11clusters.png")

all.markers <- FindAllMarkers(filtered_seurat_harmony, only.pos = T, min.pct = 0.5,
                              logfc.threshold = 0.5, 
                              max.cells.per.ident = 5000)
write.xlsx(all.markers, 'F:/IB/Proect_itiog/marker_genes_harmony.xlsx')

all.markers_sub <- FindAllMarkers(sub_obj, only.pos = T, min.pct = 0.5,
                                  logfc.threshold = all.markers_sub %>%0.5, 
                                  features = most_varibl_genes,
                                  max.cells.per.ident = 5000)
write.xlsx(all.markers_sub, 'F:/IB/Proect_itiog/marker_genes_sub_harmony.xlsx')

new.cluster.ids <- c("CD8-effector(1)", "MAIT",
                     "CD8-proliferating", "CD4-T(FH)", "CD4-Treg", 
                     "CD4-T(H-1)", "CD8-TRM(1)", "Stem-like memory", 
                     "CD8-effector(2)", "CD8-effector(3)", "CD8-TRM(2)")

names(new.cluster.ids) <- levels(sub_obj)
sub_obj <- RenameIdents(sub_obj, new.cluster.ids)
sub_obj$CellType <- Idents(sub_obj)

#Fig.10. Clustering by type of T-cell 
DimPlot(sub_obj, reduction = "umap", label = TRUE, label.size = 2.5,label.box = TRUE,
        pt.size = 0.5) + NoLegend()
ggsave("F:/IB/Proect_itiog/grafics/itog/umap_11clusters_label.png")

#Fig.11. Heatmap of differentially expressed genes
DoHeatmap(sub_obj, features = ct_markers, label = T,
          cells = sample(1:ncol(sub_obj), 5000), angle = 50, size = 3) + 
  guides(color="none") + theme(text = element_text(size = 10))
ggsave("F:/IB/Proect_itiog/grafics/itog/heatmap_11clusters_label.png")

#Fig.12. Clusters split by type of tissue
ggplot(sub_obj@meta.data, aes(x=factor(seurat_clusters,
                              labels = c("CD8-effector(1)", "MAIT",
                                         "CD8-proliferating", "CD4-T(FH)", "CD4-Treg", 
                                         "CD4-T(H-1)", "CD8-TRM(1)", "Stem-like memory", 
                                         "CD8-effector(2)", "CD8-effector(3)", "CD8-TRM(2)")), 
                              fill=Type)) + 
  geom_bar(position="dodge")+ theme_classic() + xlab("Clusters")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  scale_fill_discrete(name = "Type of tissue")
ggsave("F:/IB/Proect_itiog/grafics/itog/barplot_type_dodge_harmony.png")

#Fig.13. Clusters split by patients
ggplot(data = sub_obj@meta.data, aes(factor(seurat_clusters,
                              labels = c("CD8-effector(1)", "MAIT",
                                         "CD8-proliferating", "CD4-T(FH)", "CD4-Treg", 
                                         "CD4-T(H-1)", "CD8-TRM(1)", "Stem-like memory", 
                                         "CD8-effector(2)", "CD8-effector(3)", "CD8-TRM(2)")),
                                     fill = Patient)) +
  geom_bar(position="dodge") + theme_classic() + xlab("Clusters")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))+
  scale_fill_discrete(name = "Patient")
ggsave("F:/IB/Proect_itiog/grafics/itog/barplot_patient_dodge_harmony.png")

#Fig.14. Umap split by type of tissue
DimPlot(sub_obj, reduction = "umap", label = TRUE, label.size = 2.5,label.box = TRUE,
        pt.size = 0.5, split.by = "Type") + NoLegend()
ggsave("F:/IB/Proect_itiog/grafics/itog/umap_11clusters_typy_of_tissue.png")

#Save seurat object
saveRDS(filtered_seurat, file='F:/IB/Proect_itiog/processed_seurat_object_4s.rds')
saveRDS(filtered_seurat_harmony, file='F:/IB/Proect_itiog/processed_seurat_object_4s_harmony.rds')
saveRDS(sub_obj, file='F:/IB/Proect_itiog/processed_seurat_object_sub_4s_harmony.rds')



