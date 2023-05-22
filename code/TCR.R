install.packages("data.table")
BiocManager::install("scRepertoire")

library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(scRepertoire)

set.seed(1234)

# load seurat object
filtered_seurat <- readRDS(file='F:/IB/Proect_itiog/processed_seurat_object_4s_harmony.rds')
sub_obj <- readRDS(file='F:/IB/Proect_itiog/processed_seurat_object_sub_4s_harmony.rds')


#TCR load and trim
barcoder <- function(df, prefix){
  df$barcode <- paste0(prefix, "_", df$barcode)
  df
}

dirs_tcr <- list.dirs("F:/IB/Proect_itiog/data/tcr", recursive = F, full.names = F)
tcr <- data.frame() 
for(x in dirs_tcr[c(1:8)]){
  tcr_new <- read.csv(paste0("F:/IB/Proect_itiog/data/tcr/", x, "/filtered_contig_annotations.csv"))
  x <- gsub("-", "_", x)
  tcr_new <- barcoder(tcr_new, prefix = x)
  tcr <- rbind(tcr, tcr_new)
  rm(tcr_new)
}

head(tcr$barcode)
tail(tcr$barcode)
head(sub_obj)
tail(sub_obj)
unique(tcr$chain)

names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

data_concater <- function(x){
  x<- levels(factor(x))
  paste(x, collapse = "__")
}

tcr <- as.data.table(tcr)
uniqueN(tcr$clonotype_id) #5025

grpn = uniqueN(tcr$barcode)
pb <- txtProgressBar(min = 0, max = grpn, style = 3)
tcr_collapsed <- tcr[, {setTxtProgressBar(pb,.GRP); lapply(.SD, data_concater)} , by=barcode]

rownames(tcr_collapsed) <- tcr_collapsed$barcode
unique(tcr_collapsed$chain)
uniqueN(tcr_collapsed$clonotype_id) #5025

combined <- AddMetaData(object = sub_obj, metadata = tcr_collapsed)
head(combined@meta.data)
unique(combined@meta.data$chain)
uniqueN(combined@meta.data$clonotype_id)  #4903

#save
saveRDS(combined , file='F:/IB/Proect_itiog/processed_seurat_object_tcr_harmony.rds')

#Fig.15. TRC chain in normal and tumor sample 
ggplot(data = combined@meta.data, aes(factor(Patient), fill = chain))+ 
  geom_bar(position="fill")+ theme_bw()+
  facet_grid(~ Type) + 
  scale_x_discrete(name ="")+
  theme(axis.text.x = element_text( size = 10, hjust = 1, angle = 30))
ggsave("F:/IB/Proect_itiog/grafics/itog/tcr_chain.png")

#Fig.16. Cell with TCR
DimPlot(combined, reduction = 'umap', group.by = 'is_cell')+
  ggtitle("")+
  scale_color_discrete(labels = c("cell with TCR", "cell without TCR"))
ggsave("F:/IB/Proect_itiog/grafics/itog/tcr_cell_is.png")


#Clonotypes
for(x in dirs_tcr[c(1:8)]){
  name <- gsub("-", "_", x)
  tcr_new <- read.csv(paste0("F:/IB/Proect_itiog/data/tcr/", x, "/filtered_contig_annotations.csv"))
  tcr_new$is_cell <- toupper(tcr_new$is_cell)
  tcr_new$productive <- toupper(tcr_new$productive)
  assign(name, tcr_new)
  rm(tcr_new)
}

contig_list <- list(MD01_010_normal_1, MD01_010_tumor_1, 
                    MD043_006_normal_1, MD043_006_tumor_1,
                    MD043_008_normal_1, MD043_008_tumor_1,
                    NY016_007_normal_1, NY016_007_tumor_1)
head(contig_list[[1]])

sample_pat <- c("MD01_010_normal_1", "MD01-010_tumor_1",
             "MD043_006_normal_1", "MD043-006_tumor_1", 
             "MD043_008_normal_1", "MD043-008_tumor_1",
             "NY016_007_normal_1", "NY016-007_tumor_1")

combineT <- combineTCR(contig_list, 
                       samples = sample_pat,
                       cells ="T-AB")

#Persent unique clonotype
quantContig(combineT, cloneCall="gene+nt", scale = FALSE)+
  theme(axis.text.x = element_text( size = 10, hjust = 1, angle = 60))

#lemth of contigs
lengthContig(combineT, cloneCall="aa", chain = "both") 

#Homeostasis
clonalHomeostasis(combineT, cloneCall = "gene", 
                  cloneTypes = c(Rare = 1e-04, 
                                 Small = 0.001, 
                                 Medium = 0.01, 
                                 Large = 0.1, 
                                 Hyperexpanded = 1))+
  theme(axis.text.x = element_text( size = 10, hjust = 1, angle = 50))

#Fig. 17.  Proportion of top N clones per sample 
clonalProportion(combineT, cloneCall = "gene", split = c(10, 100, 1000, 10000))+
  theme(axis.text.x = element_text( size = 10, hjust = 1, angle = 50))
ggsave("F:/IB/Proect_itiog/grafics/itog/clonal_proportion.png")

#make a single list
sub_obj <- combineExpression(combineT, sub_obj, 
                             cloneCall="gene", 
                             group.by = "sample", 
                             proportion = FALSE, 
                             cloneTypes=c(Single=1, Small=5,
                                          Medium=20, Large=100, Hyperexpanded=500))

#Fig. 18. Distribution of the clonotypes
occupiedscRepertoire(sub_obj, x.axis = "ident", label = FALSE)+
  theme(axis.text.x = element_text( size = 10, hjust = 1, angle = 50))
ggsave("F:/IB/Proect_itiog/grafics/itog/clonal_cluster.png")

#Compare Clonotypes
compareClonotypes(combineT, numbers = 5, cloneCall = "aa", graph = "alluvial", 
                  samples = sample_pat, split.by = "Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 5))

#Visualize Gene Usage
vizGenes(combineT, gene = "V", 
         chain = "TRB", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)

#sessionInfo
clustering_tcr_sessionInfo.Rmd <- sessionInfo()
home_dir <- "F:/IB/Proect_itiog"
writeLines(capture.output(sessionInfo()), paste0(home_dir, "/session_info.", format(Sys.time(), "%Y%m%d.%H%M"), ".txt"))
