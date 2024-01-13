setwd("/media/heshidian/95B8424A542FB26C/Bioinfor_xianyu/caizhiwei")
library(data.table)

# scripts <- list.files(path = "./zellkonverter-master/R", recursive = F, full.names = F)
# for (x in scripts) {
#   source(paste0("./zellkonverter-master/R/",x))
# }

### 
library(Seurat)
library(reticulate)
reticulate::py_config()
ad <- import("anndata", convert = FALSE)
h5ad <- "GSE202051_totaldata-final-toshare.h5ad"
h5ad <- ad$read_h5ad(h5ad)
h5ad$obs$to_csv("col_data.csv")
rm(h5ad)
data_matrix <- fread("matrix.tsv", sep = "\t")
data_matrix <- as.data.frame(data_matrix)
rownames(data_matrix) <- data_matrix[,1]
data_matrix <- data_matrix[,-1]
sce <- CreateSeuratObject(data_matrix)
rm(data_matrix)
col_data <- read.csv("col_data.csv", row.names = 1)
sce <- AddMetaData(sce, col_data)

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.1) 
sce <- RunUMAP(sce, dims = 1:20)
sce <- RunTSNE(sce, dims = 1:20)

cell_types <- unique(sce@meta.data[["Level.1.Annotation"]])
Idents(sce) <- sce$Level.1.Annotation
for (i in cell_types) {
  print(i)
  temp <- subset(sce, idents = i)
  saveRDS(temp, paste0(i,".rds"))
}

saveRDS(sce, "all_seurat.rds")


###
adult_pancreas <- readRDS("Documents/zhuyi_bio/pdac-single-cell/adult_pancreas_2020.rds")
Idents(adult_pancreas) <- adult_pancreas$age
adult_pancreas <- subset(adult_pancreas, idents = c("30","46","53","59","77"))
Idents(adult_pancreas) <- adult_pancreas$Cluster
normal_Schwann <- subset(adult_pancreas, idents = "Schwann")
normal_Schwann$type <- "normal_Schwann"


sce <- readRDS("Documents/zhuyi_bio/pdac-single-cell/all_seurat.rds")
Idents(sce) <- sce$response
sce <- subset(sce, idents = "Untreated")
Idents(sce) <- sce$Level.1.Annotation
cancer_Schwann <- subset(sce, idents = "Schwann")
cancer_Schwann$type <- "cancer_Schwann"


Schwann <- merge(x = cancer_Schwann, y = list(normal_Schwann))
Idents(Schwann) <- Schwann$type




diff_markers <- FindMarkers(Schwann, ident.1 = "cancer_Schwann", ident.2 = "normal_Schwann", 
                            logfc.threshold = 0)
# volcano
DOWN_row <- intersect(which(diff_markers$p_val_adj < 0.05), which(diff_markers$avg_log2FC < -1))
UP_row <- intersect(which(diff_markers$p_val_adj < 0.05), which(diff_markers$avg_log2FC > 1))
DOWN <- diff_markers[DOWN_row,]
UP <- diff_markers[UP_row,]
write.table(DOWN, "DOWN.csv", col.names = T, sep = ',', quote = FALSE, row.names = T)
write.table(UP, "UP.csv", col.names = T, sep = ',', quote = FALSE, row.names = T)
DOWN <- DOWN[order(DOWN$avg_log2FC,decreasing = F),]
UP <- UP[order(UP$avg_log2FC,decreasing = T),]
sigdiff <- as.data.frame(rbind(UP[1:5,],
                               DOWN[1:5,]))
# sigdiff <- as.data.frame(rbind(UP,
#                                DOWN))
library(ggrepel)
library(ggplot2)
color <- rep("gray", length(diff_markers$avg_log2FC))
cex <- rep(1.2, length(diff_markers$avg_log2FC))
cex[DOWN_row] <- 2
cex[UP_row] <- 2
cex <- as.character(cex)
color[DOWN_row] <- "blue"
color[UP_row] <- "red"
pdf("volcano.pdf")
vio_plot <- ggplot(data = diff_markers, aes(x = avg_log2FC, y = -log10(p_val_adj), size = cex, colour = color, alpha = "0.3")) +
  theme_bw() +
  labs(x = "log2 FoldChange", y = "log10 adjusted P value") +
  geom_point() +
  geom_hline(aes(yintercept=-log10(0.05)), colour = "gray", linetype="dashed") +
  geom_vline(aes(xintercept= -1), colour="gray", linetype="dashed") +
  geom_vline(aes(xintercept= 1), colour="gray", linetype="dashed") +
  scale_color_manual(name = "",
                     values = c('red' = 'red', "blue" = 'blue', "gray" = "gray"),
                     labels = c("UP", 'DOWN', "Not significant")) +
  scale_size_manual(values = c('1.2' = 1.2, "2" = 2)) +
  scale_alpha_manual(values = c("0.3" = 0.3, "1" = 1)) +
  guides(size = "none", alpha = "none")
vio_plot +
  geom_text_repel(data = sigdiff, aes(x = avg_log2FC, 
                                      y = -log10(p_val_adj), 
                                      label = rownames(sigdiff),
                                      alpha = rep("1",nrow(sigdiff))),
                  size = 3,
                  box.padding = unit(0.8, "lines"),
                  point.padding = unit(0, "lines"), 
                  min.segment.length = 0,
                  segment.color = "black",
                  colour="#000000",
                  show.legend = FALSE,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100)) 
dev.off()

sigdiff2 <- as.data.frame(rbind(UP,
                               DOWN))
Schwann <- NormalizeData(Schwann)
Schwann <- FindVariableFeatures(Schwann)
Schwann <- ScaleData(Schwann)
Schwann <- RunPCA(Schwann)
Schwann <- RunUMAP(Schwann, dims = 1:15)
Schwann <- RunTSNE(Schwann, dims = 1:15)
pdf("Schwann_type_tsne.pdf", width = 8)
DimPlot(Schwann, reduction = "tsne", pt.size = 1.2)
dev.off()

pdf("diff_heatmap.pdf")
DoHeatmap(Schwann, features = rownames(sigdiff2), size = 5) + NoLegend()
dev.off()

# monocle
library(monocle)
Schwann <- Schwann[rowSums(Schwann@assays$RNA@counts)!=0,]
Schwann <- Schwann[,colSums(Schwann@assays$RNA@counts)!=0]
data <- as(as.matrix(Schwann@assays$RNA@counts), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = Schwann@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(cellData = data,
                              phenoData = pd,
                              featureData = fd)
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
library(devtools)
load_all(path = "/home/heshidian/software/monocle_2.24.0/monocle/")
disp_table <- dispersionTable(monocle_cds)
clustering_genes <- subset(disp_table, mean_expression >= 0.5)
monocle_cds <- setOrderingFilter(monocle_cds, clustering_genes$gene_id)
monocle_cds <- reduceDimension(monocle_cds, max_components = 2, method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
plot_cell_trajectory(monocle_cds, color_by = "State")
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
pdf("monocle_group.pdf")
plot_cell_trajectory(monocle_cds, color_by = "type")
dev.off()
pdf("Pseudotime.pdf")
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
dev.off()
saveRDS(monocle_cds, "monocle_cds.rds")
pdf("monocle_heatmap.pdf", width = 6, height = 8)
heatmap_plot <- plot_pseudotime_heatmap(monocle_cds[clustering_genes$gene_id,],
                        num_clusters = 3,
                        cores = 6,
                        show_rownames = F,
                        return_heatmap=T)
dev.off()
clusters <- cutree(heatmap_plot$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
write.csv(clustering, "cluster_genes.csv", quote = F)
