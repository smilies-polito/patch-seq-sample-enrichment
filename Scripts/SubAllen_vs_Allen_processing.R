# <patch-seq-sample-enrichment>
#
# Copyright © 2022 Politecnico di Torino, Control and Computer Engineering Department, SMILIES group
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


########################################################################################
#
# INSTALL REQUIRED PACKAGES IF NEEDED
# 
# N.B. This may require the installation of local libraries. Please check the README
# file of the project for a list of required packages
#
# ######################################################################################
if (!requireNamespace("httr", quietly = TRUE)) 
  install.packages("httr", dependencies = c("Depends"))

if (!requireNamespace("igraph", quietly = TRUE)) 
  install.packages("igraph", dependencies = c("Depends"))

if (!requireNamespace("systemfonts", quietly = TRUE)) 
  install.packages("systemfonts", dependencies = c("Depends"))

if (!requireNamespace("qpdf", quietly = TRUE)) 
  install.packages("qpdf", dependencies = c("Depends"))

if (!requireNamespace("textshaping", quietly = TRUE)) 
  install.packages("textshaping", dependencies = c("Depends"))

if (!requireNamespace("ragg", quietly = TRUE)) 
  install.packages("ragg", dependencies = c("Depends"))

if (!requireNamespace("pkgdown", quietly = TRUE)) 
  install.packages("pkgdown", dependencies = c("Depends"))

if (!requireNamespace("devtools", quietly = TRUE)) 
  install.packages("devtools", dependencies = c("Depends"))

if (!requireNamespace("leiden", quietly = TRUE)) 
  install.packages("leiden", dependencies = c("Depends"))

if (!requireNamespace("plotly", quietly = TRUE)) 
  install.packages("plotly", dependencies = c("Depends"))

install.packages('BiocManager', dependencies = c("Depends"))
BiocManager::install()
BiocManager::install('multtest')
install.packages('Seurat')

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("Seurat", dependencies = TRUE)

if (!requireNamespace("Seurat", quietly = TRUE)) 
  install.packages("R.utils", dependencies = c("Depends"))

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) 
  remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("dplyr", quietly = TRUE)) 
  install.packages("dplyr", dependencies = c("Depends"))

if (!requireNamespace("aricode", quietly = TRUE)) 
  install.packages("aricode", dependencies = c("Depends"))

if (!requireNamespace("data.table", quietly = TRUE)) 
  install.packages("data.table", dependencies = c("Depends"))

if (!requireNamespace("cellranger", quietly = TRUE)) 
  install.packages("cellranger", dependencies = c("Depends"))

if (!requireNamespace("readxl", quietly = TRUE)) 
  install.packages("readxl", dependencies = c("Depends"))

if (!requireNamespace("sctransform", quietly = TRUE)) 
  install.packages("sctransform", dependencies = c("Depends"))

if (!requireNamespace("class", quietly = TRUE)) 
  install.packages("class", dependencies = c("Depends"))

########################################################################################
#
# LOAD REQUIRED PACKAGES
# 
#######################################################################################

library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(aricode)
library(data.table)
library(readxl)
library(sctransform)
library(class)
library(ggplot2)


########################################################################################
#
# RI functions definition
# 
#######################################################################################

cell_similarity <- function (cell1, cell2)
{
  cell1 <- as.numeric(cell1)
  cell1[(cell1) > 0] <- 1
  
  cell2 <- as.numeric(cell2)
  cell2[(cell2) > 0] <- 1
  RI(cell1,cell2)
}

cell_similarity_hard <- function (cell2, cell1)
{
  
  cell1$gene_names <- rownames(cell1)
  colnames(cell1)[1] <- "expr"
  cell1$expr[(cell1$expr) > 0] <- 1
  cell1 <- cell1[cell1$expr > 0,]
  
  cell2$gene_names <- rownames(cell2)
  colnames(cell2)[1] <- "expr"
  cell2 <- cell2[cell2$gene_names %in% cell1$gene_names,]
  cell2$expr[(cell2$expr) > 0] <- 1
  
  RI(cell1$expr,cell2$expr)
}

cell_similarity_level <- function (cell1, cell2, level1, level2)
{
  cell1 <- as.numeric(cell1)
  cell1[cell1 < (level1$mean/4)] <- 0
  cell1[(cell1) > 0] <- 1
  
  cell2 <- as.numeric(cell2)
  cell2[cell2 < (level2$mean/4)] <- 0
  cell2[(cell2) > 0] <- 1
  RI(cell1,cell2)
  
  
}

cell_similarity_level_hard <- function (cell2, cell1, level2, level1)
{
  
  cell1$gene_names <- rownames(cell1)
  colnames(cell1)[1] <- "expr"
  cell1 <- cell1[cell1$expr > 0,]
  level1 <- level1[level1$gene_name %in% cell1$gene_names,]
  cell1$expr[(cell1$expr) < (level1$mean/4)] <- 0
  cell1$expr[(cell1$expr) > 0] <- 1
  
  
  cell2$gene_names <- rownames(cell2)
  colnames(cell2)[1] <- "expr"
  cell2 <- cell2[cell2$gene_names %in% cell1$gene_names,]
  level2 <- level2[level2$gene_name %in% cell1$gene_names,]
  cell2$expr[(cell2$expr) < (level2$mean/4)] <- 0
  cell2$expr[(cell2$expr) > 0] <- 1
  
  RI(cell1$expr,cell2$expr)
}

cell_similarity_variable <- function (cell2, cell1, features)
{
  
  cell1$gene_names <- rownames(cell1)
  colnames(cell1)[1] <- "expr"
  cell1 <- cell1[cell1$gene_names %in% features,]
  cell1$expr[cell1$expr > 0] <- 1
  
  
  #cell2 <- as.data.frame(ref_mat[,103])
  cell2$gene_names <- rownames(cell2)
  colnames(cell2)[1] <- "expr"
  cell2 <- cell2[cell2$gene_names %in% features,]
  cell2$expr[cell2$expr > 0] <- 1
  RI(cell1$expr,cell2$expr)
}

cell_similarity_variable_hard <- function (cell2, cell1, features)
{
  
  cell1$gene_names <- rownames(cell1)
  colnames(cell1)[1] <- "expr"
  cell1 <- cell1[cell1$gene_names %in% features,]
  cell1$expr[cell1$expr > 0] <- 1
  cell1 <- cell1[cell1$expr > 0,]
  
  #cell2 <- as.data.frame(ref_mat[,103])
  cell2$gene_names <- rownames(cell2)
  colnames(cell2)[1] <- "expr"
  cell2 <- cell2[cell2$gene_names %in% cell1$gene_names,]
  cell2$expr[cell2$expr > 0] <- 1
  RI(cell1$expr,cell2$expr)
}

############ loading processing of Allen #############

allen_reference <- readRDS("DATA/allen_brain.rds")
Idents(allen_reference) <- "class"

#choosing only neurons cell types#
allen_reference <- subset(allen_reference, idents = c("GABAergic", "Glutamatergic"))
sub_ref <- allen_reference
rm(allen_reference)

#DimPlot(sub_ref, pt.size = 1, group.by = "subclass") +labs(title = "Allen cells") +theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

########### creation of the subset of Allen as query and its processing ########

memory.size(30000)
reference_cells <- Cells(sub_ref)
subref_cells <- sample(reference_cells, 80)
DefaultAssay(sub_ref) <- "RNA"
sub_ref_query <- sub_ref[,subref_cells]
sub_ref <- sub_ref[,!(Cells(sub_ref) %in% subref_cells)]
sub_ref <- SCTransform(sub_ref)
DefaultAssay(sub_ref) <- "SCT"

sub_ref_query <- CreateSeuratObject(GetAssayData(object = sub_ref_query, slot = "counts"))
sub_ref_query <- FindVariableFeatures(sub_ref_query, nfeatures = 3000)
sub_ref_query <- SCTransform(sub_ref_query)
DefaultAssay(sub_ref_query) <- "SCT"
sub_ref_query <- RunPCA(sub_ref_query)
sub_ref_query <- RunUMAP(sub_ref_query, dims = 1:30, reduction.name = "umap")
sub_ref_query <- FindNeighbors(sub_ref_query, dims = 1:30)
sub_ref_query <- FindClusters(sub_ref_query, resolution = 0.01, algorithm = 3)

pdf("TMP/AA_sub_allen_cells.pdf")
DimPlot(sub_ref_query)
dev.off()

saveRDS(sub_ref,"TMP/AR_SR_REF_obj")
saveRDS(sub_ref_query,"TMP/AR_SR_query_obj")


####### Finding Anchors ########

DefaultAssay(sub_ref_query) <- "SCT"
DefaultAssay(sub_ref) <- "SCT"
features_SCT <- SelectIntegrationFeatures(object.list = c(sub_ref, sub_ref_query))
obj.list <- PrepSCTIntegration(object.list = c(sub_ref, sub_ref_query), anchor.features = features_SCT)

anchors_SCT <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT",
                                      anchor.features = features_SCT)



############# data integration ##########

integrated <- IntegrateData(anchorset = anchors_SCT, normalization.method = "SCT", dims = 1:30, k.weight = 80)
integrated <- ScaleData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "umap")
integrated <- FindNeighbors(integrated, dims = 1:30,return.neighbor = TRUE)

pdf("TMP/AA_integrated.pdf")
DimPlot(integrated, pt.size = 2, group.by = "subclass") +labs(title = "Integrated visualization") +theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + scale_shape_manual(values = c(1, 4))
dev.off()

saveRDS(integrated,"TMP/AR_SR_INT_obj")


######### anchors metrics computation #######

anchors_couples_SCT <- as.data.frame(anchors_SCT@anchors)
setorderv(anchors_couples_SCT, cols = "score", order = -1)

anchors_couples_SCT <- anchors_couples_SCT[anchors_couples_SCT$dataset1 == 1,]
anchors_couples_SCT$cell1_name <- Cells(sub_ref)[anchors_couples_SCT$cell1]
anchors_couples_SCT$cell2_name <- Cells(sub_ref_query)[anchors_couples_SCT$cell2]


features <- row.names(sub_ref)
features <- features[features %in% row.names(sub_ref_query)]
ref_mat <- sub_ref@assays[["SCT"]]@counts
ref_mat <- ref_mat[features,]

query_mat <- sub_ref_query@assays[["SCT"]]@counts
query_mat <- query_mat[features,] 

#loading list of electrophysiology genes #
elettrophys.genes <- read.table("DATA/elettrophys_genes.txt", quote="\"", comment.char="")


#creation of table with all the anchors, and  computation #
anchors_couples_SCT$label <- sub_ref@meta.data[Cells(sub_ref)[anchors_couples_SCT$cell1],][["subclass"]]
anchors_couples_SCT$label_query <- sub_ref@meta.data[anchors_couples_SCT$cell2_name,][["subclass"]]

anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI = cell_similarity(as.character(ref_mat[,cell1]),as.character(query_mat[,cell2])))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_var = cell_similarity_variable(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = features_SCT))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_eletp = cell_similarity_variable(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = elettrophys.genes$V1))

anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_strict = cell_similarity_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2])))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_var_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = features_SCT))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_eletp_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = elettrophys.genes$V1))

int_mat <- integrated@assays[["integrated"]]@data

anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_int = cell_similarity(as.character(int_mat[,cell1_name]),as.character(int_mat[,cell2_name])))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise() %>% mutate(RI_eletp_int = cell_similarity_variable(as.data.frame(int_mat[,cell1_name]),as.data.frame(int_mat[,cell2_name]), features = elettrophys.genes$V1))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise() %>% mutate(RI_var_int = cell_similarity_variable(as.data.frame(int_mat[,cell1_name]),as.data.frame(int_mat[,cell2_name]), features = features_SCT))

anchors_couples_SCT$RI_mean <- rowMeans(anchors_couples_SCT[,c("RI", "RI_var", "RI_eletp","RI_strict", "RI_var_strict", "RI_eletp_strict","RI_int", "RI_var_int", "RI_eletp_int")])
anchors_couples_SCT$RI_mean_loose <- rowMeans(anchors_couples_SCT[,c("RI", "RI_var", "RI_eletp")])
anchors_couples_SCT$RI_mean_strict <- rowMeans(anchors_couples_SCT[,c("RI_strict", "RI_var_strict", "RI_eletp_strict")])
anchors_couples_SCT$RI_mean_int <- rowMeans(anchors_couples_SCT[,c("RI_int", "RI_var_int", "RI_eletp_int")])
anchors_couples_SCT$RI_mean_var <- rowMeans(anchors_couples_SCT[,c("RI_var_strict", "RI_var_int")])
anchors_couples_SCT$RI_mean_eletp <- rowMeans(anchors_couples_SCT[,c("RI_eletp_strict", "RI_eletp_int")])


saveRDS(anchors_couples_SCT,"TMP/AR_SR_anchors_couples_SCT")

anchors_couples_top <- anchors_couples_SCT %>% group_by(cell2) %>% top_n(1, RI_mean_strict)

pdf("TMP/anchors_refrence_cells.pdf")
DimPlot(integrated,group.by = "subclass", cells = c(anchors_couples_top$cell1_name, anchors_couples_top$cell2_name) )
dev.off()

######### Nearest neighbor ########

train <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
train <- train[1:length(Cells(sub_ref)),]
train_label <- rownames(train)
test_data <- as.data.frame(integrated@reductions[["umap"]]@cell.embeddings)
test_data <- test_data[(length(Cells(sub_ref))+1):length(rownames(test_data)),]


NN_couples <- test_data
NN_couples$query_cells <- rownames(NN_couples)
NN_couples$ref_cells <- as.character(knn(train, test_data, train_label, k=1))
NN_couples$UMAP_1 <- NULL
NN_couples$UMAP_2 <- NULL

pdf("TMP/NN_refrence_cells.pdf")
DimPlot(integrated,group.by = "subclass", cells = c(NN_couples$query_cells, NN_couples$ref_cells))
dev.off()

ìNN_couples$label <- integrated@meta.data[NN_couples$ref_cells,][["subclass"]]

NN_couples <- NN_couples %>% rowwise()%>% mutate(RI = cell_similarity(as.character(ref_mat[,ref_cells]),as.character(query_mat[,query_cells])))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_var = cell_similarity_variable(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = features_SCT))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_eletp = cell_similarity_variable(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = elettrophys.genes$V1))

NN_couples <- NN_couples %>% rowwise()%>% mutate(RI_strict = cell_similarity_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells])))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_var_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = features_SCT))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_eletp_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = elettrophys.genes$V1))

NN_couples <- NN_couples %>% rowwise()%>% mutate(RI_int = cell_similarity(as.character(int_mat[,ref_cells]),as.character(int_mat[,query_cells])))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_var_int = cell_similarity_variable(as.data.frame(int_mat[,ref_cells]),as.data.frame(int_mat[,query_cells]), features = features_SCT))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_eletp_int = cell_similarity_variable(as.data.frame(int_mat[,ref_cells]),as.data.frame(int_mat[,query_cells]), features = elettrophys.genes$V1))


NN_couples$RI_mean <- rowMeans(NN_couples[,c("RI", "RI_var", "RI_eletp","RI_strict", "RI_var_strict", "RI_eletp_strict","RI_int", "RI_var_int", "RI_eletp_int")])
NN_couples$RI_mean_loose <- rowMeans(NN_couples[,c("RI", "RI_var", "RI_eletp")])
NN_couples$RI_mean_strict <- rowMeans(NN_couples[,c("RI_strict", "RI_var_strict", "RI_eletp_strict")])
NN_couples$RI_mean_int <- rowMeans(NN_couples[,c("RI_int", "RI_var_int", "RI_eletp_int")])
NN_couples$RI_mean_var <- rowMeans(NN_couples[,c("RI_var_strict", "RI_var_int")])
NN_couples$RI_mean_eletp <- rowMeans(NN_couples[,c("RI_eletp_strict", "RI_eletp_int")])



####### plotting of anchors and NN #######

ac <- anchors_couples_top[,c("RI_mean","RI_mean_loose","RI_mean_strict","RI_mean_int","RI_var_strict","RI_var_int","RI_mean_var","RI_eletp_strict","RI_eletp_int","RI_mean_eletp" )]
ac$type <- "Anchors"
nc <- NN_couples[,c("RI_mean","RI_mean_loose","RI_mean_strict","RI_mean_int","RI_var_strict","RI_var_int","RI_mean_var","RI_eletp_strict","RI_eletp_int","RI_mean_eletp" )]
nc$type <- "NN couples"

saveRDS(NN_couples,"TMP/AR_SR_NN_couples")
saveRDS(anchors_couples_SCT,"TMP/AR_SR_anchors_couples_SCT")

results <- rbind(ac , nc)

pdf("TMP/AA_RI_mean_RI_mean_s.pdf")
ggplot(results, aes (RI_mean_strict, RI_mean, color = type)) + geom_point(size = 4) + labs(x = expression("RI"[mean]^s), y = expression("RI"[mean]), color = "Couples", title = "SubAllen vs Allen")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.title=element_text(size=20), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text = element_text(size=15))
dev.off()

pdf("TMP/AA_RI_mean_RI_mean_int.pdf")
ggplot(results, aes (RI_mean_int, RI_mean, color = type)) + geom_point(size = 4) + labs(x = expression("RI"[mean]^int), y = expression("RI"[mean]), color = "Couples", title = "SubAllen vs Allen")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.title=element_text(size=20), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text = element_text(size=15))
dev.off()

pdf("TMP/AA_RI_var_RI_eph.pdf")
ggplot(results, aes (RI_mean_var, RI_mean_eletp, color = type)) + geom_point(size = 4) + labs(x = expression("RI"[var_mean]), y = expression("RI"[eph_mean]), color = "Couples", title = "SubAllen vs Allen")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.title=element_text(size=20), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text = element_text(size=15))
dev.off()

mean(ac$RI_mean)
mean(nc$RI_mean)

mean(ac$RI_mean_strict)
mean(nc$RI_mean_strict)
