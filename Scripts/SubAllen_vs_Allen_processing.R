library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(aricode)
library(data.table)
library(readxl)
library(sctransform)
library(class)
library(ggplot2)

########### RI functions definition##########
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

allen_reference <- readRDS("C:/.../allen_brain.rds")
allen_reference <- as.Seurat(allen_reference)
Idents(allen_reference) <- "class"

#choosing onlyneurons cell types#
allen_reference <- subset(allen_reference, idents = c("GABAergic", "Glutamatergic"))
sub_ref <- allen_reference
rm(allen_reference)

DimPlot(sub_ref, pt.size = 1, group.by = "subclass") +labs(title = "Allen cells") +theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 

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

DimPlot(sub_ref_query)

saveRDS(sub_ref,"C:/.../Results_obj/AR_SR_REF_obj")
saveRDS(sub_ref_query,"C:/.../Results_obj/AR_SR_query_obj")


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

DimPlot(integrated, pt.size = 2, group.by = "subclass") +labs(title = "Integrated visualization") +theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) + scale_shape_manual(values = c(1, 4))

DimPlot(integrated, cells = c(anchors_couples_SCT$cell1_name, anchors_couples_SCT$cell2_name), group.by = "subclass", label = TRUE)
saveRDS(integrated,"C:/.../Results_obj/AR_SR_INT_obj")


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

#genes_mean <- as.data.frame(rownames(ref_mat))
#colnames(genes_mean) <- "gene_name"
#genes_mean <- genes_mean %>% rowwise() %>% mutate(mean = sum(ref_mat[gene_name,])/sum(ref_mat[gene_name,]>0))
#saveRDS(genes_mean,"C:/Users/loren/BIBM/2022/gene_mean")
#genes_mean <- readRDS("C:/Users/loren/BIBM/2022/gene_mean")
#genes_mean_query <- as.data.frame(rownames(query_mat))
#colnames(genes_mean_query) <- "gene_name"
#genes_mean_query <- genes_mean_query %>% rowwise() %>% mutate(mean = sum(query_mat[gene_name,])/sum(query_mat[gene_name,]>0))
#genes_mean <- genes_mean[genes_mean$gene_name %in% genes_mean_query$gene_name, ]
#genes_mean_query <- genes_mean_query[genes_mean_query$gene_name %in% genes_mean$gene_name, ]


#loading list of electrophysiology genes #
elettrophys.genes <- read.table("C:/.../elettrophys genes.txt", quote="\"", comment.char="")


#creation of table with all the anchors, and  computation #
anchors_couples_SCT$label <- sub_ref@meta.data[Cells(sub_ref)[anchors_couples_SCT$cell1],][["subclass"]]
anchors_couples_SCT$label_query <- sub_ref@meta.data[anchors_couples_SCT$cell2_name,][["subclass"]]

anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI = cell_similarity(as.character(ref_mat[,cell1]),as.character(query_mat[,cell2])))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_var = cell_similarity_variable(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = var_features))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_eletp = cell_similarity_variable(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = elettrophys.genes$V1))
#anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_level_4 = cell_similarity_level(as.character(ref_mat[,cell1]),as.character(query_mat[,cell2]),genes_mean, genes_mean_query))

anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_strict = cell_similarity_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2])))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_var_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = var_features))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_eletp_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = elettrophys.genes$V1))
#anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_strict_level_4 = cell_similarity_level_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]),genes_mean, genes_mean_query))

#anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_var_i = cell_similarity_variable(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = var_features))
#anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_var_i_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,cell1]),as.data.frame(query_mat[,cell2]), features = var_features))

int_mat <- integrated@assays[["integrated"]]@data

anchors_couples_SCT <- anchors_couples_SCT %>% rowwise()%>% mutate(RI_int = cell_similarity(as.character(int_mat[,cell1_name]),as.character(int_mat[,cell2_name])))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise() %>% mutate(RI_eletp_int = cell_similarity_variable(as.data.frame(int_mat[,cell1_name]),as.data.frame(int_mat[,cell2_name]), features = elettrophys.genes$V1))
anchors_couples_SCT <- anchors_couples_SCT %>% rowwise() %>% mutate(RI_var_int = cell_similarity_variable(as.data.frame(int_mat[,cell1_name]),as.data.frame(int_mat[,cell2_name]), features = var_features))

anchors_couples_SCT$RI_mean <- rowMeans(anchors_couples_SCT[,c("RI", "RI_var", "RI_eletp","RI_strict", "RI_var_strict", "RI_eletp_strict","RI_int", "RI_var_int", "RI_eletp_int")])
anchors_couples_SCT$RI_mean_loose <- rowMeans(anchors_couples_SCT[,c("RI", "RI_var", "RI_eletp")])
anchors_couples_SCT$RI_mean_strict <- rowMeans(anchors_couples_SCT[,c("RI_strict", "RI_var_strict", "RI_eletp_strict")])
anchors_couples_SCT$RI_mean_int <- rowMeans(anchors_couples_SCT[,c("RI_int", "RI_var_int", "RI_eletp_int")])
anchors_couples_SCT$RI_mean_var <- rowMeans(anchors_couples_SCT[,c("RI_var_strict", "RI_var_int")])
anchors_couples_SCT$RI_mean_eletp <- rowMeans(anchors_couples_SCT[,c("RI_eletp_strict", "RI_eletp_int")])

#ggplot(anchors_couples_SCT, aes (RI_mean_strict, RI_mean)) + geom_jitter()

saveRDS(anchors_couples_SCT,"C:/.../Results_obj/AR_SR_anchors_couples_SCT")
anchors_couples_SCT <- readRDS("C:/.../Results_obj/AR_SR_anchors_couples_SCT")

anchors_couples_top <- anchors_couples_SCT %>% group_by(cell2) %>% top_n(1, RI_mean_strict)

DimPlot(integrated,group.by = "subclass", cells = c(anchors_couples_top$cell1_name, anchors_couples_top$cell2_name) )

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

DimPlot(integrated,group.by = "subclass", cells = c(NN_couples$query_cells, NN_couples$ref_cells))

NN_couples$label <- integrated@meta.data[NN_couples$ref_cells,][["subclass"]]

NN_couples <- NN_couples %>% rowwise()%>% mutate(RI = cell_similarity(as.character(ref_mat[,ref_cells]),as.character(query_mat[,query_cells])))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_var = cell_similarity_variable(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = var_features))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_eletp = cell_similarity_variable(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = elettrophys.genes$V1))
#NN_couples <- NN_couples %>% rowwise()%>% mutate(RI_rf = cell_similarity_level(as.character(ref_mat[,ref_cells]),as.character(query_mat[,query_cells]),genes_mean, genes_mean_query))

NN_couples <- NN_couples %>% rowwise()%>% mutate(RI_strict = cell_similarity_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells])))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_var_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = var_features))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_eletp_strict = cell_similarity_variable_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]), features = elettrophys.genes$V1))
#NN_couples <- NN_couples %>% rowwise()%>% mutate(RI_strict_rf = cell_similarity_level_hard(as.data.frame(ref_mat[,ref_cells]),as.data.frame(query_mat[,query_cells]),genes_mean, genes_mean_query))

NN_couples <- NN_couples %>% rowwise()%>% mutate(RI_int = cell_similarity(as.character(int_mat[,ref_cells]),as.character(int_mat[,query_cells])))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_var_int = cell_similarity_variable(as.data.frame(int_mat[,ref_cells]),as.data.frame(int_mat[,query_cells]), features = var_features))
NN_couples <- NN_couples %>% rowwise() %>% mutate(RI_eletp_int = cell_similarity_variable(as.data.frame(int_mat[,ref_cells]),as.data.frame(int_mat[,query_cells]), features = elettrophys.genes$V1))


NN_couples$RI_mean <- rowMeans(NN_couples[,c("RI", "RI_var", "RI_eletp","RI_strict", "RI_var_strict", "RI_eletp_strict","RI_int", "RI_var_int", "RI_eletp_int")])
NN_couples$RI_mean_loose <- rowMeans(NN_couples[,c("RI", "RI_var", "RI_eletp")])
NN_couples$RI_mean_strict <- rowMeans(NN_couples[,c("RI_strict", "RI_var_strict", "RI_eletp_strict")])
NN_couples$RI_mean_int <- rowMeans(NN_couples[,c("RI_int", "RI_var_int", "RI_eletp_int")])
NN_couples$RI_mean_var <- rowMeans(NN_couples[,c("RI_var_strict", "RI_var_int")])
NN_couples$RI_mean_eletp <- rowMeans(NN_couples[,c("RI_eletp_strict", "RI_eletp_int")])

saveRDS(NN_couples,"C:/.../Results_obj/AR_SR_NN_couples")

####### plotting of anchors and NN #######

ac <- anchors_couples_top[,c("RI_mean","RI_mean_loose","RI_mean_strict","RI_mean_int","RI_var_strict","RI_var_int","RI_mean_var","RI_eletp_strict","RI_eletp_int","RI_mean_eletp" )]
ac$type <- "Anchors"
nc <- NN_couples[,c("RI_mean","RI_mean_loose","RI_mean_strict","RI_mean_int","RI_var_strict","RI_var_int","RI_mean_var","RI_eletp_strict","RI_eletp_int","RI_mean_eletp" )]
nc$type <- "NN couples"

saveRDS(NN_couples,"C:/.../Results_obj/SR_SR_NN_couples")
saveRDS(anchors_couples_SCT,"C:/.../Results_obj/SR_SR_anchors_couples_SCT")

results <- rbind(ac , nc)
ggplot(results, aes (RI_mean_strict, RI_mean, color = type)) + geom_point(size = 4) + labs(x = expression("RI"[mean]^s), y = expression("RI"[mean]), color = "Couples", title = "SubAllen vs Allen")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.title=element_text(size=20), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text = element_text(size=15))

ggplot(results, aes (RI_mean_int, RI_mean, color = type)) + geom_point(size = 4) + labs(x = expression("RI"[mean]^int), y = expression("RI"[mean]), color = "Couples", title = "SubAllen vs Allen")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.title=element_text(size=20), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text = element_text(size=15))

ggplot(results, aes (RI_mean_var, RI_mean_eletp, color = type)) + geom_point(size = 4) + labs(x = expression("RI"[var_mean]), y = expression("RI"[eph_mean]), color = "Couples", title = "SubAllen vs Allen")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), legend.title=element_text(size=20), legend.text=element_text(size=18), axis.title = element_text(size=20), axis.text = element_text(size=15))

mean(ac$RI_mean)
mean(nc$RI_mean)

mean(ac$RI_mean_strict)
mean(nc$RI_mean_strict)


ggplot(results, aes (RI_mean_int, RI_mean, color = type)) + geom_jitter()
ggplot(results, aes (RI_mean_strict, RI_mean_int, color = type)) + geom_jitter()

ggplot(results, aes (RI_mean_var, RI_mean_eletp, color = type)) + geom_jitter()

ggplot(results, aes (RI_var_strict, RI_eletp_strict, color = type)) + geom_jitter()
ggplot(results, aes (RI_var_int, RI_eletp_int, color = type)) + geom_jitter()
