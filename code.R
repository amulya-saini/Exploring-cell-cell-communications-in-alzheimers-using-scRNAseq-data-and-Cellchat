# devtools::install_github("jinworks/CellChat")
# devtools::install_github("jokergoo/ComplexHeatmap")
# .libPaths()
# devtools::install_github("jokergoo/circlize")
# install.packages("digest")
# install.packages('NMF')
# install.packages("presto")
# packages <- c("curl", "glue", "htmltools", "processx", "Rcpp", "promises", "httpuv")
# install.packages(packages)
# install.packages("rlang")
# devtools::install_github("immunogenomics/presto")
# BiocManager::install("Biobase")
# BiocManager::install("BiocNeighbors")
# BiocManager::install("BiocGenerics")
# BiocManager::install('glmGamPoi')
library(CellChat)
library(patchwork)
library(Matrix)
library(Seurat)
library(dplyr)
library(tidyr)

# we are not integrating the data here because we want to look at the cell-cell communication 
# in alzheimer's 

############### FAD (Control) #######################

# Loading the familial Alzheimer's disease dataset that has not been treated with APC
# APC (Activated protein C)
fad.dir <- "C:/Users/saini/Downloads/project_systems/FAD/"
# Load the matrix file using Read10X with gzipped files
fad_cnt_mtx <- Read10X(fad.dir, 
                       gene.column = 2, 
                       cell.column = 1, 
                       unique.features = TRUE, 
                       strip.suffix = FALSE)

# creating a seurat object of fad
fad = CreateSeuratObject(counts = fad_cnt_mtx, 
                         project = "FAD", 
                         min.cells = 3, 
                         min.features = 200)

# creating a new column to store the percentage of mitochondrial counts originating from a set of features
fad[["percent.mt"]] <- PercentageFeatureSet(fad, 
                                            pattern = "^mt-")

# Show QC metrics for the first 5 cells 
head(fad@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(fad, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships.
# fad
plot1 <- FeatureScatter(fad, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(fad, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2

# filtering the dataset
fad <- subset(fad, 
              subset = nFeature_RNA > 200 & 
                nFeature_RNA < 4500 & 
                percent.mt < 25)

# Normalizing the data
fad <- NormalizeData(fad, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)
fad[["RNA"]]$data

# finding variable genes
fad <- FindVariableFeatures(fad, 
                            selection.method = "vst", 
                            nfeatures = 2000)

# Identify the 10 most highly variable genes
fad_top10 <- head(VariableFeatures(fad), 10)

# plot variable features with labels
plot1 <- VariableFeaturePlot(fad)
plot2 <- LabelPoints(plot = plot1, 
                     points = fad_top10, 
                     repel = TRUE)
plot2

# Scaling the data
fad.all.genes <- rownames(fad)
fad <- ScaleData(fad, 
                 features = fad.all.genes)

# Perform linear dimensional reduction
fad <- RunPCA(fad, 
              features = VariableFeatures(object = fad))
ElbowPlot(fad)

# finding neighbors by constructing a KNN graph based on the euclidean distance in PCA space
fad <- FindNeighbors(fad, 
                     dims = 1:20)

# grouping cells together
fad <- FindClusters(fad, 
                    resolution = 0.6)
head(Idents(fad), 5)

# Runnning non-linear dimensionality reduction (UMAP)
fad <- RunUMAP(fad, dims = 1:20, 
               # umap.method = 'umap-learn', 
               # metric = 'correlation', 
               verbose = FALSE)
DimPlot(fad, 
        reduction = "umap",
        label =T)

# find markers for every cluster compared to all remaining cells, report only the positive ones
fad_markers <- FindAllMarkers(fad, 
                                  only.pos = TRUE)

fad_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# fad_markers %>%
#   dplyr::filter(cluster == 26) %>%
#   dplyr::filter(avg_log2FC > 1)

# storing all the markers with 1 log fold change in each cluster
# into a dataframe

# Initialize an empty list to store filtered dataframes
filtered_data <- list()

# Loop through clusters 0 to 21
for (i in 0:21) {
  # Filter fad_markers for the current cluster and avg_log2FC > 1
  filtered <- fad_markers %>%
    filter(cluster == i, avg_log2FC > 1)
  
  # Store the filtered dataframe in the list
  filtered_data[[i+1]] <- filtered
}

filtered_data[[22]]

# Cell Type annotations
fad <- RenameIdents(fad, "0" = "Microglia")
fad <- RenameIdents(fad, "1" = "Microglia")
fad <- RenameIdents(fad, "2" = "Microglia")
fad <- RenameIdents(fad, "3" = "Astrocytes")
fad <- RenameIdents(fad, "4" = "Endothelial cells")
fad <- RenameIdents(fad, "5" = "Oligodendrocytes")
fad <- RenameIdents(fad, "6" = "Neurons")
fad <- RenameIdents(fad, "7" = "Ependymal cells")
fad <- RenameIdents(fad, "8" = "Astrocytes")
fad <- RenameIdents(fad, "9" = "T memory cells")
fad <- RenameIdents(fad, "10" = "Oligodendrocytes")
fad <- RenameIdents(fad, "11" = "Macrophages")
fad <- RenameIdents(fad, "12" = "Unknown")
fad <- RenameIdents(fad, "13" = "Pericytes")
fad <- RenameIdents(fad, "14" = "Unknown")
fad <- RenameIdents(fad, "15" = "Neurons")
fad <- RenameIdents(fad, "16" = "Oligodendrocytes")
fad <- RenameIdents(fad, "17" = "B cells")
fad <- RenameIdents(fad, "18" = "Pericytes")
fad <- RenameIdents(fad, "19" = "Unknown")
fad <- RenameIdents(fad, "20" = "Fibroblasts")
fad <- RenameIdents(fad, "21" = "Neutrophils")

DimPlot(fad, 
        reduction = "umap",
        label = TRUE)

###################### Cellchat for FAD sample

# preparing seurat data for cellchat
data.input_fad <- fad[["RNA"]]$data
labels_fad <- Idents(fad)
meta_fad <- data.frame(labels = labels_fad, row.names = names(labels_fad))

# creating a cellchat object
fad_cellChat <- createCellChat(object = fad, group.by = "ident", assay = "RNA")

# Load the mouse CellChat database
CellChatDB <- CellChatDB.mouse 

# Add metadata to the CellChat object
fad_cellChat <- addMeta(fad_cellChat, meta = meta_fad)

# Set cell identities in the CellChat object
fad_cellChat <- setIdent(fad_cellChat, ident.use = "group")

# Show factor levels of the cell labels (group identities)
levels(fad_cellChat@idents)

# Calculate the number of cells in each cell group
groupSize <- as.numeric(table(fad_cellChat@idents))

# Assign the CellChat database to the CellChat object
fad_cellChat@DB <- CellChatDB

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# Subset the data in the CellChat object (this step is necessary even if using the whole database)
fad_cellChat <- subsetData(fad_cellChat)

# Identify overexpressed genes in the data
fad_cellChat <- identifyOverExpressedGenes(fad_cellChat)

# Identify overexpressed interactions in the data
fad_cellChat <- identifyOverExpressedInteractions(fad_cellChat)

# Project the data onto a protein-protein interaction (PPI) network specific to mouse
fad_cellChat <- projectData(fad_cellChat, PPI.mouse)

# Compute the communication probabilities between cell types using the projected data
fad_cellChat <- computeCommunProb(fad_cellChat, raw.use = FALSE)  # Use the projected data

# Filter communication interactions based on a minimum number of cells
fad_cellChat <- filterCommunication(fad_cellChat, min.cells = 10)

# Aggregate the communication network
fad_cellChat <- aggregateNet(fad_cellChat)

# Set up the plotting environment
par(mfrow = c(1, 2), xpd = TRUE)

# Visualize the number of interactions between any two cell groups using a circle plot
netVisual_circle(fad_cellChat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")

# Visualize the interaction weights/strength between any two cell groups using a circle plot
netVisual_circle(fad_cellChat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")

# Visualize interactions for each cell group individually using circle plots
mat <- fad_cellChat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Define the signaling pathways to be visualized
pathways.show <- c("ApoE")

# Visualize the hierarchy plot for the signaling pathways
vertex.receiver = seq(1, 4)
netVisual_aggregate(fad_cellChat, signaling = pathways.show, vertex.receiver = vertex.receiver)

# Visualize the circle plot for the signaling pathways
par(mfrow = c(1, 1))
netVisual_aggregate(fad_cellChat, signaling = pathways.show, layout = "circle")

# Visualize the chord plot for the signaling pathways
par(mfrow = c(1, 1))
netVisual_aggregate(fad_cellChat, signaling = pathways.show, layout = "chord")

# Visualize the heatmap for the signaling pathways
par(mfrow = c(1, 1))
netVisual_heatmap(fad_cellChat, signaling = pathways.show, color.heatmap = "Reds")

# Compute contribution of signaling pathways to the network
netAnalysis_contribution(fad_cellChat, signaling = pathways.show)

# Extract enriched ligand-receptor pairs for the signaling pathways
pairLR.ApoE <- extractEnrichedLR(fad_cellChat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.ApoE[1,]

# Visualize hierarchy plot for the individual ligand-receptor pair
netVisual_individual(fad_cellChat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Visualize circle plot for the individual ligand-receptor pair
netVisual_individual(fad_cellChat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Visualize chord plot for the individual ligand-receptor pair
netVisual_individual(fad_cellChat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Visualize all significant interactions from certain cell groups to other cell groups using bubble plot
netVisual_bubble(fad_cellChat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# Plot the gene expression distribution for the signaling pathway using violin plot
plotGeneExpression(fad_cellChat, signaling = "ApoE", enriched.only = TRUE, type = "violin")

# Compute centrality measures for the inferred communication network
fad_cellChat <- netAnalysis_computeCentrality(fad_cellChat, slot.name = "netP")

# Analyze signaling roles in the network for the signaling pathway
netAnalysis_signalingRole_network(fad_cellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


############### FAD APC (Condition) #######################

# Loading the familial Alzheimer's disease dataset that has been treated with APC
fad_apc.dir <- "C:/Users/saini/Downloads/project_systems/FAD_APC/"
# Load the matrix file using Read10X with gzipped files
fadapc_cnt_mtx <- Read10X(fad_apc.dir, 
                       gene.column = 2, 
                       cell.column = 1, 
                       unique.features = TRUE, 
                       strip.suffix = FALSE)

# creating a seurat object of fad_apc
fad_apc = CreateSeuratObject(counts = fadapc_cnt_mtx, 
                         project = "FAD_APC", 
                         min.cells = 3, 
                         min.features = 200)

# creating a new column to store the percentage of mitochondrial counts originating from a set of features
fad_apc[["percent.mt"]] <- PercentageFeatureSet(fad_apc, 
                                            pattern = "^mt-")

# Show QC metrics for the first 5 cells 
head(fad_apc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(fad_apc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationshi
# fad_apc
plot3 <- FeatureScatter(fad_apc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot4 <- FeatureScatter(fad_apc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot3 + plot4

# filtering the dataset
fad_apc <- subset(fad_apc, 
              subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 & 
                percent.mt < 25)

# Normalizing the data
fad_apc <- NormalizeData(fad_apc, 
                    normalization.method = "LogNormalize", 
                    scale.factor = 10000)

fad_apc[["RNA"]]$data

# finding variable genes
fad_apc <- FindVariableFeatures(fad_apc, 
                            selection.method = "vst", 
                            nfeatures = 2000)

# Identify the 10 most highly variable genes
fad_apc_top10 <- head(VariableFeatures(fad_apc), 10)

# plot variable features with labels
plot3 <- VariableFeaturePlot(fad_apc)
plot4 <- LabelPoints(plot = plot3, 
                     points = fad_apc_top10, 
                     repel = TRUE)
plot4

# Scaling the data
fad_apc.all.genes <- rownames(fad_apc)
fad_apc <- ScaleData(fad_apc, 
                features = fad_apc.all.genes)

# wt[["RNA"]]$scale.data (don't print)

# Perform linear dimensional reduction
fad_apc <- RunPCA(fad_apc, 
             features = VariableFeatures(object = fad_apc))
ElbowPlot(fad_apc)

# finding neighbors by constructing a KNN graph based on the euclidean distance in PCA space
fad_apc <- FindNeighbors(fad_apc, 
                     dims = 1:17)

# grouping cells together
fad_apc <- FindClusters(fad_apc, 
                    resolution = 0.3)
head(Idents(fad_apc), 5)

# Runnning non-linear dimensionality reduction (UMAP)
fad_apc <- RunUMAP(fad_apc, dims = 1:20, 
              # umap.method = 'umap-learn', 
              # metric = 'correlation', 
              verbose = FALSE)
DimPlot(fad_apc, 
        reduction = "umap",
        label =T)

# find markers for every cluster compared to all remaining cells, report only the positive ones
fad_apc_markers <- FindAllMarkers(fad_apc, 
                     only.pos = TRUE)

fad_apc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

fad_apc_markers %>%
  dplyr::filter(cluster == 23) %>%
  dplyr::filter(avg_log2FC > 1)

# Cell Type annotations
fad_apc <- RenameIdents(fad_apc, "0" = "Microglia")
fad_apc <- RenameIdents(fad_apc, "1" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "2" = "Astrocytes")
fad_apc <- RenameIdents(fad_apc, "3" = "Endothelial cells")
fad_apc <- RenameIdents(fad_apc, "4" = "Choroid plexus cells")
fad_apc <- RenameIdents(fad_apc, "5" = "Neurons")
fad_apc <- RenameIdents(fad_apc, "6" = "Neurons")
fad_apc <- RenameIdents(fad_apc, "7" = "Ependymal cells")
fad_apc <- RenameIdents(fad_apc, "8" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "9" = "Pericytes")
fad_apc <- RenameIdents(fad_apc, "10" = "Macrophages")
fad_apc <- RenameIdents(fad_apc, "11" = "Microglia")
fad_apc <- RenameIdents(fad_apc, "12" = "Fibroblasts")
fad_apc <- RenameIdents(fad_apc, "13" = "Neurons")
fad_apc <- RenameIdents(fad_apc, "14" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "15" = "T memory cells")
fad_apc <- RenameIdents(fad_apc, "16" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "17" = "Endothelial cells")
fad_apc <- RenameIdents(fad_apc, "18" = "Pericytes")
fad_apc <- RenameIdents(fad_apc, "19" = "Endothelial cells")
fad_apc <- RenameIdents(fad_apc, "20" = "Neutrophils")
fad_apc <- RenameIdents(fad_apc, "21" = "B cells")
fad_apc <- RenameIdents(fad_apc, "22" = "Microglia")
fad_apc <- RenameIdents(fad_apc, "23" = "Unknown")

DimPlot(fad_apc, 
        reduction = "umap",
        label =T)

########################### cellchat for FAD-APC

# preparing seurat data for cellchat
data.input_fad_apc <- fad_apc[["RNA"]]$data
labels_fad_apc <- Idents(fad_apc)
meta_fad_apc <- data.frame(labels = labels_fad_apc, row.names = names(labels_fad_apc))

# creating a cellchat object
fad_apc_cellChat <- createCellChat(object = fad_apc, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# set the used database in the object
fad_apc_cellChat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
fad_apc_cellChat <- subsetData(fad_apc_cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

# Identify overexpressed genes in the data
fad_apc_cellChat <- identifyOverExpressedGenes(fad_apc_cellChat)

# Identify overexpressed interactions in the data
fad_apc_cellChat <- identifyOverExpressedInteractions(fad_apc_cellChat)

# Project the data onto a protein-protein interaction (PPI) network specific to mouse
fad_apc_cellChat <- projectData(fad_apc_cellChat, PPI.mouse)

# Compute the communication probabilities between cell types using the projected data
fad_apc_cellChat <- computeCommunProb(fad_apc_cellChat, raw.use = FALSE)  # Use the projected data

# Filter communication interactions based on a minimum number of cells
fad_apc_cellChat <- filterCommunication(fad_apc_cellChat, min.cells = 10)

# Aggregate the communication network
fad_apc_cellChat <- aggregateNet(fad_apc_cellChat)

# Set up the plotting environment
par(mfrow = c(1, 2), xpd = TRUE)

# Visualize the number of interactions between any two cell groups using a circle plot
netVisual_circle(fad_apc_cellChat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")

# Visualize the interaction weights/strength between any two cell groups using a circle plot
netVisual_circle(fad_apc_cellChat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")

# Visualize interactions for each cell group individually using circle plots
mat <- fad_apc_cellChat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Define the signaling pathways to be visualized
pathways.show <- c("ApoE")

# Visualize the hierarchy plot for the signaling pathways
vertex.receiver = seq(1, 4)
netVisual_aggregate(fad_apc_cellChat, signaling = pathways.show, vertex.receiver = vertex.receiver)

# Visualize the circle plot for the signaling pathways
par(mfrow = c(1, 1))
netVisual_aggregate(fad_apc_cellChat, signaling = pathways.show, layout = "circle")

# Visualize the chord plot for the signaling pathways
par(mfrow = c(1, 1))
netVisual_aggregate(fad_apc_cellChat, signaling = pathways.show, layout = "chord")

# Visualize the heatmap for the signaling pathways
par(mfrow = c(1, 1))
netVisual_heatmap(fad_apc_cellChat, signaling = pathways.show, color.heatmap = "Reds")

# Compute contribution of signaling pathways to the network
netAnalysis_contribution(fad_apc_cellChat, signaling = pathways.show)

# Extract enriched ligand-receptor pairs for the signaling pathways
pairLR.ApoE <- extractEnrichedLR(fad_apc_cellChat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.ApoE[1,]

# Visualize hierarchy plot for the individual ligand-receptor pair
netVisual_individual(fad_apc_cellChat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Visualize circle plot for the individual ligand-receptor pair
netVisual_individual(fad_apc_cellChat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Visualize chord plot for the individual ligand-receptor pair
netVisual_individual(fad_apc_cellChat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# Visualize all significant interactions from certain cell groups to other cell groups using bubble plot
netVisual_bubble(fad_apc_cellChat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)

# Plot the gene expression distribution for the signaling pathway using violin plot
plotGeneExpression(fad_apc_cellChat, signaling = "ApoE", enriched.only = TRUE, type = "violin")

# Compute centrality measures for the inferred communication network
fad_apc_cellChat <- netAnalysis_computeCentrality(fad_apc_cellChat, slot.name = "netP")

# Analyze signaling roles in the network for the signaling pathway
netAnalysis_signalingRole_network(fad_apc_cellChat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
