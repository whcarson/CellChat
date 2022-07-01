rm(list = ls())
library(CellChat)
library(patchwork)
library(dplyr)
library(Seurat)
options(stringsAsFactors = FALSE)

directory_string <- paste("C:/Users/19193/OneDrive - The University of Chicago",
                          "/Documents/UChicago/Internships etc/",
                          "Internship Opportunities/Aghi Lab/CD31 sorted/",
                          "LCMA01_BEVR_CD31/filtered_feature_bc_matrix/",
                          sep = "")
lcma01_bevr_cd31_mtx <- ReadMtx(
    mtx = paste(directory_string, "matrix.mtx.gz", sep = ""),
    features = paste(directory_string, "features.tsv.gz", sep = ""),
    cells = paste(directory_string, "barcodes.tsv.gz", sep = ""))

lcma01_bevr_cd31 <- CreateSeuratObject(counts = lcma01_bevr_cd31_mtx)

# lcma01_bevr_cd31[["percent.mt"]] <- PercentageFeatureSet(
#     lcma01_bevr_cd31, pattern = "^MT-")
# lcma01_bevr_cd31 <- subset(
#     lcma01_bevr_cd31, subset = nFeature_RNA > 200 &
#      nFeature_RNA < 2500 & percent.mt < 5)

lcma01_bevr_cd31 <- NormalizeData(
    lcma01_bevr_cd31, normalization.method = "LogNormalize",
    scale.factor = 10000)

lcma01_bevr_cd31 <- FindVariableFeatures(
    lcma01_bevr_cd31, selection.method = "vst", nfeatures = 2000)

all_genes <- rownames(lcma01_bevr_cd31)

lcma01_bevr_cd31 <- ScaleData(lcma01_bevr_cd31, block.size = 1000)

lcma01_bevr_cd31 <- RunPCA(
    lcma01_bevr_cd31,  features = VariableFeatures(object = lcma01_bevr_cd31))
print(lcma01_bevr_cd31[["pca"]], dims = 1:15, nfeatures = 15)

lcma01_bevr_cd31 <- FindNeighbors(lcma01_bevr_cd31, dims = 1:10)
lcma01_bevr_cd31 <- FindClusters(lcma01_bevr_cd31, resolution = 0.5)

lcma01_bevr_cd31 <- RunUMAP(lcma01_bevr_cd31, dims = 1:10)
DimPlot(lcma01_bevr_cd31, reduction = "umap")
lcma01_bevr_cd31$seurat_clusters <- as.factor(as.numeric(as.character(lcma01_bevr_cd31$seurat_clusters)) + 1)

lcma01_bevr_cd31_chat <- createCellChat(
    object = lcma01_bevr_cd31, group.by = "seurat_clusters")


"Set ligand-receptor interaction db"
cell_chat_db <- CellChatDB.human
showDatabaseCategory(cell_chat_db)
dplyr::glimpse(cell_chat_db$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(
    # CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
cell_chat_db_use <- cell_chat_db # simply use the default CellChatDB

# set the used database in the object
lcma01_bevr_cd31_chat@DB <- cell_chat_db_use

# subset the expression data of signaling genes for saving computation cost
lcma01_bevr_cd31_chat <- subsetData(lcma01_bevr_cd31_chat)
lcma01_bevr_cd31_chat <- identifyOverExpressedGenes(lcma01_bevr_cd31_chat)
lcma01_bevr_cd31_chat <- identifyOverExpressedInteractions(lcma01_bevr_cd31_chat)

# project gene expression data onto PP
# (Optional: when running it, USER should set `raw.use = FALSE`
# in the function `computeCommunProb()` in order to use the projected data)

lcma01_bevr_cd31_chat <- computeCommunProb(lcma01_bevr_cd31_chat)
# Filter out the cell-cell communication
# if there are only few number of cells in certain cell groups
lcma01_bevr_cd31_chat <- filterCommunication(
    lcma01_bevr_cd31_chat, min.cells = 10)

lcma01_bevr_cd31_df_lr <- subsetCommunication(lcma01_bevr_cd31_chat)
lcma01_bevr_cd31_df_path <- subsetCommunication(
    lcma01_bevr_cd31_chat, slot.name = "netP")

lcma01_bevr_cd31_chat <- aggregateNet(lcma01_bevr_cd31_chat)

unique(lcma01_bevr_cd31_chat@idents)