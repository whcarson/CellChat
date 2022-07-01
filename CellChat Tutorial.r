rm(list = ls())
library(CellChat)
library(patchwork)
library(dplyr)
library(Seurat)
options(stringsAsFactors = FALSE)

lcma01_bevr_cd31_mtx <- ReadMtx(
    mtx = "../../UChicago/Internships etc/Internship Opportunities/
    Aghi Lab/CD31 sorted/lcma01_bevr_cd31/filtered_feature_bc_matrix/
    matrix.mtx.gz",
    features = "../../UChicago/Internships etc/Internship Opportunities/
    Aghi Lab/CD31 sorted/lcma01_bevr_cd31/filtered_feature_bc_matrix/
    features.tsv.gz", 
    cells = "../../UChicago/Internships etc/Internship Opportunities/
    Aghi Lab/CD31 sorted/lcma01_bevr_cd31/filtered_feature_bc_matrix/
    barcodes.tsv.gz")

lcma01_bevr_cd31 <- CreateSeuratObject(counts = lcma01_bevr_cd31_mtx)

lcma01_bevr_cd31[["percent.mt"]] <- PercentageFeatureSet(lcma01_bevr_cd31, pattern="^MT-")
lcma01_bevr_cd31 <- subset(lcma01_bevr_cd31, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

lcma01_bevr_cd31 <- NormalizeData(lcma01_bevr_cd31, normalization.method = "LogNormalize", scale.factor = 10000)

lcma01_bevr_cd31 <- FindVariableFeatures(lcma01_bevr_cd31, selection.method = "vst", nfeatures=2000)

all.genes <- rownames(lcma01_bevr_cd31)

lcma01_bevr_cd31 <- ScaleData(lcma01_bevr_cd31, block.size = 1000000)

lcma01_bevr_cd31 <- RunPCA(lcma01_bevr_cd31,  features=VariableFeatures(object=lcma01_bevr_cd31))
print(lcma01_bevr_cd31[["pca"]], dims=1:15, nfeatures = 15)

lcma01_bevr_cd31 <- FindNeighbors(lcma01_bevr_cd31, dims=1:10)
lcma01_bevr_cd31 <- FindClusters(lcma01_bevr_cd31, resolution=0.5)

lcma01_bevr_cd31 <- RunUMAP(lcma01_bevr_cd31, dims=1:10)
DimPlot(lcma01_bevr_cd31, reduction="umap")

lcma01_bevr_cd31_Chat <- createCellChat(object = lcma01_bevr_cd31, group.by="seurat_clusters")


"Set ligand-receptor interaction db"
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
lcma01_bevr_cd31_Chat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
lcma01_bevr_cd31_Chat <- subsetData(lcma01_bevr_cd31_Chat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
lcma01_bevr_cd31_Chat <- identifyOverExpressedGenes(lcma01_bevr_cd31_Chat)
lcma01_bevr_cd31_Chat <- identifyOverExpressedInteractions(lcma01_bevr_cd31_Chat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

lcma01_bevr_cd31_Chat <- computeCommunProb(lcma01_bevr_cd31_Chat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
lcma01_bevr_cd31_Chat <- filterCommunication(lcma01_bevr_cd31_Chat, min.cells = 10)



