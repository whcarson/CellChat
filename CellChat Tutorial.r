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

LCMA01_BEVR_CD31_MTX <- ReadMtx(
    mtx = paste(directory_string, "matrix.mtx.gz", sep = ""),
    features = paste(directory_string, "features.tsv.gz", sep = ""),
    cells = paste(directory_string, "barcodes.tsv.gz", sep = ""))


LCMA01_BEVR_CD31 <- CreateSeuratObject(counts = LCMA01_BEVR_CD31_MTX)

LCMA01_BEVR_CD31[["percent.mt"]] <- PercentageFeatureSet(LCMA01_BEVR_CD31, pattern="^MT-")
LCMA01_BEVR_CD31 <- subset(LCMA01_BEVR_CD31, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

LCMA01_BEVR_CD31 <- NormalizeData(LCMA01_BEVR_CD31, normalization.method = "LogNormalize", scale.factor = 10000)

LCMA01_BEVR_CD31 <- FindVariableFeatures(LCMA01_BEVR_CD31, selection.method = "vst", nfeatures=2000)

all.genes <- rownames(LCMA01_BEVR_CD31)

LCMA01_BEVR_CD31 <- ScaleData(LCMA01_BEVR_CD31, block.size = 1000)

LCMA01_BEVR_CD31 <- RunPCA(LCMA01_BEVR_CD31,  features=VariableFeatures(object=LCMA01_BEVR_CD31))
print(LCMA01_BEVR_CD31[["pca"]], dims=1:15, nfeatures = 15)

LCMA01_BEVR_CD31 <- FindNeighbors(LCMA01_BEVR_CD31, dims=1:10)
LCMA01_BEVR_CD31 <- FindClusters(LCMA01_BEVR_CD31, resolution=0.5)

LCMA01_BEVR_CD31 <- RunUMAP(LCMA01_BEVR_CD31, dims=1:10)
DimPlot(LCMA01_BEVR_CD31, reduction="umap")

levels(LCMA01_BEVR_CD31$seurat_clusters) <- c(
  "A", "B", "C", "D", "E", "F", "G",
  "H", "I", "J", "K", "L", "M")
LCMA01_BEVR_CD31_Chat <- createCellChat(object = LCMA01_BEVR_CD31, group.by="seurat_clusters")





"Set ligand-receptor interaction db"
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
LCMA01_BEVR_CD31_Chat@DB <- CellChatDB.use



# subset the expression data of signaling genes for saving computation cost
LCMA01_BEVR_CD31_Chat <- subsetData(LCMA01_BEVR_CD31_Chat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
LCMA01_BEVR_CD31_Chat <- identifyOverExpressedGenes(LCMA01_BEVR_CD31_Chat)
LCMA01_BEVR_CD31_Chat <- identifyOverExpressedInteractions(LCMA01_BEVR_CD31_Chat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)




LCMA01_BEVR_CD31_Chat <- computeCommunProb(LCMA01_BEVR_CD31_Chat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
LCMA01_BEVR_CD31_Chat <- filterCommunication(LCMA01_BEVR_CD31_Chat, min.cells = 10)




