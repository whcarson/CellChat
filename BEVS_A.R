devtools::install_github("sqjin/CellChat")
install.packages("NMF")
install.packages(biocon)
install.packages("dplyr")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(SeuratObject)
BiocManager::install("Biobase")  


BEVS_A.data <- ReadMtx(mtx = "BEVS_A_matrix.mtx.gz", features = "BEVS_A_features.tsv.gz", cells = "BEVS_A_barcodes.tsv.gz")
BEVS_A <- CreateSeuratObject(counts = BEVS_A.data)

## Will's Fix

##Seurat Data Prep

BEVS_A[["percent.mt"]] <- PercentageFeatureSet(BEVS_A, pattern="^MT-")
BEVS_A <- subset(BEVS_A, subset= nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
BEVS_A <- NormalizeData(BEVS_A, normalization.method = "LogNormalize", scale.factor = 10000)
BEVS_A <- FindVariableFeatures(BEVS_A, selection.method = "vst", nfeatures=2000)
all.genes <- rownames(BEVS_A)
BEVS_A <- ScaleData(BEVS_A, block.size = 1000000)

##Create Clusters
BEVS_A <- RunPCA(BEVS_A,  features=VariableFeatures(object=BEVS_A))
print(BEVS_A[["pca"]], dims=1:15, nfeatures = 15)
BEVS_A <- FindNeighbors(BEVS_A, dims=1:10)
BEVS_A <- FindClusters(BEVS_A, resolution=0.5)
BEVS_A <- RunUMAP(BEVS_A, dims=1:10)
DimPlot(BEVS_A, reduction="umap")
BEVS_A$seurat_clusters <- as.factor(as.numeric(as.character(BEVS_A$seurat_clusters)) + 1)


##Create a cell chat

BEVS_A_Chat <- createCellChat(object = BEVS_A, group.by="seurat_clusters")




"Set ligand-receptor interaction db"
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
BEVS_A_Chat@DB <- CellChatDB.use


# subset the expression data of signaling genes for saving computation cost
BEVS_A_Chat <- subsetData(BEVS_A_Chat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
BEVS_A_Chat <- identifyOverExpressedGenes(BEVS_A_Chat)
BEVS_A_Chat <- identifyOverExpressedInteractions(BEVS_A_Chat)
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


BEVS_A_Chat <- computeCommunProb(BEVS_A_Chat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
BEVS_A_Chat <- filterCommunication(BEVS_A_Chat, min.cells = 10)
## Infer the cell-cell communication at a signaling pathway level

BEVS_A_Chat <- computeCommunProbPathway(BEVS_A_Chat)

## Calculate the aggregated cell-cell communication netowrk

BEVS_A_Chat <- aggregateNet(BEVS_A_Chat)

## Vsiualization of cell-cell communication network

print(BEVS_A_Chat@netP$pathways)

## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(BEVS_A_Chat, sources.use = c(1:10), targets.use = c(1:10), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(BEVS_A_Chat, sources.use = 4, targets.use = c(1:30), signaling = BEVS_A_Chat@netP$pathways, remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) based on user's input (defined by pairLR.use)
pairLR.use <- extractEnrichedLR(BEVS_A_Chat, signaling = BEVS_A_Chat@netP$pathways)
netVisual_bubble(BEVS_A_Chat, sources.use = c(3,4), targets.use = c(1:50), pairLR.use = pairLR.use, remove.isolate = TRUE)

## Chord Plot
netVisual_chord_gene(BEVS_A_Chat, sources.use = 4, targets.use = c(5:11), signaling = BEVS_A_Chat@netP$pathways,legend.pos.x = 8)

## Gene Expression Plot
plotGeneExpression(BEVS_A_Chat, signaling = BEVS_A_Chat@netP$pathways, enriched.only = FALSE)

BEVS_A_groupSize <- as.numeric(table(BEVS_A_Chat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(BEVS_A_Chat@net$count, vertex.weight = BEVS_A_groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(BEVS_A_Chat@net$weight, vertex.weight = BEVS_A_groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



#Visualize MHC-II pathway interactions
netVisual_circle(BEVS_A_Chat@netP$prob[,,"SPP1"])


#Generate Plot
par(mfrow=c(6,5))
par(mai=c(0.1, 0.2, 0.2, 0.1))
for (pathway in BEVS_A_Chat@netP$pathways){
  print(pathway)
  netVisual_circle(BEVS_A_Chat@netP$prob[,,pathway], title.name = pathway)
}

singaling.path <- c("SPP1")

netVisual_heatmap(BEVS_A_Chat, signaling = , color.heatmap = "Reds")
