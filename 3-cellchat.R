rm(list = ls())
library(Seurat)
library(CellChat)
library(patchwork)

load("1-total.rda")
total

rm(list = setdiff(ls(), "total"))

data.input <- GetAssayData(total, assay = "RNA", slot = "data") # normalized data matrix
meta <- total@meta.data

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "MM2")
rm(list = setdiff(ls(), "cellchat"))

#> Create a CellChat object from a data matrix
#> Set cell identities for the new CellChat object
#> The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT

cellchat <- setIdent(cellchat, ident.use = "MM2") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
table(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
CellChatDB_interaction <- CellChatDB$interaction
table(CellChatDB_interaction$annotation)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


##Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
sort(table(df.net$pathway_name))
sort(table(df.net$ligand))

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

save(cellchat,file = "1-cellchat.rda")