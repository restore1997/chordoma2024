rm(list = ls())
library(Seurat)
library(infercnv)
library(data.table)

load("1-total.rda")
table(total$sample)

expression_matrix <- GetAssayData(chordoma, slot = "counts")

cell_types <- chordoma$MM
table(cell_types)
cell_names <- colnames(chordoma)
annotations_df <- data.frame(cell_names, cell_types)
annotations_df2 <- annotations_df[,"cell_types"]
annotations_df2 <- as.data.frame(annotations_df2)
rownames(annotations_df2) <- rownames(annotations_df)
colnames(annotations_df) <- "V2"

gene_order_file <- fread("gene_order_file.txt")
gene_order_file <- as.data.frame(gene_order_file)
rownames(gene_order_file) <- gene_order_file$V1
gene_order_file <- gene_order_file[,-1]

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = expression_matrix,
  annotations_file = annotations_df2,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = c("B_cell","T_cell","Myeloid","Neutrophil","Endothelial","Fibroblast")
)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="chordoma",
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE) #out_dir=tempfile()
save(infercnv_obj,file = "1-infercnv_obj_chordoma.rda")