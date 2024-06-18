######################################################
# NORMALIZATION, DIMENSIONAL REDUCTION, CLUSTERING
#
# GOAL: For both the under tissue and piwipos datasets,
# SCTransform the data, run PCA, generate neighbors
# graph and clusters, then generate UMAP.
######################################################

setwd('/n/projects/cb2350/smed_slide-seq/repo/code/')

# ================= IMPORTS ==========================

library(Seurat)

# ================= PARAMS ===========================

nPCs <- list('undertissue'=55,'piwipos'=20)
res <- list('undertissue'=1,'piwipos'=0.8)
seeds <- list('undertissue'=42,'piwipos'=112358)

tissuelvs <- c("Phagocytic Cells","Epidermis","Intestine",
               "Muscle","Neural","Parenchymal","Pharynx",
               "Protonephridia","Stem Cells","Non-differentiated")

# ================= UNDER TISSUE ======================

meta <- openxlsx::read.xlsx('../data/GSE199348_bead_metadata.xlsx',rowNames=TRUE)

seu <- readRDS('../data/seurat_objects/under_tissue.rds')

seu <- SCTransform(seu,vst.flavor='v1',return.only.var.genes=FALSE)
seu <- RunPCA(seu,npcs=nPCs$undertissue)

# NOTE!! Despite my best efforts and vst.flavor='v1'
# I cannot get the current version of Seurat to replicate
# the SCTransform behavior of prior versions.
# This results in different clustering, so we resort
# to pulling cluster designations and UMAP embedding
# from the published metadata tables.
# If you'd like to try to sort out this minor mystery,
# just email and I'll send along the original SCTranformed
# expression data :)
seu$seurat_clusters <- droplevels(factor(meta[Cells(seu),'global_cluster'],levels=0:99))
seu$tissue <- factor(meta[Cells(seu),'tissue'],levels=tissuelvs)
umap <- as.matrix(meta[Cells(seu),c('global_umap_1','global_umap_2')])
colnames(umap) <- c('umap_1','umap_2')
seu[['umap']] <- CreateDimReducObject(umap)

saveRDS(seu,'../data/seurat_objects/under_tissue.rds')

# ================= PIWI1 POSITIVE SUBSET ============

seu <- readRDS('../data/seurat_objects/piwi1_positive_subset.rds')

seu <- SCTransform(seu)
set.seed(seeds$piwipos)
seu <- RunPCA(seu,npcs=nPCs$piwipos)

# Again with the SCTransform issue and the clustering
seu$seurat_clusters <- droplevels(factor(meta[Cells(seu),'piwipos_cluster'],levels=0:99))
seu$tissue <- factor(meta[Cells(seu),'tissue'],levels=tissuelvs)
umap <- as.matrix(meta[Cells(seu),c('stem_cell_umap_1','stem_cell_umap_2')])
colnames(umap) <- c('umap_1','umap_2')
seu[['umap']] <- CreateDimReducObject(umap)
saveRDS(seu,'../data/seurat_objects/piwi1_positive_subset.rds')
