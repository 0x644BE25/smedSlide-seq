######################################################
# CREATE SEURAT OBJECTS
#
# GOALS: Create a combined Seurat object with data from
# all four pucks. 
#
# Use provided puck coordinates to create a combined,
# non-overlapping "slide" embedding.
#
# Add manually annotated fragment designations, and
# use background and intact fragment designations
# to add timepoint information.
#
# Use distance from centroid of 60 nearest neigboring
# beads to determine beads on edges of fragments and
# categorize beads as "background", "edge", or 
# "interior."
#
# Create filtered version of Seurat object with only
# beads >= 40 nUMI that are under tissue fragments
# ("edge" or "interior").
#
# Use manually determined fragment angles to correct
# for angles and create aligned "rows" embedding.
#
# Filter by piwi-1 (SMED30007406) to generate piwi1
# positive subset.
######################################################

# ================= IMPORTS ==========================

library(Seurat)
library(openxlsx)

# ================= PARAMS ===========================

metaCols <- c('timepoint','puck','fragment','region')

timepoints <- c('L46130'='06hpa','L46131'='48hpa','L46144'='06hpa','L46145'='48hpa')

rows <- list('intact'=c('L46130-1','L46144-9','L46131-21','L46145-10'),
             '06hpa-1'=c('L46130-4','L46130-5','L46130-3','L46130-2','L46130-6','L46130-7','L46130-8'),
             '06hpa-2'=c('L46144-12','L46144-13','L46144-11','L46144-10','L46144-14','L46144-15','L46144-16','L46144-17','L46144-18','L46144-19','L46144-20'),
             '48hpa-1'=c('L46131-28','L46131-27','L46131-26','L46131-25','L46131-24','L46131-22','L46131-23','L46131-31','L46131-30','L46131-29'),
             '48hpa-2'=c('L46145-11','L46145-9','L46145-8','L46145-7','L46145-6','L46145-5','L46145-4','L46145-3','L46145-2','L46145-1'))

w <- 1100 # fragment spacing width
h <- 1900 # fragment spacing height

minUMIs <- 40

nCentroidNeighbors <- 60
edgeBeadCutoff <- 30

piwi1 <- 'SMED30007406'

# It is assumed that you have placed the metadata XLSX file
# and unzipped TAR directory into your-working-directory/data/
metaURL <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE199348&format=file&file=GSE199348_bead_metadata.xlsx'
tarURL <- 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE199348&format=file'

# ================= COMBINE INDIVIDUAL PUCKS =========

files <- dir('../data/GSE199348_RAW/')
gsms <- sort(unique(sapply(strsplit(files,'_'),function(x){ x[1] })))

# make seurat_object directory exists
if (!dir.exists('../data/seurat_objects')) { dir.create('../data/seurat_objects' ) }

# BUILD SEURAT FILE
seu <- NULL
allcoords <- NULL
fragments <- NULL
for (gsm in gsms) {
  # BATCH FILES
  countsFile <- files[startsWith(files,gsm) & endsWith(files,'matrix.csv.gz')]
  puckFile <- files[startsWith(files,gsm) & endsWith(files,'puck.txt.gz')]
  batch <- strsplit(countsFile,'_')[[1]][2]
  cat(batch,'\n')
  
  # PUCK COORDINATES
  puck <- read.table(paste0('../data/GSE199348_RAW/',puckFile),sep='\t',header=FALSE,row.names=1)
  rownames(puck) <- paste(batch,rownames(puck),sep='-')
  allcoords <- rbind(allcoords,puck)

  # COUNTS
  counts <- read.csv(paste0('../data/GSE199348_RAW/',countsFile),row.names=1)
  colnames(counts) <- paste(batch,colnames(counts),sep='-')
  counts <- counts[,intersect(colnames(counts),rownames(puck))]
  
  # SEURAT OBJECT
  curr <- suppressWarnings(CreateSeuratObject(counts=counts,project='smed_slide-seq'))
  if (is.null(seu)) { seu <- curr } else { seu <- merge(seu,curr) }
}

seu <- JoinLayers(seu)

# ================= SLIDE EMBEDDING ==================

colnames(allcoords) <- c('slide_1','slide_2')
emb <- allcoords[Cells(seu),]

# find params for slide grid
sliderange <- apply(emb,2,range) 
w <- ceiling(as.numeric(sliderange[2,'slide_1']-sliderange[1,'slide_1']))
h <- ceiling(as.numeric(sliderange[2,'slide_2']-sliderange[1,'slide_2']))
x <- ceiling(max(c(h,w))/100)*100

# move each puck into grid
emb[startsWith(rownames(emb),'L46130'),2] <- emb[startsWith(rownames(emb),'L46130'),2]+x
emb[startsWith(rownames(emb),'L46144'),] <- emb[startsWith(rownames(emb),'L46144'),]+x
emb[startsWith(rownames(emb),'L46145'),1] <- emb[startsWith(rownames(emb),'L46145'),1]+x

# dim reduc
seu[['slide']] <- CreateDimReducObject(embeddings=as.matrix(emb),key='slide_',assay='RNA')
seu$puck <- factor(substr(Cells(seu),1,6),levels=sort(unique(substr(Cells(seu),1,6))))

# ================= FRAGMENTS ========================

# more metadata
frags <- read.csv('../data/bead_fragment_designations.csv',row.names=1)
lvs <- paste(rep(levels(seu$puck),each=100),rep(0:99,4),sep='-')
seu$fragment <- droplevels(factor(frags[Cells(seu),'fragment'],levels=lvs))
seu$timepoint <- setNames(timepoints[as.character(seu$puck)],Cells(seu))
seu$timepoint[seu$fragment %in% rows$intact] <- 'intact'
seu$timepoint[!seu$fragment %in% unlist(rows)] <- NA
seu$timepoint <- factor(seu$timepoint,levels=c('intact','06hpa','48hpa'))

# identify edge beads
emb <- Embeddings(seu,reduction='slide')
UMIcounts <- NULL
for (f in levels(seu$fragment)) {
  beads <- Cells(seu)[seu$fragment==f]
  if (!f %in% unlist(rows)) {
    UMIcounts <- rbind(UMIcounts,data.frame(nUMI=seu$nCount_RNA[beads],fragment=f,group='background'))
  }
  
  dists <- as.matrix(dist(emb[beads,]))

  group <- sapply(beads,function(bead){ 
    bead.dists <- sort(dists[bead,],decreasing=FALSE)
    neighbors <- names(bead.dists[2:(nCentroidNeighbors+1)])
    centroid <- colMeans(emb[neighbors,])
    self <- emb[bead,]
    if (dist(rbind(centroid,self))[1]>edgeBeadCutoff) { 
      return('edge') 
    } else { 
      return('interior') 
    }})
  UMIcounts <- rbind(UMIcounts,data.frame(nUMI=seu$nCount_RNA[beads],fragment=f,group=group))
}
write.csv(UMIcounts,'../data/nUMI_by_fragment_and_group.csv',row.names=TRUE)
  
# save all beads Seurat object
saveRDS(seu,'../data/seurat_objects/all_pucks_all_beads.rds')

# ================= ROW EMBEDDINGS ===================

# eliminate background beads
seu <- subset(seu,cells=Cells(seu)[seu$nCount_RNA>=40 & !endsWith(as.character(seu$fragment),'-0')])
seu$fragment <- droplevels(seu$fragment)

# fragment angles
angles <- read.csv('../data/tissue_fragment_angles.csv',row.names=1)
angles$rad <- pi*(angles$deg/180)
angles$correction <- (pi/2)-angles$rad

# center at 0,0 and rotate
emb <- Embeddings(seu,reduction='slide')
coords <- NULL
for (f in levels(seu$fragment)) {
  beads <- Cells(seu)[seu$fragment==f]
  theta <- angles[f,'correction']
  curr <- emb[beads,]
  
  centroid <- colMeans(curr)
  x <- curr[,1]-mean(curr[,1])
  y <- curr[,2]-mean(curr[,2])
  new <- data.frame(row.names=beads)
  new$x <- (x*cos(theta))-(y*sin(theta))
  new$y <- (x*sin(theta))+(y*cos(theta))
  new$y <- new$y-max(new$y)
  coords <- rbind(coords,new)
}

# get 'em spaced in nice rows
up <- 1
for (row in rev(names(rows))) {
  right <- .5
  for (f in rows[[row]]) {
    beads <- Cells(seu)[seu$fragment==f]
    coords[beads,'x'] <- coords[beads,'x']+(w*right)
    coords[beads,'y'] <- coords[beads,'y']+(h*up)
    right <- right+1
  }
  up <- up+1
}

# add dim reduc opbject
colnames(coords) <- c('rows_1','rows_2')
seu[['rows']] <- CreateDimReducObject(embeddings=as.matrix(coords[Cells(seu),]),key='rows_',assay='RNA')
DimPlot(seu,reduction='rows',group.by='timepoint')

# SAVE!!!
saveRDS(seu,'../data/seurat_objects/under_tissue.rds')

# ================= PIWI1+ SUBSET ====================

piwiPosBeads <- Cells(seu)[seu[['RNA']]$counts[piwi1,]>0]
seu <- subset(seu,cells=piwiPosBeads)
saveRDS(seu,'../data/seurat_objects/piwi1_positive_subset.rds')
