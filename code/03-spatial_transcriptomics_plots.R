######################################################
# SPATIAL TRANSCRIPTOMICS PLOTS
#
# Generate plots based on slide-seq data:
#   fig1B-D tissue annotation by rows, umap, piwipos tissue, piwipos timepoint
#   figS1C-K bg/interior/edge histogram, slides nUMI >> rows, umaps with clusters, barbell plots
#   FigS2A-J "digital in situ" plots, piwi1 clusters umap and barplots, more piwi+ umaps
#   FigS7A slide-seq pearson correlation box'n'whisker
######################################################

setwd('/n/projects/cb2350/smed_slide-seq/repo/code/')

# ================= IMPORTS ==========================

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggnewscale)

# ================= PARAMS ===========================

dpi <- 300
w <- 10
h <- 10

genes <- c('piwi1'='SMED30007406','vitellogenin'='SMED30003809',
           'pbgd'='SMED30023355','ca7'='SMED30021661',
           'mmp1'='SMED30019930')

# make figures directory, if needed
if (!dir.exists('../figures/')) { dir.create('../figures/') }

rfsGenes <- read.csv('../data/mmp1_FACS-RNA-seq_genes.csv',row.names=1)

# ================= LOAD DATA ========================

# unfiltered single-cell data
all <- readRDS('../data/seurat_objects/all_pucks_all_beads.rds')

# under tissue single-cell data
ut <- readRDS('../data/seurat_objects/under_tissue.rds')

# piwi1+ single-cell data
pp <- readRDS('../data/seurat_objects/piwi1_positive_subset.rds')

# published metadata
meta <- openxlsx::read.xlsx('../data/GSE199348_bead_metadata.xlsx',rowNames=TRUE)

# ================= FIGURE 1 =========================

# fig1b: tissue annotation
p1 <- DimPlot(ut,reduction='rows',group.by='tissue',cols=ut@misc$colors$tissue) +
  coord_fixed() + theme_void() + theme(legend.position='none')
p2 <- DimPlot(ut,reduction='umap',group.by='tissue',cols=ut@misc$colors$tissue) +
  coord_fixed() + theme_void()
p <- wrap_plots(list(p1,p2),ncol=2) + plot_layout(widths=c(2,1))
ggsave('../figures/fig1b_global_tissue_annotation_rows_umap.png',p,width=1.5*w,height=.75*h,units='in',dpi=dpi)

# fig1c: piwipos UMAP by tissue
p <- DimPlot(pp,reduction='umap',group.by='tissue',cols=pp@misc$colors$tissue) +
  coord_fixed() + theme_void()
ggsave('../figures/fig1c_piwi1_positive_tissue_annotation_umap.png',p,height=.6*h,units='in',dpi=dpi)

# fig1d: piwipos UMAP by time
p <- DimPlot(pp,reduction='umap',group.by='timepoint',cols=pp@misc$colors$timepoint) +
  coord_fixed() + theme_void()
# be lazy, let Seurat position cluster labels
p$layers[[2]] <- DimPlot(pp,reduction='umap',group.by='seurat_clusters',label=TRUE,combine=FALSE,raster=FALSE,label.size=5)[[1]]$layers[[2]]
ggsave('../figures/fig1d_piwi1_positive_timepoint_umap.png',p,height=.6*h,units='in',dpi=dpi)

# ================= SUPPLEMENTAL FIGURE 1 ============

# figS1c: bg vs. edge vs. interior log10nUMI histogram
# come back to this

# figS1d: slide log10nUMI all beads, min 40 UMIs
all$log10nUMI <- log10(all$nCount_RNA)
umiRange <- range(all$log10nUMI)
p <- FeaturePlot(all,reduction='slide',features='log10nUMI',raster=FALSE) + 
  scale_color_viridis(direction=-1,name='log10 nUMI',limits=range) +
  coord_fixed() + theme_void() + theme(legend.position='bottom') +
  ggtitle('>= 10 UMIs')
ggsave('../figures/figS1d_all_beads_log10nUMI_slide.png',p,width=w,height=w,units='in',dpi=dpi)

# figS1e: slide log10nUMI min 40 nUMI
min40 <- subset(all,cells=Cells(all)[all$nCount_RNA>=40])
p <- FeaturePlot(min40,reduction='slide',features='log10nUMI',raster=FALSE) + 
  scale_color_viridis(direction=-1,name='log10 nUMI',limits=umiRange) +
  coord_fixed() + theme_void() + theme(legend.position='bottom') +
  ggtitle('>= 40 UMIs')
ggsave('../figures/figS1e_min40_nUMI_log10nUMI_slide.png',p,width=w,height=w,units='in',dpi=dpi)

# figS1f: slide fragment min 40 nUMI
# pretty colos
fragCols <- setNames(colorRampPalette(ut@misc$colors$tissue[1:8])(length(levels(all$fragment))),levels(all$fragment))
fragCols[endsWith(names(fragCols),'-0')] <- '#CCCCCC'
p <- DimPlot(min40,reduction='slide',group.by='fragment',cols=fragCols) +
  coord_fixed() + theme_void() + theme(legend.position='none')
ggsave('../figures/figS1f_min40_nUMI_fragments_slide.png',p,width=w,height=w,units='in',dpi=dpi)

# figS1g: slide timepoint
p <- DimPlot(ut,reduction='slide',group.by='timepoint',cols=ut@misc$colors$timepoint) +
  coord_fixed() + theme_void() + theme(legend.position='bottom')
ggsave('../figures/figS1g_timepoint_slide.png',p,width=w,height=w,units='in',dpi=dpi)

# figS1h: rows timepoint
p <- DimPlot(ut,reduction='rows',group.by='timepoint',cols=ut@misc$colors$timepoint) +
  coord_fixed() + theme_void() + theme(legend.position='bottom')
ggsave('../figures/figS1h_timepoint_rows.png',p,width=w,height=w,units='in',dpi=dpi)

# figS1i: umap global cluster (label clusters)
p <- DimPlot(ut,reduction='umap',group.by='seurat_clusters',label=TRUE,raster=FALSE,combine=FALSE,label.size=5)[[1]] +
  coord_fixed() + theme_void() + theme(legend.position='none')
clusterLabs <- p$layers[[2]] # save for figS1k
ggsave('../figures/figS1i_global_clusters_umap.png',p,width=.7*w,height=.7*w,units='in',dpi=dpi)

# figS1j: log10nUMI, log10nFeatures by global cluster dumbell plots
df <- data.frame(t(sapply(levels(ut$seurat_clusters),function(cl){ 
  summary(log10(ut$nCount_RNA[which(ut$seurat_clusters==cl)])) })))
colnames(df) <- c('min','q1','median','mean','q3','max')
df$cluster <- factor(rownames(df),levels=rev(levels(ut$seurat_clusters)))
df$type <- factor('log10 nUMI',levels=c('log10 nUMI','log10 nFeatures'))
df2 <- data.frame(t(sapply(levels(ut$seurat_clusters),function(cl){ 
  summary(log10(ut$nFeature_RNA[which(ut$seurat_clusters==cl)])) })))
colnames(df2) <- c('min','q1','median','mean','q3','max')
df2$type <- factor('log10 nFeatures',levels=c('log10 nUMI','log10 nFeatures'))
df2$cluster <- factor(rownames(df2),levels=rev(levels(ut$seurat_clusters)))
p <- ggplot(rbind(df,df2),aes(y=cluster,color=type),size=10) + 
  geom_segment(aes(yend=cluster,x=q1,xend=q3),color='#666666',linewidth=.5) +
  geom_point(aes(x=q1)) + geom_point(aes(x=median)) + geom_point(aes(x=q3)) +
  scale_color_manual(values=c('red3','dodgerblue4')) +
  ylab('cluster ID') +
  theme_bw() + theme(legend.position='none',axis.title.x=element_blank(), axis,axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  facet_grid(.~type,switch='x')  
ggsave('../figures/figS1j_nUMI_nFeatures_by_cluster.png',p,width=5,height=5,units='in',dpi=dpi)

# figS1k: umap timepoint (label clusters)
p <- DimPlot(ut,reduction='umap',group.by='timepoint',cols=ut@misc$colors$timepoint,raster=FALSE,combine=FALSE)[[1]] +
  coord_fixed() + theme_void() + theme(legend.position='none')
p$layers[[2]] <- clusterLabs
ggsave('../figures/figS1j_timepoint_umap.png',p,width=.7*w,height=.7*w,units='in',dpi=dpi)


# ================= SUPPLEMENTAL FIGURE 2 ============

# figS2a: digital in situs for piwi1, vitellogelin, pbgd, ca7 (rows)
ptSize <- DimPlot(ut,reduction='rows')[[1]]$plot_env$pt.size
for (gene in c('piwi1','vitellogenin','pbgd','ca7')) {
  smed <- genes[gene]
  df <- data.frame(cbind(Embeddings(ut,reduction='rows'),ut[['RNA']]$counts[smed,]))
  colnames(df) <- c('x','y','expr')
  maxq <- quantile(df$expr,.9995)
  df$expr[df$expr>maxq] <- maxq
  df2 <- data.frame(x=rep(max(df$x),11),y=rep(max(df$y),11),expr=((0:10)/10)*max(df$expr))
  
  p <- ggplot(df,aes(x=x,y=y)) + 
    geom_point(size=ptSize,color='#FFFFFF19') +
    geom_point(size=ptSize,aes(color=expr)) +
    scale_color_gradientn(colors=c('#FF44FF00','#FF44FFFF','#FF99FF'),guide='none') +
    new_scale('color') +
    geom_point(data=df2,aes(x=x,y=y,color=expr),size=2) +
    scale_color_gradientn(name='expression',colors=c('#444444','#FF44FFFF','#FF99FF'),
                          guide=guide_colorbar(frame.colour='white',frame.linewidth=1,ticks.linewidth=1,ticks.colour='white')) +
    theme_void() + coord_fixed(1) +
    geom_point(data=df2,aes(x=x,y=y),color='black',size=4) +
    theme(panel.background=element_rect(fill='black',color='black'),
          plot.background=element_rect(fill="black"),
          plot.title=element_text(colour="white"),
          legend.text=element_text(color='white')) +
    ggtitle(smed)
  ggsave(paste0('../figures/figS2a-',gene,'_',smed,'.png'),p,width=w,height=h,units='in',dpi=dpi)
}

# figS2b: piwi1 feature plot (rows)
p <- FeaturePlot(ut,reduction='rows',features=genes['piwi1']) + 
  theme_void() + coord_fixed()
ggsave('../figures/figS2b_piwi1_SMED30007406_expression_rows.png',p,width=w,height=h,units='in',dpi=dpi)

# figS2c: piwi1+ plot (rows)
ut$piwi1pos <- ut[['RNA']]$counts[genes['piwi1'],]>0
p <- FeaturePlot(ut,reduction='rows',features='piwi1pos',cols=c('#CCCCCC','black'),order=TRUE) + 
  theme_void() + coord_fixed()
ggsave('../figures/figS2c_piwi1_SMED30007406_positive_beads_rows.png',p,width=w,height=h,units='in',dpi=dpi)

# figS2d: piwi1 subclusters (piwi+ umap, label clusters w/ box)
p <- DimPlot(pp,reduction='umap',group.by='seurat_clusters',label=TRUE,label.box=TRUE) +
  coord_fixed() + theme_void() + theme(legend.position='none') +
  ggtitle('piwi1+ subset subclusters')
ggsave('../figures/figS2d_piwipos_subclusters_umap.png',p,width=.6*w,height=.6*h,units='in',dpi=dpi)

# figS2e: barplot of piwi+ cluster by tissue (order per figS2f)
# cluster ordering
lvs <- sapply(levels(pp$seurat_clusters),function(cl){
  x <- table(pp$timepoint[pp$seurat_clusters==cl])[c('06hpa','48hpa')]/table(pp$timepoint)[c('06hpa','48hpa')]
  return(x/sum(x))
})
lvs <- names(sort(apply(lvs,2,function(x){x['06hpa']/x['48hpa']})))
# now tissue
df <- data.frame(sapply(levels(pp$seurat_clusters),function(cl){
  x <- table(pp$tissue[pp$seurat_clusters==cl])/sum(pp$seurat_clusters==cl)
  return(x)
}))
df$tissue <- factor(rownames(df),levels=levels(pp$tissue))
df <- reshape2::melt(df)
df$cluster <- factor(gsub('X','',as.character(df$variable)),levels=lvs)
p <- ggplot(df,aes(x=cluster,y=value,fill=tissue)) + geom_col() +
  scale_fill_manual(values=pp@misc$colors$tissue,name='tissue annotation') +
  theme_bw() + xlab('Piwi1+ cluster ID') + ylab('Proportion of assigned tissue') +
  scale_y_continuous(expand=c(0,0))
ggsave('../figures/figS2e_tissue_by_piwipos_cluster.png',p,width=w,height=.4*w,units='in',dpi=dpi)  

# figS2f: barplot of piwi+ cluster by 6hpa/48hpa
df <- data.frame(sapply(levels(pp$seurat_clusters),function(cl){
  x <- table(pp$timepoint[pp$seurat_clusters==cl])[c('06hpa','48hpa')]/table(pp$timepoint)[c('06hpa','48hpa')]
  return(x/sum(x))
}))
df$timepoint <- factor(rownames(df),levels=c('06hpa','48hpa'))
df <- reshape2::melt(df)
df$cluster <- factor(gsub('X','',as.character(df$variable)),levels=lvs)
p <- ggplot(df,aes(x=cluster,y=value,fill=timepoint)) + geom_col() +
  scale_fill_manual(values=pp@misc$colors$timepoint,name='timepoint') +
  theme_bw() + xlab('Piwi1+ cluster ID') + ylab('Proportion 6hpa vs. 48hpa') +
  scale_y_continuous(expand=c(0,0))
ggsave('../figures/figS2e_6hpa_vs_48hpa_by_piwipos_cluster.png',p,width=w,height=.4*h,units='in',dpi=dpi) 

# figS2g: secretory prediction score (piwi+ umap)
pp$temp <- meta[Cells(pp),'score.Parenchymal']
p <- FeaturePlot(pp,reduction='umap',features='temp',cols=c('#CCCCCC','black')) +
  coord_fixed() + theme_void() + theme(legend.position='none') + 
  ggtitle('prediction.score.Secretory.Cells')
ggsave('../figures/figS2g_secretory_score_piwipos_umap.png',p,width=.6*w,height=.6*h,units='in',dpi=dpi)

# figS2h: mmp1 expression (piwi+ umap)
p <- FeaturePlot(pp,reduction='umap',features=genes['mmp1']) +
  coord_fixed() + theme_void() + theme(legend.position='none')
ggsave('../figures/figS2h_mmp1_SMED30019930_expression_piwipos_umap.png',p,width=.6*w,height=.6*h,units='in',dpi=dpi)

# figS2i: intestine prediction score (piwi+ umap)
pp$temp <- meta[Cells(pp),'score.Intestine']
p <- FeaturePlot(pp,reduction='umap',features='temp',cols=c('#CCCCCC','black')) +
  coord_fixed() + theme_void() + theme(legend.position='none') + 
  ggtitle('prediction.score.Intestine')
ggsave('../figures/figS2i_secretory_score_piwipos_umap.png',p,width=.6*w,height=.6*h,units='in',dpi=dpi)

# figS2j: vitellogenin expression (piwi+ umap)
p <- FeaturePlot(pp,reduction='umap',features=genes['vitellogenin']) +
  coord_fixed() + theme_void() + theme(legend.position='none')
ggsave('../figures/figS2j_vitellogenin_SMED30003809_piwipos_umap.png',p,width=.6*w,height=.6*h,units='in',dpi=dpi)

# ================= SUPPLEMENTAL FIGURE 7 ============

# figS7a: slide-seq pearson correlation box'n'whisker
expr <- ut[['RNA']]$counts[c(rownames(rfsGenes),genes['mmp1']),]

# this will give us rows = fragments, columns = genes
df <- do.call(rbind,lapply(levels(ut$fragment),function(f){
  cells <- Cells(ut)[ut$fragment==f]
  mmp1.expr <- expr[genes['mmp1'],cells]
  frag <- data.frame(t(data.frame(sapply(rownames(rfsGenes),function(g){
    curr <- expr[g,cells]
    pc <- cor.test(x=mmp1.expr,y=curr)
    return(pc$estimate)
  }))))
  colnames(frag) <- gsub('.cor','',colnames(frag))
  rownames(frag) <- f
  return(frag)
}))
df <- reshape2::melt(df)
# order ascending by median
lvs <- rfsGenes[names(sort(sapply(rownames(rfsGenes),function(x){
  curr <- df[df$variable==x,'value']
  median(curr[!is.na(curr)])
}))),'alias']
df$gene <- factor(rfsGenes[as.character(df$variable),'alias'],levels=lvs)
# the styling's not spot-on, but I didn't make the original
p <- ggplot(df,aes(x=gene,y=value)) + 
  geom_jitter(color='#CCCCCC') +
  stat_summary(fun.min=function(x) { quantile(x,.25) },
               fun.max=function(x) { quantile(x,.75) },
               fun=median,
               linewidth=1,pch=3) +
  theme_bw() + theme(panel.grid.major=element_blank(),axis.text.x=element_text(angle=45,hjust=1),
                     axis.title.x=element_blank()) +
  labs(y='Pearson Correlation',title='Correlation in Slide-seq')
ggsave('../figures/figS7a_mmp1_RNA-FAQ-seq_gene_correlation.png',p,width=w,height=.5*h,units='in',dpi=dpi)

