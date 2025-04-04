# a-s-omics 2024-05-22
setwd("/Users/msbr/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/a-s-omics")
setwd("~/Dropbox/genetics/a-s-omics")
require(glue)
require(Seurat)
require(ggplot2)
require(cowplot)

#####
# aarons  code here
## combine all the mean estimates from an output dir

o.d <- file.path(getwd(), "data-outputs/2024-04-17/")

a.f.p <- list.files(o.d, full.names = T)
a.f   <- list.files(o.d)

# subset to mean estimates
m.f.p <- grep("P-resp-est.csv", a.f.p, value = T)
m.f   <- grep("P-resp-est.csv", a.f, value = T)

# make a big storage matrix, grabbing dims from one file
d.t <- read.csv(m.f.p[1])

a.r <- data.frame(matrix(nrow = nrow(d.t), ncol = length(m.f.p)))

# loop through files, read in, combine
for(i in 1:length(m.f.p)){
  if(i %% 250 == 0){ message(glue("on iter {i}"))}
  d.t <- read.csv(m.f.p[i])
  a.r[, i] <- d.t[, 1]
  colnames(a.r)[i] <- gsub("-est.csv", "", m.f[i])

  ## TODO combine the DIC results
  # glue("{lr.n}-DIC-comp.csv")

  ## TODO combine the means of preds/obs results
  # glue("{lr.n}obs-est-means-comp.csv")
}

# save output
write.csv(a.r, file.path(o.d, "P-all-smooth-ests.csv"), row.names = F)

#####
# step 1: turn into a seurat obj
test <- CreateAssay5Object(counts = t(a.r),
                           data = t(a.r)) # put the same data in both slots
test <- CreateSeuratObject(test,assay = 'BayesianNICHES') # label
test # confirm
#31572 spots

#fix rownames
temp <- unlist(strsplit(rownames(test),split = '-'))
temp <- temp[temp!='P']
temp <- temp[temp!='resp']
rownames(test) <- temp # aaron syas use gsub for this

# add xy coords
xy <- read.csv("~/Large Files/a-s-omics_transfers_temp/2024-04-17/xy-pred-locs-P.csv")
rownames(xy) <- colnames(test)
test <- AddMetaData(test,metadata = xy)

# Add image data for spatial plotting
xy.fixed <- xy
xy.fixed <- xy[,c('y','x')]
xy.fixed$y <- -xy.fixed$y

test@images$image =  new(
  Class = 'SlideSeq',
  assay = "BayesianNICHES",
  key = "image_",
  coordinates = xy.fixed)

#step 2: load original niches object
load("~/Library/CloudStorage/GoogleDrive-michasam.raredon@yale.edu/My Drive/Large_Files_Public/brain2/brain2.NICHES.Robj")

# step 3: limit bayesian niches to only those rows which havent deviated in range beyond a threshold

temp <- test@assays$BayesianNICHES@layers$data
temp <- as.matrix(temp)
temp[temp>100000] <- 100000 #set max
temp[temp<-100000] <- -100000 #set min
new.maxes  <- qlcMatrix::rowMax(temp)# slow
new.mins  <- qlcMatrix::rowMin(temp)# slow

orig.maxes  <- qlcMatrix::rowMax(brain2@assays$NeighborhoodToCell_ALRA@data[rownames(test@assays$BayesianNICHES),])
orig.mins  <- qlcMatrix::rowMin(brain2@assays$NeighborhoodToCell_ALRA@data[rownames(test@assays$BayesianNICHES),])

thresh.factor  <- 0.1

range.maxes <- orig.maxes*(1+thresh.factor)
range.mins <- orig.mins*(1-thresh.factor)


filter.table <- data.frame(orig.maxes = as.numeric(orig.maxes),
                           orig.mins = as.numeric(orig.mins),
                           new.maxes = as.numeric(new.maxes),
                           new.mins = as.numeric(new.mins),
                           range.maxes = as.numeric(range.maxes),
                           range.mins = as.numeric(range.mins))
rownames(filter.table) <- rownames(test)
filter.table$in.bottom <- filter.table$new.mins > filter.table$range.mins
filter.table$in.top <- filter.table$new.maxes < filter.table$range.maxes
filter.table$positive <- filter.table$new.mins > 0

keep <- rownames(filter.table[filter.table$in.top == 'TRUE' &
                                filter.table$positive == 'TRUE',])

test <- test[keep,]

# step 4: standard single cell analysis pipeline

# scale data by row
test <- ScaleData(test)

# find variable features
test <- FindVariableFeatures(test)

# look at metadata and variable features
test@assays$BayesianNICHES@meta.data

# run PC
test <- RunPCA(test,npcs = 100,features = rownames(test)) # using all features?? why not, only 1550...
ElbowPlot(test,ndims = 100)
DimHeatmap(test,dims = 1:9,balanced = T,cells = 300)
DimHeatmap(test,dims = 10:18,balanced = T,cells = 300)
DimHeatmap(test,dims = 19:27,balanced = T,cells = 300)
DimHeatmap(test,dims = 28:36,balanced = T,cells = 300)
DimHeatmap(test,dims = 37:45,balanced = T,cells = 300)
DimHeatmap(test,dims = 46:54,balanced = T,cells = 300)
DimHeatmap(test,dims = 55:63,balanced = T,cells = 300)
DimHeatmap(test,dims = 64:72,balanced = T,cells = 300)
DimHeatmap(test,dims = 73:81,balanced = T,cells = 300)
DimHeatmap(test,dims = 82:90,balanced = T,cells = 300)
DimHeatmap(test,dims = 91:99,balanced = T,cells = 300)

# embed
test <- RunUMAP(test,dims = 1:80)
test <- RunTSNE(test,dims = 1:80)
Seurat::DimPlot(test,reduction = 'umap')+ NoLegend()
Seurat::DimPlot(test,reduction = 'tsne')+ NoLegend()

# cluster
test <- FindNeighbors(test,dims = 1:80)
test <- FindClusters(test,resolution = 2)

# quick plots to check cluster resolution
Seurat::DimPlot(test,reduction = 'umap')+ NoLegend()
Seurat::DimPlot(test,reduction = 'tsne')+ NoLegend()

# determine available x and y coordinate values
table(test$x)
table(test$y)

# place a line on the spatial plot to see what clusters are intersected
y.coordinate <- -39.4974874371859
ggplot(data=test@meta.data,
       aes(x=x,y=y,color=seurat_clusters))+geom_point()+theme_minimal()+
  geom_hline(yintercept =y.coordinate)+NoLegend()

# define barcodes for this line
barcodes <- which(test@meta.data$y==y.coordinate)

# see how it looks
p0 <- SpatialDimPlot(test,label  = T,label.size = 2)+NoLegend()
p0.1 <- ggplot(data=test@meta.data,
               aes(x=x,y=y,color=seurat_clusters))+geom_point()+theme_minimal()+
  geom_hline(yintercept =y.coordinate)+NoLegend()
p1 <- Seurat::DimPlot(test,reduction = 'umap',label = T)+ NoLegend()
p2 <- Seurat::DimPlot(test,reduction = 'tsne',label = T)+ NoLegend()
p3 <- Seurat::DimPlot(test,reduction = 'umap',cells.highlight = barcodes,label = T)+ NoLegend()
p4 <- Seurat::DimPlot(test,reduction = 'tsne',cells.highlight = barcodes,label = T)+ NoLegend()
cowplot::plot_grid(p0,p1,p2,p0.1,p3,p4,nrow = 2)

# determine x locations for line divisions
meta <- test@meta.data[barcodes,]
meta$division <- NA
meta$division.point <- NA
for(i in 2:nrow(meta)){
  if(meta[i,]$seurat_clusters != meta[i-1,]$seurat_clusters){
  meta[i,]$division  <- 'division_point'
  meta[i,]$division.point <- sum(meta[i,]$x,meta[i-1,]$x)/2 }else{}
}
division.locations <- meta[!is.na(meta$division),]$division.point
#& idea could loop and get cluster boundaries in space?

# define clusters
clusters <- unique(meta$seurat_clusters)
clusters
# find some markers
mark <- FindAllMarkers(test,only.pos = T,slot = 'scale.data',logfc.threshold = 1)
mark$ratio <- mark$pct.1/mark$pct.2
mark$power <- mark$ratio*mark$avg_log2FC

# plot a demo
moi <- 'Kitl—Kit'
moi <- 'Fgf1—Fgfr1'
moi <- 'Csf3—Csf1r'
moi <- 'Wnt3—Fzd5'
moi <- 'Wnt3—Fzd8'
moi <- 'Igf1—Igf1r'
moi <- 'Fgf7—Fgfr1'
moi <- 'Shh—Boc'
moi <- 'Pdgfc—Flt4'
moi <- 'Ntn1—Neo1'
moi <- 'Cgn—Tgfbr1'
moi <- 'Tgfb2—Tgfbr1'
moi <- 'Fgf10—Fgfr3'
p1 <- SpatialDimPlot(test,label =T,label.size = 1.5)  + NoLegend()
p2 <- SpatialFeaturePlot(test,moi,interactive = F)
cowplot::plot_grid(p1,p2,ncol=2)

# Experiment 2024-05-16 line diagrams
# features <- c('Nlgn2—Nrxn3','Apoe—Ldlr','Bdnf—Ntrk2','Bmp5—Acvr1',
#               'Bmp6—Acvr2a','Calm1—Trpc3','Calm2—Pde1a','Lamb1—Itgb1',
#               'Kitl—Kit','Fgf1—Fgfr1','Hbegf—Egfr','Igf1—Igf1r')

features <- c('Csf3—Csf1r','Wnt3—Fzd5','Wnt3—Fzd8','Igf1—Igf1r','Kitl—Kit','Fgf7—Fgfr1','Shh—Boc','Pdgfc—Flt4','Ntn1—Neo1','Cgn—Tgfbr1')
colors <- c('#FFCAB1','#55DDE0','#33658A','#695958','#540D6E',
            '#F72585','#05299E','#000004','#FFC857','#2A9134','#CC2936','magenta')
names(colors) <- gsub('-','.',features)
names(colors) <- gsub('—','.',names(colors))
table(test@meta.data$y)
table(test@meta.data$x)
data.to.plot <- t(data.frame(test@assays$BayesianNICHES$scale.data[features,barcodes]))
data.to.plot <- data.frame(data.to.plot)

cluster.label <- test@meta.data[barcodes,]$seurat_clusters
data.to.plot$x <- meta$x
data.to.plot$cluster <-  cluster.label
data.to.plot2 <- reshape2::melt(data.to.plot,c('x','cluster'))
ggplot(data.to.plot2,
       aes(x=x,y=value,group=variable,color=variable))+
  geom_point()+theme_bw()+
  scale_color_manual(values = colors)

# add lines separating clusters
table(droplevels(cluster.label))

p2 <- ggplot(data.to.plot2,
             aes(x=x,y=value,group=variable,color=cluster))+geom_point()+
  theme_bw()+
  geom_vline(xintercept =division.locations,color = "black", size=0.5)

p1 <- ggplot(data.to.plot2,
             aes(x=x,y=value,group=variable,color=variable))+geom_point()+geom_line()+
  theme_bw()+
  geom_vline(xintercept =division.locations,color = "black", size=0.5)+
  scale_color_manual(values = colors)
cowplot::plot_grid(p1,p2,nrow=2,align = T)


# look at feature graph?
sub <- subset(test,cells = WhichCells(test,idents='3'))
sub <- ScaleData(sub)
sub <- FindVariableFeatures(sub)
FOI <- VariableFeatures(sub)
data.to.crush <- sub@assays$BayesianNICHES$data[FOI,]
crushed.data <- rowMeans(data.to.crush)
crushed.data <- crushed.data[crushed.data>0 & crushed.data<1]
feature.dist <- dist(crushed.data,method = 'manhattan') #distance matrix between feature rows
feature.dist.log <- log(feature.dist) # this helps make the distances all about the same size, there are a few with enormous values (highly different pairs of cells)

fit <- cmdscale(feature.dist, k=2) # MDS
to.plot <- data.frame(fit)
labels <- names(crushed.data)
to.plot$feature <- labels
ggplot(to.plot,
       aes(x=X1,
           y=X2,
           #color=labels,
           label=labels))+
  geom_point()+
  ggrepel::geom_text_repel(max.overlaps = 25)+
  theme_minimal()#+scale_x_log10()+scale_y_log10()





SpatialFeaturePlot(test,'Ctgf—Lrp6')

p1 <- SpatialDimPlot(test,label =T,label.size = 3)+ NoLegend()
p2 <- DimPlot(test,label=T,reduction = 'umap')+ NoLegend()
p3 <- DimPlot(test,label=T,reduction = 'tsne')+ NoLegend()
plot_grid(p1,p2,p3,ncol = 3)

# ~~~ legacy code below~~~
# Really weird experiment - what happens if we merge both the smoothed and raw data and look at together?
# load up starting data that I gave to Aaron
brain$Measure <- 'Raw'
test$Measure <- 'Smoothed'

# assuming Aaron is using the alra data - check
brain.data <- CreateAssay5Object(counts = brain@assays$NeighborhoodToCell_ALRA@data,
                                 data = brain@assays$NeighborhoodToCell_ALRA@data)
test.data <- CreateAssay5Object(counts = test@assays$BayesianNICHES$data,
                                data = test@assays$BayesianNICHES$data)
test.data

# Make Seurat Object with a common assay tag
brain.obj <- CreateSeuratObject(brain.data,assay = 'CompiledAssay')
test.obj <- CreateSeuratObject(test.data,assay = 'CompiledAssay')

# Fix rownames before going further
rownames(test.obj)
temp <- strsplit(rownames(test.obj),split='-')
ReformatStrSplit <- function(strsplit.output){
  temp <- strsplit.output
  grab <- c()
  for(i in 1:length(temp)){
    grab[i] <- temp[[i]][[1]]
  }
  return(grab)
}
temp <- ReformatStrSplit(temp)
temp
rownames(test.obj) <- temp
rownames(test.obj)

sum(!(rownames(test.obj) %in% rownames(brain.obj)))
sum(!(rownames(brain.obj) %in% rownames(test.obj))) # Because only working with half the data

# Add all metadata back
test.obj <- AddMetaData(test.obj,metadata = test@meta.data)
brain.obj <- AddMetaData(brain.obj,metadata = brain@meta.data)

# now we can merge
merge <- merge(brain.obj,test.obj)

table(merge$Measure)

# scale
merge <- JoinLayers(merge)
merge <- ScaleData(merge)
# variable features
merge <- FindVariableFeatures(merge)
# PCA
merge <- RunPCA(merge,npcs = 100)
ElbowPlot(merge,ndims = 100)
DimHeatmap(merge,dims = 1:9,balanced = T,cells = 300)
DimHeatmap(merge,dims = 10:18,balanced = T,cells = 300)
DimHeatmap(merge,dims = 19:27,balanced = T,cells = 300)
DimHeatmap(merge,dims = 28:36,balanced = T,cells = 300)
DimHeatmap(merge,dims = 37:45,balanced = T,cells = 300)
DimHeatmap(merge,dims = 46:54,balanced = T,cells = 300)

merge <- RunUMAP(merge,dims = 1:25)

merge <- FindNeighbors(merge,dims = 1:25)

merge <- FindClusters(merge,resolution = 1)

# experiment
DimPlot(merge,label=T)
DimPlot(merge,label=T,group.by = 'Measure')
DimPlot(merge,label=T,group.by = 'seurat_clusters')
