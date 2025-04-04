# This script pulls public spatial data and extracts the ligand and receptor information

# SetWD
setwd("/Users/msbr/GitHub/a-s-omics/code-spatial-smoothing")

# Load packages
require(Seurat)
require(SeuratData)
require(ggplot2)

# Load local functions
source('LRExtract.R')

# Load data
brain1 <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "posterior1")

# Normalize data (here we will use log normalization)
brain1 <- NormalizeData(brain1, assay = "Spatial", verbose = FALSE)
brain2 <- NormalizeData(brain2, assay = "Spatial", verbose = FALSE)

# Merge into one object
brain.merge <- merge(brain1, brain2)

# Scale data (all together)
brain.merge <- ScaleData(brain.merge)

# Test out plotting
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))
SpatialFeaturePlot(brain.merge, features = c("Fgf1", "Fgfr2"))

# Format the x and y coordinates for downstream
brain1@meta.data$x <- brain1@images$anterior1@coordinates$row
brain1@meta.data$y <- brain1@images$anterior1@coordinates$col
brain2@meta.data$x <- brain2@images$posterior1@coordinates$row
brain2@meta.data$y <- brain2@images$posterior1@coordinates$col

# Extract regions
LR.output <- list()
LR.output$brain1 <- LRExtract(brain1)
LR.output$brain2 <- LRExtract(brain2)
View(LR.output)

# Save for later
save(LR.output,file = 'LR.output.2023-11-08.Robj')
