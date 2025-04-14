############
## set IO ##
############

# setwd("~/Dropbox/genetics/a-s-omics/")
i.d <- file.path(getwd(), "/data-inputs/mouse-embryo-raw")
o.d <- file.path(getwd(), "/data-outputs/mouse-embryo-raw", Sys.Date())
c.d <- file.path(getwd(), "/code-spatial-smoothing/raw-counts")
dir.create(o.d, recursive = T, showWarnings = FALSE)


##########
## pkgs ##
##########
library(ggplot2)
library(INLA)
## INLA:::inla.binary.install("Rocky Linux-8") ## seems like calling this once per INLA install is ok
## to fix the following error: /home/aeoz001/R/x86_64-pc-linux-gnu-library/4.2/INLA/bin/linux/64bit/inla.mkl: /lib64/libm.so.6: version `GLIBC_2.29' not found (required by /home/aeoz001/R/x86_64-pc-linux-gnu-library/4.2/INLA/bin/linux/64bit/libicuuc.so.66)
## TODO ask bisonnet if this can be fixed on the cluster
library(splancs)
library(fmesher)
library(fields)
require(purrr)
require(Seurat)
require(SeuratObject)
require(glue)
require(data.table)
require(dplyr)
require(inlabru)
require(ggplot2)
source(file.path(c.d, '../make.pred.grid.v2.R'))


################
## set params ##
################
# in order of appearance in scripts

## ---------
## SCRIPT 00
## ---------

# bin size of raw data
agg.size <- 50

# subset data to some region
sub.bound <- NULL
sub.bound <- splancs::as.points(matrix(c(15000, 25000, 25000, 15000, 15000,
                                         11500, 11500, 18000, 18000, 11500), ncol = 2))

# features to model

## OPTION 1: hand-coded selection
## ## initial selection chosen by aoz in dec 2024
## l.to.mod <- c("Mdk", "Igf2", "Dlk1", "Col1a2", "Col1a1", "Fn1",
##               "Vcan", "Col14a1", "Slit2", "Rtn4", "Slit3", "Sema3a",
##               "Cxcl12")
## r.to.mod <- c("Itgb1", "Igf2r", "Igf1r", "Lrp1", "Notch2", "Cd44",
##               "Robo2", "Tnfrsf19", "Robo1", "Nrp1")
## m.to.mod <- c("MDK_ITGB1","IGF2_IGF2R", "IGF2_IGF1R", "MDK_LRP1",
##               "DLK1_NOTCH2", "COL1A2_CD44", "COL1A1_CD44", "FN1_CD44",
##               "VCAN_CD44", "COL14A1_CD44", "SLIT2_ROBO2",
##               "RTN4_TNFRSF19", "SLIT2_ROBO1", "SLIT3_ROBO2",
##               "SEMA3A_NRP1", "CXCL12_ITGB1")
## specific.l.r.to.mod <- c('Igf2', 'Igf2r', 'Igf1r', 'Mdk', 'Itgb1', 'Lrp1')



## OPTION 2:
## large (but tailored) list of mechanisms msr selected in jan 2025
## m.to.mod <- list.files("~/Dropbox/a-s-omics/mechanisms of high interest") |> strsplit(, split = " ") |> lapply(FUN = function(x){x[1]}) |> unlist() |> gsub(pattern = "—", replacement = "-")
## l.to.mod <- strsplit(m.to.mod, "-") |> lapply(FUN = function(x){x[1]}) |> unlist() |> unique()
## r.to.mod <- strsplit(m.to.mod, "-") |> lapply(FUN = function(x){x[2]}) |> unlist() |> unique()

## OPTION 3:
## list of mechanisms that relate to ligands of interest
# l.of.interest <- c('Wnt4', 'Wnt5a', 'Bmp4', 'Il17', 'Il17b','Angptl4','Igf1')
## ligands for grant, grab all receptors from hierarchy, subset to those available in data
# load(file.path(i.d, "fantom-hierarchy.Robj"))
# hier <- as.data.table(output)
# hier.s <- subset(hier, LIGAND %in% toupper(l.of.interest))

## load data from sam
load(file.path(c.d, "features.to.use.2025-04-04.Robj"))
f.to.mod <- to.use$gene

f.to.drop <- c("a", "Calcb", "Cckar") # TODO, look into failure. better priors? also should add trycatch to both fit and predict
f.to.mod <- setdiff(f.to.mod, f.to.drop)


# prediction grid resolution (same in x and y)
# TODO, different x and y resolutions
# Not currently used
pred.grid.res <- 200
