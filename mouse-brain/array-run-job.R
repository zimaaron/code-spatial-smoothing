## this script is intended to be used to within an array job setting
##
## each job will receive a different integer index and each job should
## use that index to process a different slice of data

##################
## process args ##
##################
args <- commandArgs(trailingOnly=TRUE)

## date.dir <- as.character(args[1]) ## for output organization
## start.i <- as.integer(args[2])   ## slurm only goes up to 1001. add a start index to offset past the limit
## j.i     <- start.i + as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) ## array job specific index

date.dir <- "for-yale-pres-plots" ## for output organization
start.i <- 0   ## slurm only goes up to 1001. add a start index to offset past the limit
j.i     <- 1


###############
## libraries ##
###############
require(glue)
require(Seurat)


##############
## setup io ##
##############

setwd("~/a-s-omics")
setwd("~/Dropbox/genetics/a-s-omics")

i.d <- "./data-inputs/mouse-brain"
# i.d <- "~/Dropbox/genetics/a-s-omics/data-inputs"
o.d <- file.path("./data-outputs/mouse-brain", date.dir)
c.d <- "./code-spatial-smoothing" ## other code

dir.create(o.d, recursive = T, showWarnings = FALSE)

## TODO write readme that saves with each run
## or append the readme to one (table?) file? this seems like a better idea

## specify the names of lr combos that we'll save draw outputs for.
lr.to.save.draws <- c("Hgf—Met", "Shh—Boc", "Ndp—Fzd4",
                      "Fgf1—Fgfr1", "Bmp2—Eng", "Gdnf—Ret",
                      "Lpl—Lrp1", "Ndp—Fzd4", "Agrn—Musk", "Apoe—Ldlr")

###############
## load data ##
###############


## if(j.i <= 2251){ ## do left side


##   side <- "A" # anterior

##   load(file.path(i.d, "brain.NICHES.Robj")) # loads brain obj

##   # find data at:
##   # brain@assays$NeighborhoodToCell_ALRA@data

##   lr.n <- brain@assays$NeighborhoodToCell_ALRA@data@Dimnames[[1]][j.i] ## name
##   lr.d <- brain@assays$NeighborhoodToCell_ALRA@data[eval(lr.n), ] |> as.numeric() ## data

##   xy.obs <- cbind(brain@images$anterior1@coordinates$col,
##                   -brain@images$anterior1@coordinates$row) ## TODO why are rows/cols labeled like this?
## }else if(j.i > 2251){

##   j.i <- j.i - 2251

##   side <- "P" # posterior

##   load(file.path(i.d, "brain2.NICHES.Robj")) # loads brain2 obj

##   # find data at:
##   # brain@assays$NeighborhoodToCell_ALRA@data

##   lr.n <- brain2@assays$NeighborhoodToCell_ALRA@data@Dimnames[[1]][j.i] ## name
##   lr.d <- brain2@assays$NeighborhoodToCell_ALRA@data[eval(lr.n), ] |> as.numeric() ## data

##   xy.obs <- cbind(brain2@images$posterior1@coordinates$col,
##                   -brain2@images$posterior1@coordinates$row) ## TODO why are rows/cols labeled like this?
## }

side <- "A" # anterior
load(file.path(i.d, "brain.NICHES.Robj")) # loads brain obj
# find data at:
# brain@assays$NeighborhoodToCell_ALRA@data


for(lr.n in sort(lr.to.save.draws[c(5, 6, 7, 10)])){

  lr.d <- brain@assays$NeighborhoodToCell_ALRA@data[eval(lr.n), ] |> as.numeric() ## data

  xy.obs <- cbind(brain@images$anterior1@coordinates$col,
                  -brain@images$anterior1@coordinates$row) ## TODO why are rows/cols labeled like this?




colnames(xy.obs) <- c("x", "y")
if(j.i == 1){
  write.csv(xy.obs, file.path(o.d, glue("xy-obs-locs-{side}.csv")), row.names = F)
}

###############
## run model ##
###############

  source(file.path(c.d, "mouse-brain/array-run-model.R"))

}


###############
## save outs ##
###############

## done, for now, in array-run-model.R
