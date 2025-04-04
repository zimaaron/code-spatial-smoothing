## this script is intended to be used to within an array job setting
##
## each job will receive a different integer index and each job should
## use that index to process a different slice of data



##################
## process args ##
##################
args <- commandArgs(trailingOnly=TRUE)

date.dir <- as.character( args[1] ) ## for output organization
time.i   <- as.integer(   args[2] ) ## which embryo timepoint to process (1-8)
start.i  <- as.integer(   args[3] )  ## slurm only goes up to 3000. add a start index to offset past the limit
j.i      <- start.i + as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) ## array job specific index

# for profile timing
#time.mat <- matrix(ncol = 4, nrow = 3)
#for(time.i in 1:4){
#  rm(list = setdiff(ls(), c("time.i", "time.mat", "tic")))

# start the timer
tic <- function () {
    now <- proc.time()
    function () {
        proc.time() - now
    }
}
toc <- tic()

# for testing and building prep objects
## date.dir <- "testing"
## start.u <- 0
## j.i <- 1
## time.i <- 1
## setwd("~/Dropbox/genetics/a-s-omics/")


###############
## libraries ##
###############
require(glue)
require(Seurat)

##############
## setup io ##
##############

## setwd("~/a-s-omics")

i.d <- "./data-inputs/mouse-embryo"
# i.d <- "~/Dropbox/genetics/a-s-omics/data-inputs"
o.d <- file.path("./data-outputs/mouse-embryo", date.dir)
c.d <- "./code-spatial-smoothing" ## other code

dir.create(o.d, recursive = T, showWarnings = FALSE)

## TODO write readme that saves with each run
## or append the readme to one (table?) file? this seems like a better idea

## specify the names of lr combos that we'll save draw outputs for.
lr.to.save.draws <- c(
"Vegfa—Kdr","Angpt1—Tie1",
"Bmp2—Bmpr1","Ccl2—Ccr5",
"Col4a1—Itga1","Cthrc1—Fzd6",
"Efna4—Epha5","Efna4—Epha6",
"Efna4—Epha7","Egf—Erbb3")

# TODO is 5% too high? it drops about half of the data
drop.pct <- 5 # drop LR features that have less than drop.pct% non-zero cts

###############
## load data ##
###############

## data organization
# niches.alra.list is a list with one timepoint per entry
load(file.path(i.d, "niches.alra.list.Robj")) # loads data processed by Sam ## TODO should this be the norm object?

# proceed to preprocess objects for each/all time points here
# using prep-data-inputs.R
# this only needs to be done once

## Make a list of mechanisms to filter out
## Filter out mechanisms that:
# * don't appear more than 5% of time across all timepoints IN THE RAW COUNTS

total.spots <- 0
non.zero.ct <- rep(0, dim(niches.alra.list[[1]]$NeighborhoodToCell@assays$NeighborhoodToCell@layers$counts)[1])
non.zero.ct.mat <- matrix(0, ncol = length(niches.alra.list),
                          nrow = dim(niches.alra.list[[1]]$NeighborhoodToCell@assays$NeighborhoodToCell@layers$counts)[1])
for(i in 1:length(niches.alra.list)){

  time.spots <- dim(niches.alra.list[[i]]$NeighborhoodToCell@assays$NeighborhoodToCell@layers$counts)[2]

  total.spots <- total.spots + time.spots

  time.non.zero.ct <- rowSums(niches.alra.list[[i]]$NeighborhoodToCell@assays$NeighborhoodToCell@layers$counts)

  non.zero.ct <- non.zero.ct + time.non.zero.ct

  non.zero.ct.mat[, i] <- time.non.zero.ct / time.spots
}
non.zero.ct.mat <- cbind(non.zero.ct.mat, non.zero.ct / total.spots)

## how many will a 1% thresh drop? 5%?
# sum(non.zero.ct / total.spots < .05)
## 417 out of 2036 dropped at 5% (101 dropped at 1% )

# all lr feature names
all.lr.n <- niches.alra.list[[1]]$NeighborhoodToCell@assays$NeighborhoodToCell@features |> rownames()

# index of those we will use
keep.lr.idx <- which(non.zero.ct / total.spots >= .05)
# run 1619 LR features at 5% thresh

# select one of the "keepers" based on job index (they are constant across time)
lr.n <-  all.lr.n[keep.lr.idx][j.i]

print(glue("processing feature {lr.n} at time {time.i}"))

# grab single vector of data for the selected time  LR feature
lr.d <- niches.alra.list[[time.i]]$NeighborhoodToCell@assays$NeighborhoodToCell@layers$data[which(all.lr.n == lr.n), ] |> as.numeric()

xy.obs <- cbind( niches.alra.list[[time.i]]$NeighborhoodToCell@meta.data$x |> as.numeric(),
                 niches.alra.list[[time.i]]$NeighborhoodToCell@meta.data$y |> as.numeric())


#}else if(j.i > 2251){

  ## j.i <- j.i - 2251

  ## side <- "P" # posterior

  ## load(file.path(i.d, "brain2.NICHES.Robj")) # loads brain2 obj

  ## # find data at:
  ## # brain@assays$NeighborhoodToCell_ALRA@data

  ## lr.n <- brain2@assays$NeighborhoodToCell_ALRA@data@Dimnames[[1]][j.i] ## name
  ## lr.d <- brain2@assays$NeighborhoodToCell_ALRA@data[eval(lr.n), ] |> as.numeric() ## data

  ## xy.obs <- cbind(brain2@images$posterior1@coordinates$col,
  ##                 -brain2@images$posterior1@coordinates$row) ## TODO why are rows/cols labeled like this?
#}

colnames(xy.obs) <- c("x", "y")
if(j.i == 1){
  write.csv(xy.obs,
            #file.path(o.d, glue("xy-obs-locs-time-{sprintf('%02d',time.i)}.csv")),
            file.path(o.d, glue("xy-obs-locs-time-{time.i}.csv")),
            row.names = F)
}

rm(niches.alra.list)
for(i in 1:3){gc()}

###############
## run model ##
###############

source(file.path(c.d, "mouse-embryo/array-run-model.R"))

# for profile timing
# }

###############
## save outs ##
###############

## done, for now, in array-run-model.R
