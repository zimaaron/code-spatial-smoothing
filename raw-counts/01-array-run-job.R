## this script is intended to be used to within an array job setting
##
## each job will receive a different integer index and each job should
## use that index to process a different slice of data

## source("~/Dropbox/genetics/a-s-omics/code-spatial-smoothing/raw-counts/01-array-run-job.R")

setwd("~/Dropbox/genetics/a-s-omics/")

##################
## process args ##
##################
## args <- commandArgs(trailingOnly=TRUE)
## date.dir <- as.character( args[1] ) ## for output organization
## time.i   <- as.integer(   args[2] ) ## which embryo timepoint to process (1-8)
## start.i  <- as.integer(   args[3] )  ## slurm only goes up to 3000. add a start index to offset past the limit
## j.i      <- start.i + as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) ## array job specific index
# for testing and building prep objects
## date.dir <- "testing"
## start.u <- 0
## j.i <- 1
## time.i <- 1

##########################################################
## set user selected params, in/out dirs, and load pkgs ##
##########################################################
source("./code-spatial-smoothing/raw-counts/00-load-pkg-set-io-set-params.R")

# pick one
poi.mod <- TRUE#FALSE
zip.mod <- TRUE
zap.mod <- TRUE#FALSE

#######################
## load data, subset ##
#######################
## agg size set in set-params
load(file.path(i.d, glue("pre-fit-obj-binned-{agg.size}.Rdata")))
## subset, make boundary domain
if(!is.null(sub.bound)){
  in.sb <- which(splancs::inout(all.dat[, .(x, y)], sub.bound))
  sub.dat <- all.dat[in.sb, ]
  boundary.coords <- sub.bound
}else{
  sub.dat <- all.dat
  boundary.coords <- NULL #TODO
}
domain <- boundary.coords |> fmesher::fm_segm()

## subset to hierarchy mechanisms
sub.dat <- sub.dat[feat %in% f.to.mod, ]

########################################
## make mesh (params from set-params) ##
########################################

# set mesh params, tuned by agg.size
mesh.params <- data.table(size = c(25, 50, 75, 100),
                          max.edge.inner.mult = c(.1, .15, .35, .5))
max.edge = diff(range(boundary.coords[,1]))/(3*5)
bound.outer = diff(range(boundary.coords[,1]))/3

mesh <- fm_mesh_2d(
  #sub.dat[feat == 'Igf2', .(x, y)],
  boundary = list(fm_segm(c(domain))),
  max.edge = c(mesh.params[size == agg.size,
                           max.edge.inner.mult], 2) * max.edge,
  offset = c(max.edge, bound.outer / 1.5),
  cutoff = max.edge / 10,
  min.angle = 21)
#plot(mesh)
# check which mesh vertices are inside the domain
ds1 <- splancs::as.points(boundary.coords)
#  ds1.poly <- ds1[chull(ds1),]
in.dom <- which(splancs::inout(mesh$loc[, 1:2], ds1))

########################################
## run model across selected features ##
########################################

for(lr.n in f.to.mod[305:length(f.to.mod)]){

  cat('\n\n\n')
  for(i in 1:3){print(glue('ON FEAT: {lr.n}: {which(f.to.mod==lr.n)} of {length(f.to.mod)}\n'))}
  cat('\n\n\n')

  # make data for this run
  run.dat <- sub.dat[feat == lr.n, ]
  set(run.dat, j = "feat.present", value = as.integer(run.dat[, feat.present]))
  set(run.dat, j = "total.present", value = as.integer(run.dat[, total.present]))

  if(FALSE){
    fields::quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
                       nx = run.dat[, length(unique(x))], ny = run.dat[, length(unique(y))])
  }

  # run the model
  source(file.path(c.d, "02-array-run-model.R"))
}
