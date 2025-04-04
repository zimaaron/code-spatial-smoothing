
## this script is intended to be used to within an array job setting
##
## each job will receive a different integer index and each job should
## use that index to process a different slice of data
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

run.lr.joint <- TRUE
run.lr.indep <- !run.lr.joint

run.lr.joint.nonlinear <- TRUE
if(run.lr.joint.nonlinear){
  run.lr.joint <- F
  run.lr.indep <- F
}

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
sub.dat[, feat := toupper(feat)]
keep.l <- hier.s[, .(LIGAND)] |> unlist()
keep.l <- which(keep.l %in% sub.dat[, feat])
keep.r <- hier.s[, .(RECEPTOR)] |> unlist()
keep.r <- which(keep.r %in% sub.dat[, feat])
keep.idx <- intersect(keep.l, keep.r) # keep mechanisms with both lig and rec in dataset
hier.s <- hier.s[keep.idx, ]

m.to.mod <- hier.s[, MECHANISM]
m.to.mod <- gsub("_", "-", m.to.mod)

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
if(run.lr.indep){
  for(lr.n in c(l.to.mod, r.to.mod)){
    #for(lr.n in c('Mdk', 'Itgb1', 'Igf2', 'Igf2r', 'Igf1r', 'Lrp1')){

    cat('\n\n\n')
    for(i in 1:3){print(glue('ON FEAT: {lr.n}: {which(c(l.to.mod, r.to.mod) == lr.n)} of {length(c(l.to.mod, r.to.mod))}\n'))}
    cat('\n\n\n')

    # make data for this run
    run.dat <- sub.dat[feat == lr.n, ]
    set(run.dat, j = "feat.present", value = as.integer(run.dat[, feat.present]))
    set(run.dat, j = "total.present", value = as.integer(run.dat[, total.present]))

    ## # plot feature and total
    ## par(mfrow = c(1, 2))
    ## fields::quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
    ##                    main = glue('{lr.n}'))
    ## fields::quilt.plot(run.dat[, x], run.dat[, y], run.dat[, total.count],
    ##                    main = 'Total')
    ## par(mfrow = c(1, 1))

    # run the model
    source(file.path(c.d, "02-array-run-model.R"))
  }
}



if(run.lr.joint){

  ## check if the mechanism has data for both lig and recep, otherwise
  ## remove it from the list
  bad.m <- as.character();
  for(lr.n in m.to.mod){;
    l.n <- strsplit(lr.n, "-")[[1]][1];
    r.n <- strsplit(lr.n, "-")[[1]][2];
    l.run.dat.n <- sub.dat[feat == l.n, .N];
    r.run.dat.n <- sub.dat[feat == r.n, .N];
    if(l.run.dat.n == 0 | r.run.dat.n == 0){;
      bad.m <- c(bad.m, lr.n);
    };
  };
  m.to.mod <- setdiff(m.to.mod, bad.m);


  # subset to relevant data and run
  for(lr.n in m.to.mod){


    cat('\n\n\n');
    for(i in 1:3){print(glue('ON FEAT: {lr.n}: {which(m.to.mod == lr.n)} of {length(m.to.mod)}\n'))};
    cat('\n\n\n');

    l.n <- strsplit(lr.n, "-")[[1]][1];
    r.n <- strsplit(lr.n, "-")[[1]][2];

    l.run.dat <- sub.dat[feat == l.n, ];
    set(l.run.dat, j = "feat.present", value = as.integer(l.run.dat[, feat.present]));
    set(l.run.dat, j = "total.present", value = as.integer(l.run.dat[, total.present]));
    setnames(l.run.dat,
             c("feat", "feat.count", "feat.observed", "feat.present"),
             c("feat.l", "feat.l.count", "feat.l.observed", "feat.l.present"));
    r.run.dat <- sub.dat[feat == r.n, .(x, y, feat, feat.count, feat.observed, feat.present)];
    set(r.run.dat, j = "feat.present", value = as.integer(r.run.dat[, feat.present]));
    setnames(r.run.dat,
             c("feat", "feat.count", "feat.observed", "feat.present"),
             c("feat.r", "feat.r.count", "feat.r.observed", "feat.r.present"));
    run.dat <- merge(l.run.dat, r.run.dat, by = c("x", "y"));

    source(file.path(c.d, "02-array-run-model.R"))
  }
}

if(run.lr.joint.nonlinear){

  ## check if the mechanism has data for both lig and recep, otherwise
  ## remove it from the list
  bad.m <- as.character();
  for(lr.n in m.to.mod){;
    l.n <- strsplit(lr.n, "-")[[1]][1];
    r.n <- strsplit(lr.n, "-")[[1]][2];
    l.run.dat.n <- sub.dat[feat == l.n, .N];
    r.run.dat.n <- sub.dat[feat == r.n, .N];
    if(l.run.dat.n == 0 | r.run.dat.n == 0){;
      bad.m <- c(bad.m, lr.n);
    };
  };
  m.to.mod <- setdiff(m.to.mod, bad.m);


  # subset to relevant data and run
  for(lr.n in m.to.mod){


    cat('\n\n\n');
    for(i in 1:3){print(glue('ON FEAT: {lr.n}: {which(m.to.mod == lr.n)} of {length(m.to.mod)}\n'))};
    cat('\n\n\n');

    l.n <- strsplit(lr.n, "-")[[1]][1];
    r.n <- strsplit(lr.n, "-")[[1]][2];

    l.run.dat <- sub.dat[feat == l.n, ];
    set(l.run.dat, j = "feat.present", value = as.integer(l.run.dat[, feat.present]));
    set(l.run.dat, j = "total.present", value = as.integer(l.run.dat[, total.present]));
    setnames(l.run.dat,
             c("feat", "feat.count", "feat.observed", "feat.present"),
             c("feat.l", "feat.l.count", "feat.l.observed", "feat.l.present"));
    r.run.dat <- sub.dat[feat == r.n, .(x, y, feat, feat.count, feat.observed, feat.present)];
    set(r.run.dat, j = "feat.present", value = as.integer(r.run.dat[, feat.present]));
    setnames(r.run.dat,
             c("feat", "feat.count", "feat.observed", "feat.present"),
             c("feat.r", "feat.r.count", "feat.r.observed", "feat.r.present"));
    run.dat <- merge(l.run.dat, r.run.dat, by = c("x", "y"));

    source(file.path(c.d, "02-array-run-model-nonlinear.R"))
  }
}


## # plot feature and total
## par(mfrow = c(1, 2))
## fields::quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
##                    main = glue('{lr.n}'))
## fields::quilt.plot(run.dat[, x], run.dat[, y], run.dat[, total.count],
##                    main = 'Total')
## par(mfrow = c(1, 1))

# run the model

###############
## save outs ##
###############
## outputs saved in array-run-model.R
