
## this file contains functions to perform spatial kriging using INLA and SPDE

## notes:
#
# 1) input takes in:
#    - xy location matrix (named cols x and y)
#    - a vector of values at each location
#
# 2) outputs the same (but with an updated xy loc mat)
#

## steps:
# 0) set user params
# 1) load data
# 2) make FEM triangulation
# 3) prep obj for model fit
# 4) fit model
# 5) predict
# 6) save output

########################################
####  set user-specified params here
########################################

library(rSPDE)
library(ggplot2)
library(INLA)
library(splancs)
library(fmesher)

# dirs
code.dir     <- "./code-spatial-smoothing"
out.top.dir  <- "./data-outputs"
in.dir       <- "./data-inputs"
out.subdir <- NULL ## subdir to store output, if NULL, will use current date-string

########################################
## load data
########################################
load(file.path(in.dir, "LR.output.2023-11-08.Robj"))

xy.obs  <- LR.output$brain1$location.info
lig.obs <- LR.output$brain1$ligand.info
lig.obs <- lig.obs["Fgf1", ]
rec.obs <- LR.output$brain1$receptor.info
rec.obs <- rec.obs["Fgfr2", ]

# set data to smooth
tx.obs <- lig.obs



#############################
## make mesh
############################
coords <- as.matrix(xy.obs)
prdomain <- fm_nonconvex_hull(as.matrix(xy.obs), -0.03, -0.05, resolution = c(100, 100))
prmesh <- fm_mesh_2d(boundary = prdomain, max.edge = c(3, 5), cutoff = 0.2)
plot(prmesh, asp = 1, main = "")
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")


#Create the observation matrix
Abar <- rspde.make.A(mesh = prmesh, loc = coords)

#Create the rspde model object
rspde_model <- rspde.matern(mesh = prmesh)

#Create the index and inla.stack object
mesh.index <- rspde.make.index(name = "field", mesh = prmesh)
stk.dat <- inla.stack(
  data = list(y = tx.obs), A = list(Abar, 1), tag = "est",
  effects = list(c(mesh.index),
                 list(long = inla.group(coords[, 1]),
                      lat = inla.group(coords[,2]),
                      #seaDist = inla.group(seaDist),
                      Intercept = 1)))

#Create the formula object and fit the model
f.s <- y ~ -1 + Intercept +  #f(seaDist, model = "rw1") +
  f(field, model = rspde_model)

# fit with INLA
rspde_fit <- inla(f.s, family = "Gamma", data = inla.stack.data(stk.dat),
            verbose = FALSE,
            control.inla=list(int.strategy='eb'),
            control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE))






# mesh params. tweak to increase resolution of approx
mesh.edge.params <- c(2, 10) # inner max edge, outer max edge

# prior matern params
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
matern.pri <- c(2, .95, 100, .05) ## a, b, c, d


# predict resolution (number of pixels in x and y)
pred.grid.res <- 200


########################################
## set environment
########################################
setwd(code.dir)

require(data.table)
require(glue)
require(fields) # require(rSPDE)
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
require(INLA)

source('make.pred.grid.R')

if(is.null(out.subdir)){
  out.subdir <- Sys.Date()
}

out.dir <- file.path(out.top.dir, out.subdir)
dir.create(out.dir, recursive = T)

########################################
## FEM triangulation
########################################

# get locs for array to be analyzed
xy.obs <- LR.output$brain1$location.info

# HACK -> TODO remove this when we've fixed xy
xy.obs <- as.matrix(cbind(x = xy.obs$y, y = -xy.obs$x))

mesh_s <- inla.mesh.2d(loc.domain = xy.obs,
                       max.edge = mesh.edge.params)
plot(mesh_s);points(xy.obs, col = 2, pch = 16)

########################################
## fit model
########################################

# make inla input objects

design.mat <- data.frame(int = rep(1, length(tx.obs)))#,
                            #dt[, inla.covs, with=F],
                            #clust.id = 1:nrow(dt))

A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = xy.obs
)

spde <- inla.spde2.pcmatern(mesh=mesh_s, alpha =2,
                            ## constr = TRUE, integrate-to-zero constraint
                            prior.range = matern.pri[1:2],
                            prior.sigma = matern.pri[3:4])

space   <- inla.spde.make.index("space",
                                n.spde = spde$n.spde)

stack.obs <- inla.stack(tag='est',
                        data=list(Y=as.numeric(tx.obs)), ## response
                        A=list(A, 1), ## proj matrix, not sure what the 1 is for
                        effects=list(
                          space,
                          design.mat))

inla.f <- formula('Y ~ -1 + int + f(space, model=spde)')

fit <- inla(inla.f, data = inla.stack.data(stack.obs),
            control.predictor = list(A = inla.stack.A(stack.obs),
                                     compute = FALSE),
            #control.fixed = list(expand.factor.strategy = 'inla'),
            control.inla = list(strategy = "gaussian",
                                int.strategy = "eb"),
            control.compute=list(config = TRUE),
            #control.family = list(hyper = list(prec = list(prior="pc.prec",
            #                                               param = norm.prec.pri))),
            family = "normal",
            num.threads = 4, #
            # scale = dt$N,
            verbose = FALSE, ## this must be false to get the logfile
            keep = FALSE)

########################################
## predict
########################################


xy.pred <- make.pred.grid(xy.obs = xy.obs,
                          retain.shape = T,
                          res = pred.grid.res)

A.pred <- inla.spde.make.A(
  mesh = mesh_s,
  loc = xy.pred)

sfield_nodes <- fit$summary.random$space['mean']
inla.int <- as.numeric(fit$summary.fixed['mean'])
field <- as.numeric((A.pred %*% as.data.frame(sfield_nodes)[, 1]) + inla.int)

quilt.plot(xy.pred, field)
quilt.plot(xy.obs, tx.obs , nx = 50, ny = 50)

# TODO - why is the amplitude getting so severely flattened?
