
## this file contains functions to perform spatial kriging using INLA
## and implementing a dual conditional likelihood where we model:
## occurence v non-occurence with a binomial
## conditional on occurence, the standardized (positive) response is modeled with a gamma dist.

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

# library(rSPDE)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(ggplot2)
library(INLA)
library(splancs)
library(fmesher)
library(fields)

# dirs
code.dir     <- "./code-spatial-smoothing"
out.top.dir  <- "./data-outputs"
in.dir       <- "./data-inputs"
out.sub.dir <- NULL ## sub.dir to store output, if NULL, will use current date-string

source(file.path(code.dir, 'make.pred.grid.R'))

if(is.null(out.sub.dir)){
  out.sub.dir <- Sys.Date()
}

out.dir <- file.path(out.top.dir, out.sub.dir)
dir.create(out.dir, recursive = T)

########################################
## load data
########################################
load(file.path(in.dir, "LR.output.2023-11-08.Robj")) ## for locations
load(file.path(in.dir, "fgf1.fgfr2.matrix.Robj"))    ## for signalling response

xy.obs.orig  <- LR.output$brain1$location.info

# HACK -> TODO remove this when we've fixed xy
xy.obs <- as.matrix(cbind(x = xy.obs.orig$y, y = -xy.obs.orig$x))

lig.obs <- LR.output$brain1$ligand.info
lig.obs <- lig.obs["Fgf1", ]
rec.obs <- LR.output$brain1$receptor.info
rec.obs <- rec.obs["Fgfr2", ]

# set data to smooth
tx.obs <- fgf1.fgfr2.matrix |> as.numeric()

# fgf1-fgfr1
# tx.obs <- lig.obs


#############################
## make dual observations
## z_i: occurrence~binomial(pi_i, n_i=1)
## y_i: response~gamma(a_i, b_i)
##
## logit(pi_i) = b0^z + b1*zeta_i + u_i
##
## zeta_i: shared spatial component
## u_i binom specific spatial
##
##
## E(y_i) = a_i/b_i = mu_i
## Var(y_i) = a_i/b_i^2 = mu_i^2/tau (tau=precision param)
##
## log(mu_i)=b0^y + zeta_i
##
############################

z <- as.numeric(tx.obs > 0)
y <- ifelse(tx.obs > 0, tx.obs, NA)

# check we did things correctly
mean(is.na(y) == (z == 0))

############################
## make mesh, spde, A
############################

# mesh params. tweak to increase resolution of approx
# mesh.edge.params <- c(2, 10) # inner max edge, outer max edge

# prior matern params
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
matern.pri <- c(10, .95, 10, .05) ## a, b, c, d


# predict resolution (number of pixels in x and y)
pred.grid.res <- 200




coords <- as.matrix(xy.obs)

# since the shape may not be convex, use the nonconvex.hull to avoid many small triangles
domain <- fm_nonconvex_hull(coords, -0.03, -0.05, resolution = c(100, 100))
# TODO - make boundary extend less far to min nodes and computation
mesh <- fm_mesh_2d(boundary = domain, max.edge = c(3, 10), cutoff = 1)
plot(mesh, asp = 1, main = "")
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")

spde <- inla.spde2.pcmatern(mesh=mesh, alpha =2,
                            ## constr = TRUE, integrate-to-zero constraint
                            prior.range = matern.pri[1:2],
                            prior.sigma = matern.pri[3:4])

A.est <- inla.spde.make.A(mesh = mesh, loc = coords)

save(file = file.path(out.top.dir, out.sub.dir, "pre-fit-obj.Rdata"), list = ls())
load(file.path(out.top.dir, out.sub.dir, "pre-fit-obj.Rdata"))

############################
## create data stacks
############################
stk.y <- inla.stack(data = list(amount = y, # for single model
                                alldata = cbind(y, NA)), # for joint model
                    A = list(A.est, 1),
                    effects = list(
                      list(y.field = 1:spde$n.spde),
                      list(y.intercept = rep(1, length(y)))),
                    tag = "est.y"
                    )


stk.z <- inla.stack(data = list(occurrence = z, # for single model
                                alldata = cbind(NA, z)), # for joint model
                    A = list(A.est, 1),
                    effects = list(
                      list(z.field = 1:spde$n.spde,
                           zc.field = 1:spde$n.spde),
                      list(z.intercept = rep(1, length(z)))),
                    tag = "est.z"
                    )

# spatial response only
f.y <- amount ~ -1 + y.intercept + f(y.field, model = spde)
out.y <- inla(f.y, family = "gamma",
              data = inla.stack.data(stk.y),
              control.predictor = list(compute = TRUE, A = inla.stack.A(stk.y)),
              control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
              control.inla = list(strategy = "laplace"))

# spatial binomial only
f.z <- occurrence ~ -1 + z.intercept + f(z.field, model = spde)
out.z <- inla(f.z, family = "binomial",
              data = inla.stack.data(stk.z),
              control.predictor = list(A = inla.stack.A(stk.z), compute = TRUE),
              control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
              control.inla = list(strategy = "laplace"))

#######################
# spatial joint model #
# the real deal       #
#######################
stk.yz <- inla.stack(stk.y, stk.z)
f.yz <- alldata ~ -1 + y.intercept + z.intercept +
  f(y.field, model = spde) + f(z.field, model = spde) +
  f(zc.field, copy = "y.field", fixed = FALSE)
out.yz <- inla(f.yz, family = c("gamma", "binomial"),
               data = inla.stack.data(stk.yz),
               control.predictor = list(A = inla.stack.A(stk.yz), compute = TRUE),
               control.compute = list(dic = TRUE, config = TRUE, return.marginals.predictor=TRUE),
               control.inla = list(strategy = "laplace"))

# non-spatial joint
f.nospatial <- alldata ~ -1 + z.intercept + y.intercept
out.nospatial <- inla(f.nospatial, family=c("gamma", "binomial"),
                      data=inla.stack.data(stk.yz),
                      control.predictor=list(A=inla.stack.A(stk.yz)),
                      control.compute=list(dic=TRUE),
                      control.inla=list(strategy="laplace"))

# joint but not shared spatial component
f.noshared <- alldata ~ -1 + y.intercept + z.intercept +
  f(y.field, model=spde) + f(z.field, model=spde)
out.noshared <- inla(f.noshared, family=c("gamma", "binomial"),
                     data=inla.stack.data(stk.yz),
                     control.predictor=list(A=inla.stack.A(stk.yz)),
                     control.compute=list(dic=TRUE),
                     control.inla=list(strategy="laplace"))

save(file = file.path(out.top.dir, out.sub.dir, "dual-likelihood-test-fit.Rdata"), list = ls())
load(file.path(out.top.dir, "2024-02-15/dual-likelihood-test-fit.Rdata"))

## perform model comparison between the three models using DIC
# note: need to be careful with the DIC when using multiple likelihoods (joint)
#       must sum local DIC for each obs
# we want lower DIC which indicates better fit

rbind(nospatial=tapply(out.nospatial$dic$local.dic, out.nospatial$dic$family, sum),
      noshared=tapply(out.noshared$dic$local.dic, out.noshared$dic$family, sum),
      separate=c(out.y$dic$dic, out.z$dic$dic),
      joint=tapply(out.yz$dic$local.dic, out.yz$dic$family, sum))


## calculate the expected value of the posterior marginal
## distributions on the response scale
idy <- which(z == 1)
exp.y <- sapply(out.yz$marginals.linear.predictor[idy],
                function(m){inla.emarginal(exp, m)})

# is this which call right? it's what is in the INLA book pg 276
idz <- which(!is.na(z))
exp.z <- sapply(out.yz$marginals.linear.pred[length(tx.obs) + idz],
                function(m){inla.emarginal(inla.link.invlogit, m)})

# compare these values to the observed mean for the signal occurrence
# and amount -- a nice simple sanity check to make sure everything is
# as expected
c(yPositive.obs = mean(tx.obs[which(tx.obs > 0)]),
  yPositive.pred = mean(exp.y))

c(z.obs = mean(z, na.rm = T), z.pred = mean(exp.z))

## plot the posterior marginal distributions for y.int, z.int,
## tau for separate and joint models, and beta_1 for joint
## as in figure 8.10 pg 277
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=2:0)
plot(inla.smarginal(out.yz$marginals.fixed[[1]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.fixed[[1]][,2], out.y$marginals.fixed[[1]][,2])), xlab=expression(b[0]^y))
lines(inla.smarginal(out.y$marginals.fixed[[1]]), lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.fixed[[2]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.fixed[[2]][,2],out.z$marginals.fixed[[1]][,2])), xlab=expression(b[0]^z))
lines(inla.smarginal(out.z$marginals.fixed[[1]]), lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.hyperpar[[1]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.hy[[1]][,2], out.y$marginals.hy[[1]][,2])), xlab=expression(tau))
lines(inla.smarginal(out.y$marginals.hyperpar[[1]]),lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.hyperpar[[6]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.hy[[6]][,2])),  xlab=expression(beta[1]))

## to check if the occurrence and amount of signaling are spatially
## dependent, we can test the significance of the spatial random
## effects zeta(s) and u(s)
boundary.coords <- as(domain, "Spatial")@polygons[[1]]@Polygons[[1]]@coords # convert type
in.d <- which(inout(mesh$loc[, 1:2], boundary.coords))
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,1,0))
ordxy <- intersect(order(out.yz$summary.random$y.field$mean), in.d)
plot(out.yz$summary.random$y.field$mean[ordxy], type="l", ylab=expression(xi[i]), ylim=range(out.yz$summary.random$y.field[, 4:6]))
for (i in c(4,6)) lines(out.yz$summary.random$y.field[ordxy,i], lty=2)
abline(h=0, lty=3)
ordxz <- intersect(order(out.yz$summary.random$z.field$mean), in.d)
plot(out.yz$summary.random$z.field$mean[ordxz], type="l", ylab=expression(u[i]), ylim=range(out.yz$summary.random$z.field[, 4:6]))
for (i in c(4,6)) lines(out.yz$summary.random$z.field[ordxz,i], lty=2)
abline(h=0, lty=3)

#--- Significance of the spatial effects ---#
c(n.nodes.in.d=length(in.d),
  expected.out=0.05*length(in.d),
  observed.y.out=sum(out.yz$summary.random$y.field[in.d,4]>0 |
                       out.yz$summary.random$y.field[in.d,6]<0),
  observed.z.out=sum(out.yz$summary.random$z.field[in.d,4]>0 |
                       out.yz$summary.random$z.field[in.d,6]<0))


## extract the spatial effects
y.field <- inla.spde2.result(out.yz, "y.field", spde)
z.field <- inla.spde2.result(out.z, "z.field", spde)

par(mfrow=c(2,2))
plot.default(y.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
plot.default(y.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
plot.default(z.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
plot.default(z.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')

#### plot spatial random fields

## make projection grid

# predict resolution (number of pixels in x and y)
pred.grid.res <- 200

xy.pred <- make.pred.grid(xy.obs = xy.obs,
                          retain.shape = T,
                          res = pred.grid.res)

A.pred <- inla.spde.make.A(
  mesh = mesh,
  loc = xy.pred)

rf.grid <- list(xi.mean = fm_evaluate(mesh, field = out.y$summary.random$y.field$mean, loc = xy.pred),
                xi.sd = fm_evaluate(mesh, field = out.y$summary.random$y.field$sd, loc = xy.pred),
                obs = tx.obs,
                u.mean = fm_evaluate(mesh, field = out.yz$summary.random$z.field$mean, loc = xy.pred),
                u.sd = fm_evaluate(mesh, field = out.yz$summary.random$z.field$sd, loc = xy.pred),
                occ.obs = z)

xy.in <- inout(xy.pred, boundary.coords)
for (i in c(1:2, 4:5)) rf.grid[[i]][!xy.in] <- NA

## TODO make this plot with intercept

pdf(file.path(out.top.dir, out.sub.dir, "two-spatial-field-effects.pdf"),
    width = 12, height = 8)
par(mfrow=c(2,3), mar=c(0,0,0,0))
for (i in 1:6) {
  if(grepl("obs", names(rf.grid)[i])){
    quilt.plot(xy.obs, rf.grid[[i]],
               nx = 60, ny = 50,
               nlevel=100, col=viridis(100),
               legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
               xaxt = 'n', yaxt = 'n')
  }else{
    quilt.plot(xy.pred, rf.grid[[i]],
               nx = pred.grid.res * .9, ny = pred.grid.res * .9,
               nlevel=101, col=viridis(100),
               legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
               xaxt = 'n', yaxt = 'n')
  }
  lines(boundary.coords)
}
dev.off()

#####################################
#--- Prediction of the responses ---#
#####################################

s1 <- inla.posterior.sample(n=1,result=out.yz)
names(s1[[1]])
grep("y.intercept", rownames(s1[[1]]$latent), fixed=TRUE)

ids <- lapply(c('y.intercept', 'z.intercept', 'y.field', 'z.field', 'zc.field'), function(x) grep(x, rownames(s1[[1]]$latent), fixed=TRUE))
pred.y.f <- function(s) exp(s$latent[ids[[1]], 1] + s$latent[ids[[3]], 1])
pred.z.f <- function(s) 1/(1 + exp(-(s$latent[ids[[2]], 1] + s$latent[ids[[4]], 1] + s$latent[ids[[5]], 1]))) ## NOTE! the latent zc.field is already beta.zc*y.field

s1000 <- inla.posterior.sample(1000, out.yz)

prd.y <- sapply(s1000, pred.y.f)
prd.z <- sapply(s1000, pred.z.f)

prd <- list(y.mean=fm_evaluate(mesh, field=rowMeans(prd.y), loc = xy.pred))
prd$y.sd <- fm_evaluate(mesh, field=apply(prd.y, 1, sd), loc = xy.pred)
prd$obs <- tx.obs
prd$z.mean <- fm_evaluate(mesh, field=rowMeans(prd.z), loc = xy.pred)
prd$z.sd <- fm_evaluate(mesh, field=apply(prd.z, 1, sd), loc = xy.pred)
prd$occ.obs <- z
for (j in c(1:2, 4:5)) prd[[j]][!xy.in] <- NA

# *** Code for Figure 8.14
png(file.path(out.top.dir, out.sub.dir, "fgf1-response-occurrence-field-ests.png"),
    width = 12, height = 8, units = "in", res = 500)
par(mfcol=c(2,3), mar=c(0,0,0,0))
for (i in 1:6) {
  if(grepl("obs", names(prd)[i])){
    quilt.plot(xy.obs, prd[[i]],
               nx = 60, ny = 50,
               nlevel=100, col=viridis(100),
               legend.mar=7, legend.args=list(text=names(prd)[i], line=0),
               xaxt = 'n', yaxt = 'n')
  }else{
    quilt.plot(xy.pred, prd[[i]],
               nx = pred.grid.res * .9, ny = pred.grid.res * .9,
               nlevel=101, col=viridis(100), axes = FALSE,
               legend.mar=7, legend.args=list(text=names(prd)[i], line=0),
               xaxt = 'n', yaxt = 'n')
  }
  lines(boundary.coords)
}
dev.off()
