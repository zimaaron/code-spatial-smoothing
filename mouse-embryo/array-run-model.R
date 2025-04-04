
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

##########
## pkgs ##
##########

# library(rSPDE)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(ggplot2)
library(INLA)
## INLA:::inla.binary.install("Rocky Linux-8") ## seems like calling this once per INLA install is ok
## to fix the following error: /home/aeoz001/R/x86_64-pc-linux-gnu-library/4.2/INLA/bin/linux/64bit/inla.mkl: /lib64/libm.so.6: version `GLIBC_2.29' not found (required by /home/aeoz001/R/x86_64-pc-linux-gnu-library/4.2/INLA/bin/linux/64bit/libicuuc.so.66)
## TODO ask bisonnet if this can be fixed on the cluster
library(splancs)
library(fmesher)
library(fields)
require(purrr)
require(SeuratObject)

source(file.path(c.d, 'make.pred.grid.v2.R'))

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

z <- as.numeric(lr.d > 0)
y <- ifelse(lr.d > 0, lr.d, NA)

## check we did things correctly
#mean(is.na(y) == (z == 0)) ## should equal 1

############################
## make mesh, spde, A
############################

# mesh params. tweak to increase resolution of approx
# mesh.edge.params <- c(2, 10) # inner max edge, outer max edge

# prior matern params
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
matern.pri <- c(40, .95, 10, .05) ## a, b, c, d


# predict resolution (number of pixels in x and y)
pred.grid.res <- 200

# load:
# domain: nonconvex hull
# mesh: SPDE mesh
# boundary.coords: boundary coords for plotti
# in.d: which coords are inside the domain
load(file.path(i.d, glue("pre-fit-obj-{time.i}.Rdata")))


spde <- inla.spde2.pcmatern(mesh=mesh, alpha =2,
                            ## constr = TRUE, integrate-to-zero constraint
                            prior.range = matern.pri[1:2],
                            prior.sigma = matern.pri[3:4])

A.est <- inla.spde.make.A(mesh = mesh, loc = coords)


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

## make a quiet version of INLA to parse logs
#inla.q <- quietly(inla)

# spatial response only


spat.start.time <- toc()[3] / 60
f.y <- amount ~ -1 + y.intercept + f(y.field, model = spde)
out.y <- inla(f.y, family = "gamma",
              data = inla.stack.data(stk.y),
              control.predictor = list(compute = TRUE, A = inla.stack.A(stk.y)),
              control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
              control.inla = list(strategy = "laplace"))

spat.start.elapsed.time <- glue("time to fit spatial response only at time {time.i} is: {round(toc()[3]/60 - spat.start.time,2)} min")
# time.mat[1, time.i] <- spat.start.elapsed.time # for profile timing
print(spat.start.elapsed.time)

## # spatial binomial only
## f.z <- occurrence ~ -1 + z.intercept + f(z.field, model = spde)
## out.z <- inla(f.z, family = "binomial",
##               data = inla.stack.data(stk.z),
##               control.predictor = list(A = inla.stack.A(stk.z), compute = TRUE),
##               control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
##               control.inla = list(strategy = "laplace"))

## print(glue("time to fit spatial binom only is: {toc()[3]}"))

#######################
# spatial joint model #
# the real deal       #
#######################

full.start.time <- toc()[3] / 60
stk.yz <- inla.stack(stk.y, stk.z)
f.yz <- alldata ~ -1 + y.intercept + z.intercept +
  f(y.field, model = spde) + f(z.field, model = spde) +
  f(zc.field, copy = "y.field", fixed = FALSE)
out.yz <- inla(f.yz, family = c("gamma", "binomial"),
               data = inla.stack.data(stk.yz),
               control.predictor = list(A = inla.stack.A(stk.yz), compute = TRUE),
               control.compute = list(dic = TRUE, config = TRUE, return.marginals.predictor=TRUE),
               control.inla = list(strategy = "laplace"),
               verbose = F)

full.start.elapsed.time <- (glue("time to fit full model at time {time.i} is: {round((toc()[3]/60 - full.start.time) ,2)} min"))
# time.mat[2, time.i] <- full.start.elapsed.time # for profile timing
print(full.start.elapsed.time)

## # non-spatial joint
## f.nospatial <- alldata ~ -1 + z.intercept + y.intercept
## out.nospatial <- inla(f.nospatial, family=c("gamma", "binomial"),
##                       data=inla.stack.data(stk.yz),
##                       control.predictor=list(A=inla.stack.A(stk.yz)),
##                       control.compute=list(dic=TRUE),
##                       control.inla=list(strategy="laplace"),
##                       verbose = F)

## # joint but not shared spatial component
## f.noshared <- alldata ~ -1 + y.intercept + z.intercept +
##   f(y.field, model=spde) + f(z.field, model=spde)
## out.noshared <- inla(f.noshared, family=c("gamma", "binomial"),
##                      data=inla.stack.data(stk.yz),
##                      control.predictor=list(A=inla.stack.A(stk.yz)),
##                      control.compute=list(dic=TRUE),
##                      control.inla=list(strategy="laplace"))

##save(file = file.path(out.top.dir, out.sub.dir, "dual-likelihood-test-fit.Rdata"), list = ls())
##load(file.path(out.top.dir, "2024-02-15/dual-likelihood-test-fit.Rdata"))

## perform model comparison between the three models using DIC
# note: need to be careful with the DIC when using multiple likelihoods (joint)
#       must sum local DIC for each obs
# we want lower DIC which indicates better fit

## TODO rewrite all numerical outputs into one output row -> csv
## write.csv(rbind(nospatial=tapply(out.nospatial$dic$local.dic, out.nospatial$dic$family, sum),
##                 noshared=tapply(out.noshared$dic$local.dic, out.noshared$dic$family, sum),
##                 separate=c(out.y$dic$dic, out.z$dic$dic),
##                 joint=tapply(out.yz$dic$local.dic, out.yz$dic$family, sum)),
##           file.path(o.d, glue("{lr.n}-DIC-comp.csv")))

## calculate the expected value of the posterior marginal
## distributions on the response scale
idy <- which(z == 1)
exp.y <- sapply(out.yz$marginals.linear.predictor[idy],
                function(m){inla.emarginal(exp, m)})

## TODO is this which call right? it's what is in the INLA book pg 276
idz <- which(!is.na(z))
## NOTE 28aug2024 some of the marginals for z (pre transform) have Inf height
## this breaks inla.emarginal, this the tryCatch. But it's not a big deal since we
## exp.z is only used to store mean(exp.z) as a summary metric
exp.z <- tryCatch({sapply(out.yz$marginals.linear.predictor[length(lr.d) + idz],
                          function(m){inla.emarginal(inla.link.invlogit, m)})
                  },
                  error = function(err){
                    NA
                  })


## # compare these values to the observed mean for the signal occurrence
## # and amount -- a nice simple sanity check to make sure everything is
## # as expected
## means.comp <- matrix(c(yPositive.obs = mean(lr.d[which(lr.d > 0)]),
##                        yPositive.pred = mean(exp.y),
##                        z.obs = mean(z, na.rm = T), z.pred = mean(exp.z)),
##                      nrow = 2, byrow = T)
## colnames(means.comp) <- c("obs", "pred")
## rownames(means.comp) <- c("yPositive", "z")
## write.csv(means.comp, file.path(o.d, glue("{lr.n}-obs-est-means-comp.csv")))

## plot the posterior marginal distributions for y.int, z.int,
## tau for separate and joint models, and beta_1 for joint
## as in figure 8.10 pg 277
#print("starting first plot")

png(file.path(o.d, glue("{lr.n}-{time.i}-post-marg-dist.png")),
    width = 8, height = 8, units = "in", res = 300)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=2:0)
plot(inla.smarginal(out.yz$marginals.fixed[[1]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.fixed[[1]][,2], out.y$marginals.fixed[[1]][,2])),
     xlab=expression(b[0]^y))
lines(inla.smarginal(out.y$marginals.fixed[[1]]), lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.fixed[[2]]), type='l', ylab='Density',
    #ylim=c(0,max(out.yz$marginals.fixed[[2]][,2], out.z$marginals.fixed[[1]][,2])),
                  xlab=expression(b[0]^z))
#lines(inla.smarginal(out.z$marginals.fixed[[1]]), lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.hyperpar[[1]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.hy[[1]][,2], out.y$marginals.hy[[1]][,2])),
     xlab=expression(tau))
lines(inla.smarginal(out.y$marginals.hyperpar[[1]]),lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.hyperpar[[6]]), type='l', ylab='Density',
     ylim=c(0,max(out.yz$marginals.hy[[6]][,2])),  xlab=expression(beta[1]))
dev.off()

## to check if the occurrence and amount of signaling are spatially
## dependent, we can test the significance of the spatial random
## effects zeta(s) and u(s)
png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-effects-ordered.png")),
    width = 8, height = 8, units = "in", res = 300)
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,1,0))
ordxy <- intersect(order(out.yz$summary.random$y.field$mean), in.d) # in.d(omain)
plot(out.yz$summary.random$y.field$mean[ordxy], type="l", ylab=expression(xi[i]), ylim=range(out.yz$summary.random$y.field[, 4:6]))
for (i in c(4,6)) lines(out.yz$summary.random$y.field[ordxy,i], lty=2)
abline(h=0, lty=3)
ordxz <- intersect(order(out.yz$summary.random$z.field$mean), in.d)
plot(out.yz$summary.random$z.field$mean[ordxz], type="l", ylab=expression(u[i]), ylim=range(out.yz$summary.random$z.field[, 4:6]))
for (i in c(4,6)) lines(out.yz$summary.random$z.field[ordxz,i], lty=2)
abline(h=0, lty=3)
dev.off()

## #--- Significance of the spatial effects ---#
spat.sig <- data.frame(n.nodes.in.d=length(in.d),
                       expected.out=0.05*length(in.d),
  observed.y.out=sum(out.yz$summary.random$y.field[in.d,4]>0 |
                       out.yz$summary.random$y.field[in.d,6]<0),
  observed.z.out=sum(out.yz$summary.random$z.field[in.d,4]>0 |
                       out.yz$summary.random$z.field[in.d,6]<0))
#write.csv(spat.sig, file.path(o.d, glue("{lr.n}-spat-sig.csv")), row.names = F)

out.yz.summ <- data.frame(lr.n = lr.n,
                          time = time.i,
                          var = c("local.dic",
                                  "family.dic",
                                  "yPos-obs",
                                  "yPos-pred",
                                  "z-obs",
                                  "z-pred",
                                  colnames(spat.sig)
                                  ),
                          val = c(tapply(out.yz$dic$local.dic, out.yz$dic$family, sum),
                                  # the tapply returns both local and family dic
                                  mean(lr.d[which(lr.d > 0)]),
                                  mean(exp.y),
                                  mean(z, na.rm = T),
                                  mean(exp.z),
                                  as.numeric(spat.sig[1, ]))
                         )

write.csv(out.yz.summ, file.path(o.d, glue("{lr.n}-{time.i}-output-summary.csv")), row.names = F)


## extract the spatial effects
y.field <- inla.spde2.result(out.yz, "y.field", spde)
z.field <- inla.spde2.result(out.yz, "z.field", spde)

png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-effects-marginals.png")),
    width = 8, height = 8, units = "in", res = 300)
par(mfrow=c(2,2))
plot.default(y.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
plot.default(y.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
plot.default(z.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
plot.default(z.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')
dev.off()

#### plot spatial random fields

## make projection grid

# predict resolution (number of pixels in x and y)
xy.pred <- make.pred.grid(xy.obs = xy.obs,
                          retain.shape = T,
                          res = pred.grid.res,
                          boundary.poly = boundary.coords)
if(j.i  == 1){
  write.csv(xy.pred, file.path(o.d, glue("xy-pred-locs-{time.i}.csv")), row.names = F)
}

A.pred <- inla.spde.make.A(
  mesh = mesh,
  loc = xy.pred)

rf.grid <- list(xi.mean = fm_evaluate(mesh, field = out.y$summary.random$y.field$mean, loc = xy.pred),
                xi.sd = fm_evaluate(mesh, field = out.y$summary.random$y.field$sd, loc = xy.pred),
                obs = lr.d,
                u.mean = fm_evaluate(mesh, field = out.yz$summary.random$z.field$mean, loc = xy.pred),
                u.sd = fm_evaluate(mesh, field = out.yz$summary.random$z.field$sd, loc = xy.pred),
                occ.obs = z)

xy.in <- inout(xy.pred, boundary.coords)
for (i in c(1:2, 4:5)) rf.grid[[i]][!xy.in] <- NA


png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-fields.png")),
    width = 12, height = 8, units = "in", res = 500)
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
  lines(boundary.coords, col = "orange", lwd = 2)
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
prd.y.cond.z <- matrix(rbinom(n = length(prd.z), size = 1, p = prd.z),
                       ncol = ncol(prd.z))* prd.y

# project all draws to pred grid
prd.y.proj <- apply(prd.y, 2,
                    function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
prd.z.proj <- apply(prd.z, 2,
                    function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
prd.y.cond.z.proj <- apply(prd.y.cond.z, 2,
                           function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})


# always write the means
write.csv(as.matrix(rowMeans(prd.y.proj), ncol = 1),
          file.path(o.d, glue("{lr.n}-{time.i}-resp-est.csv")), row.names = F)
write.csv(as.matrix(rowMeans(prd.y.cond.z.proj), ncol = 1),
          file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-est.csv")), row.names = F)
write.csv(as.matrix(apply(prd.y.cond.z.proj, 1, sd), ncol = 1),
          file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-sd.csv")), row.names = F)

# only save some of the complete draws -- these are big files (.5Gb)
# and only save the "final" field draws
if(lr.n %in% lr.to.save.draws){
  write.csv(prd.y.cond.z.proj,
            file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-draws.csv")), row.names = F)
}

# prepare to plot obs, est, sd for y, z, y|z

prd <- list(resp.obs = lr.d)
prd$y.mean <- fm_evaluate(mesh, field=rowMeans(prd.y), loc = xy.pred)
prd$y.sd <- fm_evaluate(mesh, field=apply(prd.y, 1, sd), loc = xy.pred)

prd$occ.obs <- z
prd$z.mean <- fm_evaluate(mesh, field=rowMeans(prd.z), loc = xy.pred)
prd$z.sd <- fm_evaluate(mesh, field=apply(prd.z, 1, sd), loc = xy.pred)

prd$cond.resp.obs <- lr.d
prd$y.cond.z.mean <- fm_evaluate(mesh, field=rowMeans(prd.y.cond.z), loc = xy.pred)
prd$y.cond.z.sd <- fm_evaluate(mesh, field=apply(prd.y.cond.z, 1, sd), loc = xy.pred)

# knock out places where we don't want to plot/predict
for (j in c(2:3, 5:6, 8:9)) prd[[j]][!xy.in] <- NA

titles <- c("Y Obs", "Y Est", "Y SD",
            "Z Obs", "Z Est", "Z SD",
            "Y|Z Obs", "Y|Z Est", "Y|Z SD")

# *** Code for Figure 8.14
png(file.path(o.d, glue("{lr.n}-{time.i}-response-occurrence-fields.png")),
    width = 12, height = 12, units = "in", res = 500)
par(mfcol=c(3,3), mar=c(0,0,0,1))
for (i in 1:length(prd)) {

  # set xlim and ylim for consistency
  xl <- range(xy.obs[, 1]) + c(-.05, .05) * (max(xy.obs[, 1]) - min(xy.obs[, 1]))
  yl <- range(xy.obs[, 2]) + c(-.05, .05) * (max(xy.obs[, 2]) - min(xy.obs[, 2]))

  ## first, set z/color ranges and flag large predicted values in estimates
  if(i %in% c(1, 2, 7, 8)){

    zl <- c(0, max(prd$resp.obs) * 1.01) ## mean and data obs plots
    zcols <- c(viridis(100), "#b33609")

    # this should be a shrinkage/smoothing model. outputs above max observed are truncated and
    # plotted as dark red
    prd[[i]][prd[[i]] > max(prd$resp.obs) * 1.01] <- max(prd$resp.obs) * 1.01

  }else if(i %in% c(4, 5)){## prob of non-zero response mean/est and obs

    zl <- 0:1
    if(i == 4) zcols <- viridis(2)
    if(i == 5) zcols <- viridis(100)

  }else{ ## sd plots, no scaling

    zl <- c(0, max(prd[[i]], na.rm = T))
    zcols <- turbo(100)

  }

  ## second, plot the model outputs differently than the raw data

  if(grepl("obs", names(prd)[i])){
    bubblePlot(xy.obs[, 1], xy.obs[, 2], prd[[i]],
               col=zcols,
               legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
               legend.args=list(text="", line=0),
               xaxs = 'i', yaxs = 'i', # o.w. R adds 4% to pad axes and
               # then plot this won't align with quilt.plot
               xaxt = "n", yaxt = "n", bty = "n",
               xlim = xl, ylim = yl, zlim = zl)
  }else{
    quilt.plot(xy.pred, prd[[i]],
               nx = pred.grid.res * .9, ny = pred.grid.res * .9,
               nlevel=length(zcols), col=zcols,
               legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
               legend.args=list(text="", line=0),
               xaxt = 'n', yaxt = 'n', bty = "n",
               xlim = xl, ylim = yl, zlim = zl)
  }

  ## add title and boundary lines
  text(max(xy.obs[, 1]), max(xy.obs[, 2]), adj = c(1, 0.5), titles[i], cex = 2)
  lines(boundary.coords, lwd = 2)
}
dev.off()


## png(file.path(o.d, "test.png"), width = 5, height = 5, units = "in", res = 300)
## DefaultAssay(brain)<- "NeighborhoodToCell_ALRA"
## SpatialFeaturePlot(brain, features = lr.n, slot = "data") +
##   scale_fill_viridis_c(limits = c(0.2, 1))
## dev.off()

final.elapsed.time <- (glue("final elapsed time at time {time.i} : {round(toc()[3]/60,2)} min"))
# time.mat[3, time.i] <- final.elapsed.time # for profile timing
print(final.elapsed.time)
