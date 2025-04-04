
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
## if(j.i == 1){INLA:::inla.binary.install("Rocky Linux-8")} ## seems like calling this once per INLA install is ok
## to fix the following error: /home/aeoz001/R/x86_64-pc-linux-gnu-library/4.2/INLA/bin/linux/64bit/inla.mkl: /lib64/libm.so.6: version `GLIBC_2.29' not found (required by /home/aeoz001/R/x86_64-pc-linux-gnu-library/4.2/INLA/bin/linux/64bit/libicuuc.so.66)
## TODO ask bisonnet if this can be fixed on the cluster
library(splancs)
library(fmesher)
library(fields)
require(purrr)
require(Seurat)
require(SeuratObject)

source(file.path(c.d, 'make.pred.grid.R'))

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
matern.pri <- c(10, .95, 10, .05) ## a, b, c, d


# predict resolution (number of pixels in x and y)
pred.grid.res <- 200

coords <- as.matrix(xy.obs)

# since the shape may not be convex, use the nonconvex.hull to avoid many small triangles
## domain <- fm_nonconvex_hull(coords, -0.03, -0.05, resolution = c(100, 100))
## # TODO - make boundary extend less far to min nodes and computation
## mesh <- fm_mesh_2d(boundary = domain, max.edge = c(3, 10), cutoff = 1)
## plot(mesh, asp = 1, main = "")
## points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")

# convert type for later plotting
boundary.coords <- as(domain, "Spatial")@polygons[[1]]@Polygons[[1]]@coords
# check which mesh vertices are inside the domain
in.d <- which(inout(mesh$loc[, 1:2], boundary.coords))

spde <- inla.spde2.pcmatern(mesh=mesh, alpha =2,
                            ## constr = TRUE, integrate-to-zero constraint
                            prior.range = matern.pri[1:2],
                            prior.sigma = matern.pri[3:4])

A.est <- inla.spde.make.A(mesh = mesh, loc = coords)

# save(file = file.path(out.top.dir, out.sub.dir, "pre-fit-obj.Rdata"), list = ls())
# load(file.path(out.top.dir, out.sub.dir, "pre-fit-obj.Rdata"))

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
## f.y <- amount ~ -1 + y.intercept + f(y.field, model = spde)
## out.y <- inla(f.y, family = "gamma",
##               data = inla.stack.data(stk.y),
##               control.predictor = list(compute = TRUE, A = inla.stack.A(stk.y)),
##               control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
##               control.inla = list(strategy = "laplace"))

## # spatial binomial only
## f.z <- occurrence ~ -1 + z.intercept + f(z.field, model = spde)
## out.z <- inla(f.z, family = "binomial",
##               data = inla.stack.data(stk.z),
##               control.predictor = list(A = inla.stack.A(stk.z), compute = TRUE),
##               control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
##               control.inla = list(strategy = "laplace"))

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
               control.inla = list(strategy = "laplace"),
               verbose = F)

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
## distributions on the response scalex
idy <- which(z == 1)
exp.y <- sapply(out.yz$marginals.linear.predictor[idy],
                function(m){inla.emarginal(exp, m)})

# is this which call right? it's what is in the INLA book pg 276
idz <- which(!is.na(z))
exp.z <- sapply(out.yz$marginals.linear.pred[length(lr.d) + idz],
                function(m){inla.emarginal(inla.link.invlogit, m)})

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

## png(file.path(o.d, glue("{lr.n}-{side}-post-marg-dist.png")),
##     width = 8, height = 8, units = "in", res = 300)
## par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=2:0)
## plot(inla.smarginal(out.yz$marginals.fixed[[1]]), type='l', ylab='Density',
##      ylim=c(0,max(out.yz$marginals.fixed[[1]][,2], out.y$marginals.fixed[[1]][,2])),
##      xlab=expression(b[0]^y))
## lines(inla.smarginal(out.y$marginals.fixed[[1]]), lty=2)
## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## plot(inla.smarginal(out.yz$marginals.fixed[[2]]), type='l', ylab='Density',
##     ylim=c(0,max(out.yz$marginals.fixed[[2]][,2], out.z$marginals.fixed[[1]][,2])),
##                   xlab=expression(b[0]^z))
## lines(inla.smarginal(out.z$marginals.fixed[[1]]), lty=2)
## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## plot(inla.smarginal(out.yz$marginals.hyperpar[[1]]), type='l', ylab='Density',
##      ylim=c(0,max(out.yz$marginals.hy[[1]][,2], out.y$marginals.hy[[1]][,2])),
##      xlab=expression(tau))
## lines(inla.smarginal(out.y$marginals.hyperpar[[1]]),lty=2)
## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## plot(inla.smarginal(out.yz$marginals.hyperpar[[6]]), type='l', ylab='Density',
##      ylim=c(0,max(out.yz$marginals.hy[[6]][,2])),  xlab=expression(beta[1]))
## dev.off()

## ## to check if the occurrence and amount of signaling are spatially
## ## dependent, we can test the significance of the spatial random
## ## effects zeta(s) and u(s)
## png(file.path(o.d, glue("{lr.n}-{side}-spatial-effects-ordered.png")),
##     width = 8, height = 8, units = "in", res = 300)
## par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,1,0))
## ordxy <- intersect(order(out.yz$summary.random$y.field$mean), in.d) # in.d(omain)
## plot(out.yz$summary.random$y.field$mean[ordxy], type="l", ylab=expression(xi[i]), ylim=range(out.yz$summary.random$y.field[, 4:6]))
## for (i in c(4,6)) lines(out.yz$summary.random$y.field[ordxy,i], lty=2)
## abline(h=0, lty=3)
## ordxz <- intersect(order(out.yz$summary.random$z.field$mean), in.d)
## plot(out.yz$summary.random$z.field$mean[ordxz], type="l", ylab=expression(u[i]), ylim=range(out.yz$summary.random$z.field[, 4:6]))
## for (i in c(4,6)) lines(out.yz$summary.random$z.field[ordxz,i], lty=2)
## abline(h=0, lty=3)
## dev.off()

## ## #--- Significance of the spatial effects ---#
## spat.sig <- data.frame(n.nodes.in.d=length(in.d),
##                        expected.out=0.05*length(in.d),
##   observed.y.out=sum(out.yz$summary.random$y.field[in.d,4]>0 |
##                        out.yz$summary.random$y.field[in.d,6]<0),
##   observed.z.out=sum(out.yz$summary.random$z.field[in.d,4]>0 |
##                        out.yz$summary.random$z.field[in.d,6]<0))
## #write.csv(spat.sig, file.path(o.d, glue("{lr.n}-spat-sig.csv")), row.names = F)

## out.yz.summ <- data.frame(lr.n = lr.n,
##                           side = side,
##                           var = c("local.dic",
##                                   "family.dic",
##                                   "yPos-obs",
##                                   "yPos-pred",
##                                   "z-obs",
##                                   "z-pred",
##                                   colnames(spat.sig)
##                                   ),
##                           val = c(tapply(out.yz$dic$local.dic, out.yz$dic$family, sum),
##                                   # the tapply returns both local and family dic
##                                   mean(lr.d[which(lr.d > 0)]),
##                                   mean(exp.y),
##                                   mean(z, na.rm = T),
##                                   mean(exp.z),
##                                   as.numeric(spat.sig[1, ]))
##                          )

## write.csv(out.yz.summ, file.path(o.d, glue("{lr.n}-{side}-output-summary.csv")), row.names = F)


## extract the spatial effects
y.field <- inla.spde2.result(out.yz, "y.field", spde)
z.field <- inla.spde2.result(out.yz, "z.field", spde)

## png(file.path(o.d, glue("{lr.n}-{side}-spatial-effects-marginals.png")),
##     width = 8, height = 8, units = "in", res = 300)
## par(mfrow=c(2,2))
## plot.default(y.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
## plot.default(y.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
## plot.default(z.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
## plot.default(z.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')
## dev.off()

#### plot spatial random fields

## make projection grid

# predict resolution (number of pixels in x and y)
xy.pred <- make.pred.grid(xy.obs = xy.obs,
                          retain.shape = T,
                          res = pred.grid.res)
if(j.i  == 1){
  write.csv(xy.pred, file.path(o.d, glue("xy-pred-locs-{side}.csv")), row.names = F)
}

A.pred <- inla.spde.make.A(
  mesh = mesh,
  loc = xy.pred)

## rf.grid <- list(xi.mean = fm_evaluate(mesh, field = out.y$summary.random$y.field$mean, loc = xy.pred),
##                 xi.sd = fm_evaluate(mesh, field = out.y$summary.random$y.field$sd, loc = xy.pred),
##                 obs = lr.d,
##                 u.mean = fm_evaluate(mesh, field = out.yz$summary.random$z.field$mean, loc = xy.pred),
##                 u.sd = fm_evaluate(mesh, field = out.yz$summary.random$z.field$sd, loc = xy.pred),
##                 occ.obs = z)

## xy.in <- inout(xy.pred, boundary.coords)
## for (i in c(1:2, 4:5)) rf.grid[[i]][!xy.in] <- NA


## png(file.path(o.d, glue("{lr.n}-{side}-spatial-fields.png")),
##     width = 12, height = 8, units = "in", res = 500)
## par(mfrow=c(2,3), mar=c(0,0,0,0))
## for (i in 1:6) {
##   if(grepl("obs", names(rf.grid)[i])){
##     quilt.plot(xy.obs, rf.grid[[i]],
##                nx = 60, ny = 50,
##                nlevel=100, col=viridis(100),
##                legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
##                xaxt = 'n', yaxt = 'n')
##   }else{
##     quilt.plot(xy.pred, rf.grid[[i]],
##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
##                nlevel=101, col=viridis(100),
##                legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
##                xaxt = 'n', yaxt = 'n')
##   }
##   lines(boundary.coords)
## }
## dev.off()

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
#prd.y.proj <- apply(prd.y, 2,
#                    function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
#prd.z.proj <- apply(prd.z, 2,
#                    function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})

mean.y.cond.z.proj <- fm_evaluate(mesh, field=rowMeans(prd.y.cond.z), loc = xy.pred)

line.coords <- data.frame(x = seq(21, 115, length = 250), y = rep(-35, 250)) |> as.matrix()

prd.z.line.proj <- apply(prd.y.cond.z, 2,
                         function(s){fm_evaluate(mesh, field=s, loc = line.coords)})

line.res <- data.frame(feat = lr.n,
                       x = line.coords[, 1],
                       y = apply(prd.z.line.proj, 1, quantile, probs = .5),
                       l = apply(prd.z.line.proj, 1, quantile, probs = .25),
                       u = apply(prd.z.line.proj, 1, quantile, probs = .75)
                       )

write.csv(line.res, file.path(o.d, glue("{lr.n}-{side}-line-summ.csv")), row.names = F)
## plot(line.res$x, line.res$y)
## lines(line.res$x, line.res$u)
## lines(line.res$x, line.res$l)

top.line.coords <- bot.line.coords <- line.coords
top.line.coords[, 2] <- top.line.coords[, 2] + .5
bot.line.coords[, 2] <- bot.line.coords[, 2] - .5
line.bbox <- rbind(top.line.coords,
                   bot.line.coords[nrow(bot.line.coords):1, ],
                   top.line.coords[1, ])

png(file.path(o.d, glue("{lr.n}-{side}-boxed-mean.png")),
    width = 8, height = 8, units = "in", res = 300)
par(mar = c(0, 0, 2, 0))
quilt.plot(xy.pred, mean.y.cond.z.proj,
           nx = pred.grid.res * .9, ny = pred.grid.res * .9,
           nlevel=101, col=viridis(100), axes = FALSE,
           legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
           legend.args=list(text="", line=0),
           xaxt = 'n', yaxt = 'n', main = lr.n)
lines(line.bbox, col = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")[which((sort(lr.to.save.draws[c(5, 6, 7, 10)])) == lr.n)], lwd = 5)
dev.off()



# pull back in relevant csvs and stitch together ggplot

## rm(pd)
## for(lr.n in lr.to.save.draws[c( 5, 6, 7, 10)]){

##   if(exists("pd")){
##     new.pd <- read.csv(file.path(o.d, glue("{lr.n}-{side}-line-summ.csv")))
##     pd <- rbind(pd, new.pd)
##   }else{
##     pd <- read.csv(file.path(o.d, glue("{lr.n}-{side}-line-summ.csv")))
##   }

## }

## lp <- ggplot(pd, aes(x = x, y = y, group = feat, col = feat)) +
##   geom_line() +
##   geom_ribbon(aes(ymin = l, ymax = u, fill = feat), alpha = 0.25) +
##   xlab("") + ylab("NICHES Cell Signaling")
## #  scale_fill_manual(values = rainbow(10)[c(4, 5, 6, 7, 10)]) +
## #    scale_color_manual(values = rainbow(10)[c(4, 5, 6, 7, 10)])
## lp

## ggsave(file.path(o.d, glue("boxed-multi-lines.png")),lp,
##     width = 15, height = 5, units = "in")


## prd.y.cond.z.proj <- apply(prd.y.cond.z, 2,
##                            function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})


## #### make plots showing different posterior draws



## for(i in 1:20){
## png(file.path(o.d, glue("{lr.n}-{side}-response-occurrence-fields-draws-{i}.png")),
##     width = 16, height = 6, units = "in", res = 500)

## par(mfrow=c(1,3), mar=c(0,0,0,0))

## zl <- range(c(mean.y.cond.z.proj,
##               lr.d,
##               prd.y.cond.z.proj[, 1:20]))

##     quilt.plot(xy.pred, mean.y.cond.z.proj,
##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
##                nlevel=101, col=viridis(100), axes = FALSE,
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##                legend.args=list(text="", line=0),
##                xaxt = 'n', yaxt = 'n', zlim = zl)#, main = titles[i])
## text(25, -15, glue("{lr.n} Mean "), cex = 2)
## lines(boundary.coords)

## xy.obs.plot <- xy.obs


## bubblePlot(xy.obs, lr.d,
##        #        nx = 60, ny = 50,
##                nlevel=100, col=viridis(100),
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##            legend.args=list(text="", line=0), size = 1.5,
##            axes = F, lab.breaks = NULL, xlab = "", ylab = "",
##             zlim = zl)#, main = titles[i])
## text(25, -13, glue("{lr.n} NICHES"), cex = 2)
##   lines(boundary.coords)

## quilt.plot(xy.pred, prd.y.cond.z.proj[, i],
##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
##                nlevel=101, col=viridis(100), axes = FALSE,
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##                legend.args=list(text="", line=0),
##                xaxt = 'n', yaxt = 'n', zlim = zl)#, main = titles[i])
## text(25, -15, glue("{lr.n} Draw {i}"), cex = 2)
##   lines(boundary.coords)
## dev.off()

## }

## png(file.path(o.d, glue("{lr.n}-{side}-raw-data-dot-sizes.png")),
##     width = 16, height = 6, units = "in", res = 500)

## par(mfrow=c(1,3), mar=c(0,0,0,0))

## bubblePlot(xy.obs, lr.d,
##        #        nx = 60, ny = 50,
##                nlevel=100, col=viridis(100),
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##            legend.args=list(text="", line=0), size = 1,
##            axes = F, lab.breaks = NULL, xlab = "", ylab = "",
##             zlim = zl)#, main = titles[i])
## text(25, -13, glue("{lr.n} NICHES"), cex = 2)
##   lines(boundary.coords)

## bubblePlot(xy.obs, lr.d,
##        #        nx = 60, ny = 50,
##                nlevel=100, col=viridis(100),
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##            legend.args=list(text="", line=0), size = 1.5,
##            axes = F, lab.breaks = NULL, xlab = "", ylab = "",
##             zlim = zl)#, main = titles[i])
## text(25, -13, glue("{lr.n} NICHES"), cex = 2)
##   lines(boundary.coords)

## bubblePlot(xy.obs, lr.d,
##        #        nx = 60, ny = 50,
##                nlevel=100, col=viridis(100),
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##            legend.args=list(text="", line=0), size = 2,
##            axes = F, lab.breaks = NULL, xlab = "", ylab = "",
##             zlim = zl)#, main = titles[i])
## text(25, -13, glue("{lr.n} NICHES"), cex = 2)
##   lines(boundary.coords)


## dev.off()


## # average over areas

## png(file.path(o.d, glue("{lr.n}-{side}-areal-boxes.png")),
##     width = 6, height = 6, units = "in", res = 500)
## quilt.plot(xy.pred, mean.y.cond.z.proj,
##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
##                nlevel=101, col=viridis(100), axes = FALSE,
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##                legend.args=list(text="", line=0),
##                xaxt = 'n', yaxt = 'n', zlim = zl)#, main = titles[i])
## box1 <- matrix(c(68, -36, 72, -36, 72, -34, 68, -34, 68, -36), ncol = 2, byrow = T)
## box2 <- matrix(c(95, -31, 99, -31, 99, -29, 95, -29, 95, -31), ncol = 2, byrow = T)
## lines(box1, col = "deeppink")
## lines(box2, col = "cyan")
## dev.off()


## in.b1 <- which(inout(xy.pred[, 1:2], box1))
## in.b2 <- which(inout(xy.pred[, 1:2], box2))


## b1.avg <- colMeans(prd.y.cond.z.proj[in.b1, ])
## b2.avg <- colMeans(prd.y.cond.z.proj[in.b2, ])

## df <- data.frame(avg = c(b1.avg, b2.avg),
##                  box = c(rep("box1", 1000),
##                          rep("box2", 1000)))

## # Use semi-transparent fill
## p<-ggplot(df, aes(x=avg, fill=box, color=box)) +
##   geom_histogram(position="identity", alpha=0.5) + ggtitle("Comparison of the Distribution of Average Communication within two Areas") + theme(legend.position="none")
## p
## ggsave(p, filename = file.path(o.d, glue("{lr.n}-{side}-areal-boxes-distributions.png")),
##        height = 6, width = 8)



## prd.y.cond.z.line.proj <- fm_evaluate(mesh, field=rowMeans(prd.y.cond.z), loc = line.coords)


## # always write the means
## write.csv(as.matrix(rowMeans(prd.y.proj), ncol = 1),
##           file.path(o.d, glue("{lr.n}-{side}-resp-est.csv")), row.names = F)
## write.csv(as.matrix(rowMeans(prd.y.cond.z.proj), ncol = 1),
##           file.path(o.d, glue("{lr.n}-{side}-cond-resp-est.csv")), row.names = F)

## write.csv(as.matrix(rowMeans(prd.y.cond.z.line.proj), ncol = 1),
##           file.path(o.d, glue("{lr.n}-{side}-cond-resp-line-est.csv")), row.names = F)
## write.csv(as.matrix(apply(prd.y.cond.z.proj, 1, sd), ncol = 1),
##           file.path(o.d, glue("{lr.n}-{side}-cond-resp-line-sd.csv")), row.names = F)


## # TODO always write the SD
## # TODO ASK SAM: better to save as separate file, or as second column in one file?
## #write.csv(as.matrix(rowMeans(prd.y.proj), ncol = 1),
## #          file.path(o.d, glue("{lr.n}-{side}-resp-est.csv")), row.names = F)
## #write.csv(as.matrix(rowMeans(prd.y.cond.z.proj), ncol = 1),
## #          file.path(o.d, glue("{lr.n}-{side}-cond-resp-est.csv")), row.names = F)
## # only save some of the complete draws -- these are big files (.5Gb)
## ## if(lr.n %in% lr.to.save.draws){
## ##   write.csv(prd.y.cond.z.proj,
## ##             file.path(o.d, glue("{lr.n}-{side}-cond-resp-draws.csv")), row.names = F)
## ##   write.csv(prd.y.proj,
## ##             file.path(o.d, glue("{lr.n}-{side}-resp-draws.csv")), row.names = F)
## ## }


## prd <- list(y.cond.z.mean=fm_evaluate(mesh, field=rowMeans(prd.y.cond.z), loc = xy.pred))
## prd$y.cond.z.sd <- fm_evaluate(mesh, field=apply(prd.y.cond.z, 1, sd), loc = xy.pred)
## prd$cond.resp.obs<- lr.d

## prd$y.mean <- fm_evaluate(mesh, field=rowMeans(prd.y), loc = xy.pred)
## prd$y.sd <- fm_evaluate(mesh, field=apply(prd.y, 1, sd), loc = xy.pred)
## prd$resp.obs <- lr.d

## prd$z.mean <- fm_evaluate(mesh, field=rowMeans(prd.z), loc = xy.pred)
## prd$z.sd <- fm_evaluate(mesh, field=apply(prd.z, 1, sd), loc = xy.pred)
## prd$occ.obs <- z
## for (j in c(1:2, 4:5, 7:8)) prd[[j]][!xy.in] <- NA

## titles <- c("Y|Z Estimates", "Y|Z sd", "Y|Z obs",
##             "Y Estimates", "Y sd", "Y obs",
##             "Z Estimates", "Z sd", "Z obs")

## # *** Code for Figure 8.14
## png(file.path(o.d, glue("{lr.n}-{side}-response-occurrence-fields.png")),
##     width = 12, height = 12, units = "in", res = 500)
## par(mfcol=c(3,3), mar=c(0,0,0,0))
## for (i in 1:length(prd)) {

##   if(i %in% c(1, 3, 4, 6)){
##     zl <- c(0, max(prd$resp.obs))
##   }else if(i %in% c(7, 9)){
##     zl <- c(0, 1)
##   }else{
##     zl <- range(prd[[i]])
##   }



##   if(grepl("obs", names(prd)[i])){
##     ## TODO - move to Seurat plotting for raw data.
##     ## I think I need to move what I want to plot into brain to do this
##     ## starter code:
##     ## DefaultAssay(brain)<- "NeighborhoodToCell_ALRA"
##     ## SpatialFeaturePlot(brain, features = lr.n, slot = "data") +
##     ##   scale_fill_viridis_c(limits = zl)
##     quilt.plot(xy.obs, prd[[i]],
##                nx = 60, ny = 50,
##                nlevel=100, col=viridis(100),
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##                legend.args=list(text="", line=0),
##                xaxt = 'n', yaxt = 'n', zlim = zl)#, main = titles[i])
##     text(25, -15, titles[i], cex = 2)
##   }else{
##     quilt.plot(xy.pred, prd[[i]],
##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
##                nlevel=101, col=viridis(100), axes = FALSE,
##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
##                legend.args=list(text="", line=0),
##                xaxt = 'n', yaxt = 'n', zlim = zl)#, main = titles[i])
##     text(25, -15, titles[i], cex = 2)
##   }
##   lines(boundary.coords)
## }
## dev.off()


## ## png(file.path(o.d, "test.png"), width = 5, height = 5, units = "in", res = 300)
## ## DefaultAssay(brain)<- "NeighborhoodToCell_ALRA"
## ## SpatialFeaturePlot(brain, features = lr.n, slot = "data") +
## ##   scale_fill_viridis_c(limits = c(0.2, 1))
## ## dev.off()
