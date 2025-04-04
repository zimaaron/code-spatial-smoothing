
## this file contains functions to perform spatial kriging

## notes:
#
# 0) input takes the form of a matrix with x,y pairs, and
#    a second matrix with data1, data2, ... columns corresponding to the xy mat
#    output has the same form with different/denser x,y
#
# 1) the user must supply the parameters for the covariance function
# --> IMPORTANT --> this means these are not estimated from the data
#
# 2) input

require(data.table)
# require(glue)
require(fields) # require(rSPDE)

## matern cov from fields
## Matern(d , range = 1,alpha=1/range, smoothness = 0.5,
##        nu= smoothness, phi=1.0)
## x <- seq(0,10,by=.01)
## plot(x = x, y = Matern(d = x, range = 1, smoothness=.5, phi = 1))
## lines(x = x, y = Matern(d = x, range = 2, smoothness=.5, phi = 1), lty = 1)
## lines(x = x, y = Matern(d = x, range = 2, smoothness=2, phi = 1), lty = 2)
## lines(x = x, y = Matern(d = x, range = 2, smoothness=2, phi = 2), lty = 3)

# makes a grid by bounding the input locs in a square and making a res X res grid
make.pred.grid <- function(xy.obs, res = grid.res, retain.shape = T){

  # box
  min.x <- min(xy.obs$x)
  max.x <- max(xy.obs$x)

  min.y <- min(xy.obs$y)
  max.y <- max(xy.obs$y)

  # new axes
  y.grid <- seq(min.y, max.y, length = res)
  x.grid <- seq(min.x, max.x, length = res)

  # grid
  xy.pred <- expand.grid(x.grid, y.grid)
  colnames(xy.pred) <- c("x", "y")

  if(retain.shape){

    xy.pred$keep <- FALSE

    for(xx in sort(unique(xy.obs$x))){

      # find nearest new columns
      xx.dist <- abs(x.grid - xx)
      xxp1.dist <- abs(x.grid - (xx + 1))
      xxm1.dist <- abs(x.grid - (xx - 1))
      new.xx <- x.grid[which(xx.dist < xxp1.dist & xx.dist < xxm1.dist)]

      # get y bounds
      r.y <- range(xy.obs[xy.obs$x == xx, ]$y)

      xy.pred$keep[which(xy.pred$x %in% new.xx &
                     xy.pred$y >= r.y[1] & xy.pred$y <= r.y[2])] <- T
    }
    xy.pred <- subset(xy.pred, keep == T)
    xy.pred$keep <- NULL

  }

  return(xy.pred)

}

# fit smoothing
fit.smooth.field <- function(xy.data, data.vec, mat.range, mat.smooth){

  Krig(x = xy.data,
       Y = data.vec,
       Covariance = "Matern", theta = mat.range, smoothness = mat.smooth)

}

# predict smooth
pred.smooth.field <- function(smooth.fit, xy.pred){

  predict.Krig(smooth.fit, x = xy.pred)

}


### test

load("./code-spatial-smoothing/LR.output.2023-11-08.Robj")
# str(LR.output)

xy.obs <- LR.output$brain1$location.info
lig.obs <- LR.output$brain1$ligand.info
lig.obs <- lig.obs["Fgf1", ]
rec.obs <- LR.output$brain1$receptor.info
lig.obs <- rec.obs["Fgfr2", ]

# initialize grid
grid.res <- 200
xy.pred <- make.pred.grid(xy.obs = xy.obs, retain.shape = F, res = grid.res)

# ligand smoothing
lig.smooth.fit <- fit.smooth.field(xy.data = xy.obs,
                                   data.vec = lig.obs, mat.range = 5, mat.smooth = .5)
lig.smooth <- pred.smooth.field(lig.smooth.fit, xy.pred = xy.pred)
quilt.plot(y=-xy.pred$x, x=xy.pred$y, z=lig.smooth,nx=grid.res, ny=grid.res)


# Receptor smoothing
rec.smooth.fit <- fit.smooth.field(xy.data = xy.obs,
                                   data.vec = lig.obs, mat.range = 1, mat.smooth = 5)
rec.smooth <- pred.smooth.field(rec.smooth.fit, xy.pred = xy.pred)
quilt.plot(y=-xy.pred$x, x=xy.pred$y, z=rec.smooth,nx=grid.res, ny=grid.res )

quilt.plot(y=-xy.pred$x, x=xy.pred$y, z=rec.smooth*lig.smooth,nx=grid.res, ny=grid.res )


# Niches data smoothing
niche.smooth.fit <- fit.smooth.field(xy.data = xy.obs,
                                     data.vec = as.numeric(fgf1.fgfr2.matrix),
                                     mat.range = 0.5, mat.smooth = 5)
niche.smooth <- pred.smooth.field(niche.smooth.fit, xy.pred = xy.pred)

png(file = 'niches.smooth.1000test.png',width = 12,height = 10,units = 'in',res=1000)
quilt.plot(y=-xy.pred$x, x=xy.pred$y, z=niche.smooth,nx=grid.res, ny=grid.res,col = Seurat:::SpatialColors(100))
dev.off()
