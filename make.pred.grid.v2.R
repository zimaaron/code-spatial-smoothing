# makes a grid by bounding the input locs in a square and making a res X res grid
make.pred.grid <- function(xy.obs, res = grid.res, retain.shape = T, boundary.poly){

  # box
  min.x <- min(xy.obs[,1])
  max.x <- max(xy.obs[,1])

  min.y <- min(xy.obs[,2])
  max.y <- max(xy.obs[,2])

  # new axes
  y.grid <- seq(min.y, max.y, length = res)
  x.grid <- seq(min.x, max.x, length = res)

  # grid
  xy.pred <- expand.grid(x.grid, y.grid)
  colnames(xy.pred) <- c("x", "y")

  if(retain.shape){

    # then only predict inside initial domain

    xy.keep <- which(splancs::inout(as.matrix(xy.pred), boundary.poly))
    xy.pred <- xy.pred[xy.keep, ]

  }

  return(as.matrix(xy.pred))

}
