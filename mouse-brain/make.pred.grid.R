# makes a grid by bounding the input locs in a square and making a res X res grid
make.pred.grid <- function(xy.obs, res = grid.res, retain.shape = T){

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

    xy.keep <- rep(FALSE, nrow(xy.pred))

    for(xx in sort(unique(xy.obs[, 1]))){

      # find nearest new columns
      xx.dist <- abs(x.grid - xx)
      xxp1.dist <- abs(x.grid - (xx + 1))
      xxm1.dist <- abs(x.grid - (xx - 1))
      new.xx <- x.grid[which(xx.dist < xxp1.dist & xx.dist < xxm1.dist)]

      # get y bounds
      r.y <- range(xy.obs[xy.obs[, 1] == xx, 2])

      xy.keep[which(xy.pred[, 1] %in% new.xx &
                      xy.pred[, 2] >= r.y[1] & xy.pred[, 2] <= r.y[2])] <- T
    }
    xy.pred <- xy.pred[xy.keep, ]
  }

  return(as.matrix(xy.pred))

}
