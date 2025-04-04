require(ReacTran)

# a 1d finite difference grid example from the pkg vignette
test.grid.1d <- setup.grid.1D(L = 12000, , N = 240, dx.1 = 50)
plot(test.grid.1d)

n = 100
dy = dx = 100 / n
Dy = Dx = 5
r = -0.02
Bc = 300
irr <- 20
vx <- 1
y <- matrix(nrow = n, ncol = n, 0)


Diff2D <- function (t, y, parms, N) {

  CONC <- matrix(nrow = N, ncol = N, y)
  # Transport
  Tran <-tran.2D(CONC, D.x = Dx, D.y = Dy, C.y.down = Bc,
                 dx = dx, dy = dy, v.x = vx)

  # transport  reaction
  dCONC <- Tran$dC + r*CONC

  # Bioirrigation in a central spot
  mid <- N/2
  dCONC[mid, mid] <- dCONC[mid, mid] + irr*(Bc - CONC[mid, mid])

  return (list(dCONC))
}

times <- seq(0, 100, 5)
print(system.time(
  out2 <- ode.2D(func = Diff2D, y = as.vector(y), times = times, N = n,
                 parms = NULL, lrw = 10000000, dimens = c(n, n))
))

for(i in seq(0, 95, 5)){
  image(out2, method = "filled.contour", main = "")
}
