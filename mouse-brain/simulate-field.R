require(fields)

pred.xy <- data.frame(1:500, 1:500)

obj<- circulantEmbeddingSetup( pred.xy, Covariance="Matern", range = 2, smoothness = .5, phi = 1)
# range - dist before spatial correlation decays to 10%, >0
# smoothness - mean square differentiability, >0
# phi - amplitude, >0
set.seed( 223)
# simulate 1 realization
look <- circulantEmbedding( obj)
look2 <- circulantEmbedding( obj)
# take a look at first two
set.panel(3, 3)

for(i in 1:9){

  # simulate 1 realization
  look <- circulantEmbedding( obj)

  # plot
  image.plot( sort(unique(pred.xy[, 1])), sort(unique(pred.xy[, 2])), look)
  title("simulated gaussian fields")
}
