## to save time, we make some objects in advance.
## in particular, we make:

# 1) we make the nonconvex hull, mesh

require(sf)
require(splancs)
require(fmesher)
require(concaveman)
require(ggplot2)

mesh.list <- list(NULL)
gg.list <- list(NULL)
for(time.i in 1:4){ # timepoints

  message(glue("on time {time.i}"))

  xy.obs <- cbind( niches.alra.list[[time.i]]$NeighborhoodToCell@meta.data$x |> as.numeric(),
                  niches.alra.list[[time.i]]$NeighborhoodToCell@meta.data$y |> as.numeric())

  coords <- as.matrix(xy.obs)

  # concaveman returns coords in CLOCKWISE fashion
  # INLA interprets interior to be defined relative to COUTERCLOCKWISE order!
  boundary.coords <- concaveman(coords, concavity = 1.6, length_threshold = 0)
  domain <- boundary.coords[nrow(boundary.coords):1, ] |> fm_segm()
   plot(domain)

  # tune parameters per time point

  max.edge = diff(range(coords[,1]))/(3*5)
  bound.outer = diff(range(coords[,1]))/3

  m.e.inner <- c(.35, .25, .2, .1)
  m.e.cutoff <- c(10, 10, 10, 12)


  # make mesh

  mesh <- fm_mesh_2d(boundary = list(fm_segm(c(domain))),
                      max.edge = c(m.e.inner[time.i], 2) * max.edge,
                      offset = c(max.edge, bound.outer / 1.5),
                      cutoff = max.edge / m.e.cutoff[time.i],
                     min.angle = 21)

  mesh.list[[time.i]] <- mesh

  # plot mesh

  gg.list[[time.i]] <- ggplot() + inlabru::gg(mesh) + ggtitle(glue("time: {time.i} -- mesh nodes: {mesh$n} -- data obs: {nrow(coords)}"))
  ggsave(gg.list[[time.i]], filename = file.path(i.d, glue("mesh-time-{time.i}.png")))


  # check which mesh vertices are inside the domain
  ds1 <- as.points(boundary.coords)
  ds1.poly <- ds1[chull(ds1),]
  in.d <- which(splancs::inout(mesh$loc[, 1:2], ds1))

  save("coords", "domain", "mesh", "boundary.coords", "in.d",
       file = file.path(i.d, glue("pre-fit-obj-{time.i}.Rdata")))

}


mesh.list2 <- list(NULL)
for(i in 1:3){
  mesh.list2[[i]] <- fm_as_sfc(mesh.list[[i]])
}

mesh.list2[[2]] <- mesh.list2[[2]] + c(75, 0)
mesh.list2[[3]] <- mesh.list2[[3]] + c(167, 0)
mesh.sfc <- c(mesh.list2[[1]], mesh.list2[[2]], mesh.list2[[3]])

png(file.path(i.d, "meshes-same-axes.png"),
    width = 20, height = 12, units = "in", res = 300)
plot(mesh.sfc)
dev.off()



## comb.mesh <- fm_mesh_2d(loc = mesh.sfc, boundary = boundary)
## ggplot() + inlabru::gg(comb.mesh)
## #plot(comb.mesh)

## all.mesh <- c(mesh.list[[1]], mesh.list[[2]])
## ggplot() + inlabru::gg()

## plot 3
gg.mesh.1.3 <- ggpubr::ggarrange(plotlist = gg.list, ncol = 3, nrow = 1)
## print(gg.mesh.1.3)
ggsave(gg.mesh.1.3, filename = file.path(i.d, glue("meshes-times-1-3.png")),
    width = 22, height = 12, units = 'in')

## ## plot 4
## gg.all.mesh <- ggpubr::ggarrange(plotlist = gg.list, ncol = 2, nrow = 2)
## print(gg.all.mesh)
## ggsave(gg.all.mesh, filename = file.path(i.d, glue("meshes-times-1-4.png")),
##     width = 16, height = 16, units = 'in')
