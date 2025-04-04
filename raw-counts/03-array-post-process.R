## plot the preditions
## TODO: combine all the mean estimates from an output dir

setwd("~/Dropbox/genetics/a-s-omics/")
i.d <- file.path(getwd(), "/data-inputs/mouse-embryo-raw")
o.d <- "~/Dropbox/genetics/a-s-omics//data-outputs/mouse-embryo-raw/2025-01-27" # initial output location

require(data.table)
require(deSolve)
require(ReacTran)
require(glue)
require(fields)

# make plots?
plot.preds <- T

# get files in output dir
a.f.p <- list.files(file.path(o.d, "prediction-objects"), full.names = T)
a.f   <- list.files(file.path(o.d, "prediction-objects"))

# subset to pred objects
p.f.p <- grep("rds", a.f.p, value = T)
p.f   <- grep("rds", a.f, value = T)

# get list of all LR mechanism features
feats <- grep(".rds", a.f, value = T) |>
  strsplit("\\.") |>
  lapply(function(x){x[[1]]}) |>
  unlist() |>
  strsplit("-") |>
  lapply(function(x){paste0(tail(x, 2), collapse = '-')}) |>
  unlist() |> unique() |> sort()

# and get the individual L and R names
lr.feats <- strsplit(feats, "-") |> unlist() |> unique() |> sort()

####################################################
## Collate all putputs into one convenient object ##
####################################################

# TODO

####################################
## plot output prediction objects ##
####################################

if(plot.preds){
  # loop through the mechanism preds and process
  dir.create(file.path(o.d, "convolution-plots"))
  for(pred.fp in p.f.p){
    cat(glue("on {which(p.f.p==pred.fp)} of {length(p.f.p)}\n"))
    cat("\n")
    preds <- readRDS(pred.fp)
    pred.t <- preds$total
    pred.l <- preds$feat.l
    pred.r <- preds$feat.r

    m.n <- pred.fp |>
      strsplit("\\.") |>
      lapply(function(x){x[[1]]}) |>
      unlist() |>
      strsplit("-") |>
      lapply(function(x){paste0(tail(x, 2), collapse = '-')}) |>
      unlist() |> unique() |> sort()
    l.n <- strsplit(m.n, "-")[[1]][1]
    r.n <- strsplit(m.n, "-")[[1]][2]

    ## plot

    qp.res.x <- pred.t[, length(unique(x))]
    qp.res.y <- pred.t[, length(unique(y))]

    png(file.path(o.d, "convolution-plots", glue('data-vs-model-pois-pois-CONVOLVED-{l.n}-{r.n}.png')),
        width = (qp.res.x / qp.res.y) * 13 + 8, height = 13, units = 'in', res = 300)
    par(mfrow = c(3, 3),
        mai = c(.62, 0.82, .62, 1.22))
    fields.style();quilt.plot(pred.l[, x],
                              pred.l[, y],
                              pred.l[, mean],
                              main = glue('Estimated {l.n} Counts' ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.l[, x],
                              pred.l[, y],
                              pred.l[, mean] / pred.t[, mean],
                              main = glue('Estimated {l.n} Counts per Total Counts' ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.t[, x],
                              pred.t[, y],
                              pred.t[, mean],
                              main = glue('Estimated Total Counts'),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)

    fields.style();quilt.plot(pred.r[, x],
                              pred.r[, y],
                              pred.r[, mean],
                              main = glue('Estimated {r.n} Counts' ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.r[, x],
                              pred.r[, y],
                              pred.r[, mean] / pred.t[, mean],
                              main = glue('Estimated {r.n} Counts per Total Counts' ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.t[, x],
                              pred.t[, y],
                              pred.t[, mean],
                              main = glue('Estimated Total Counts'),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)


    fields.style();quilt.plot(pred.l[, x],
                              pred.l[, y],
                              pred.l[, mean] * pred.r[, mean],
                              main = glue('{l.n}-{r.n} Convolved Counts' ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.l[, x],
                              pred.l[, y],
                              pred.l[, mean] * pred.r[, mean] /(pred.t[, mean] ^ 2),
                              main = glue('{l.n}-{r.n} Convolved Density' ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    fields.style();quilt.plot(pred.t[, x],
                              pred.t[, y],
                              pred.t[, mean] ^ 2,
                              main = expression('(Estimated Total Counts)' ^ 2),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
    dev.off()
  }
}
