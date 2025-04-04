## combine all the mean estimates from an output dir

o.d <- "~/a-s-omics/data-outputs/mouse-embryo/2024-08-28-100test"

a.f.p <- list.files(o.d, full.names = T)
a.f   <- list.files(o.d)

# get list of all features
feats <- grep("—", a.f, value = T) |>
  strsplit("-") |>
  lapply(function(x){x[[1]]}) |>
  unlist() |> unique() |> sort()

# get list of all the time points run
times <- grep("—", a.f, value = T) |>
  strsplit("-") |>
  lapply(function(x){x[[2]]}) |>
  unlist() |> unique() |> sort()

# subset to mean estimates
m.f.p <- grep("est.csv", a.f.p, value = T)
m.f   <- grep("est.csv", a.f, value = T)

# subset to sd
s.f.p <- grep("sd.csv", a.f.p, value = T)
s.f   <- grep("sd.csv", a.f, value = T)


# make a big storage matrix for all of the results for each of the times,
# grabbing dims from one of each of the pred loc files
a.m.list <- list(NULL) # means
a.s.list <- list(NULL) # sds
for(i in as.numeric(times)){

}

d.t <- read.csv(m.f.p[1])

a.r <- data.frame(matrix(nrow = nrow(d.t), ncol = length(m.f.p)))

# loop through files, read in, combine
for(i in 1:length(m.f.p)){
  if(i %% 250 == 0){ message(glue("on iter {i}"))}
  d.t <- read.csv(m.f.p[i])
  a.r[, i] <- d.t[, 1]
  colnames(a.r)[i] <- gsub("-est.csv", "", m.f[i])

  ## TODO combine the DIC results
  # glue("{lr.n}-DIC-comp.csv")

  ## TODO combine the means of preds/obs results
  # glue("{lr.n}obs-est-means-comp.csv")
}

# save output
write.csv(a.r, file.path(o.d, "0000-all-smooth-ests.csv"), row.names = F)
