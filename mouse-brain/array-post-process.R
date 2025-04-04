## combine all the mean estimates from an output dir

o.d <- "~/a-s-omics/data-outputs/2024-03-14/"

a.f.p <- list.files(o.d, full.names = T)
a.f   <- list.files(o.d)

# subset to mean estimates
m.f.p <- grep("est.csv", a.f.p, value = T)
m.f   <- grep("est.csv", a.f, value = T)

# make a big storage matrix, grabbing dims from one file
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
