## aoz dec 2024
##########################################
## set params that can't be set earlier ##
##########################################

# predict resolution (number of pixels in x and y)
qp.res.x <- run.dat[, length(unique(x))]
qp.res.y <- run.dat[, length(unique(y))]

dir.create(file.path(o.d, "prediction-objects")) # dir for model fit objects
## dir.create(file.path(o.d, "fitted-models")) # dir for prediction outputs


############################
## setup spde matern
############################

## prior matern params
## c(a, b, c, d), where
## P(sp.range < a) = b
## P(sp.sigma > c) = d
#matern.pri.total <- c(2000, .95, 1, .05) ## a, b, c, d
matern.pri.feat.count <- c(500, .95, 1, .05) ## a, b, c, d
matern.pri.feat.present <- c(1000, .95, 10, .05) ## a, b, c, d

## matern.remain <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
##                                      constr = TRUE, # integrate-to-zero constraint
##                                      prior.range = matern.pri.total[1:2],
##                                      prior.sigma = matern.pri.total[3:4])
matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                         constr = TRUE, # integrate-to-zero constraint
                                         prior.range = matern.pri.feat.count[1:2],
                                         prior.sigma = matern.pri.feat.count[3:4])
## TODO, tune this prior? it's hard with the 0s and 1s areas in space
## and the logit transform
matern.feat.present <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.present[1:2],
                                           prior.sigma = matern.pri.feat.present[3:4])


#########################
## run INLA bru model(s)
#########################

# We want to obtain CPO data from the estimations
bru_options_set(control.compute = list(cpo = TRUE))
n.samp <- 1000


##~~~~~~~~~~~~~~~~~
## poisson model
##~~~~~~~~~~~~~~~~~
if(poi.mod){

  matern.pri.feat.count <- c(500, .95, .1, .05) ## a, b, c, d
  matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.count[1:2],
                                           prior.sigma = matern.pri.feat.count[3:4])

  comps <- ~
    feat.count.int(1) +
    feat.count.field(cbind(x, y), model = matern.feat.count)

  poi.lik <- bru_obs(
    family = "poisson",
    data = run.dat,
    formula = feat.count ~ feat.count.int + feat.count.field,
    E = total.count
  )

  try(
  {
    cat('\n')
    for(i in 1){print(glue('{lr.n}: fitting poisson model'))}
    cat('\n')

    fit.poi <- bru(comps,
                   poi.lik,
                   # options
                   options = list(bru_verbose = 4,
                                  bru_max_iter = 1)
                   )

  }, silent = T)

  if(exists("fit.poi")){
    try(
    {
      cat('\n')
      for(i in 1){print(glue('{lr.n}: predicting from fitted poisson model'))}
      cat('\n')

      pred.poi <- predict(
        fit.poi, run.dat,
        ~ {
          lambda <- exp( feat.count.int + feat.count.field )
          expect <- lambda * total.count
          list(
            lambda = lambda,
            expect = expect,
            obs_prob = dpois(feat.count, expect)
          )
        },
        n.samples = n.samp
      )
    }, silent = T)
  }

  if(exists("fit.poi")){
    poi_pit <- fit.poi$cpo$pit * c(NA_real_, 1)[1 + (run.dat$feat.count > 0)]
  }else{
    poi_pit <- NA
  }

  if(exists("pred.poi")){
    # For Poisson, the posterior conditional variance is equal to
    # the posterior conditional mean, so no need to compute it separately.
    expect_poi <- pred.poi$expect
    expect_poi$pred_var <- expect_poi$mean + expect_poi$sd^2
    expect_poi$log_score <- -log(pred.poi$obs_prob$mean)
  }else{
    expect_poi$expect <- run.dat
    expect_poi$mean <- NA
    expect_poi$median <- NA
    expect_poi$pred_var <- NA
    expect_poi$log_score <- NA
  }
}

##~~~~~~~~~~~~~~~~~
## zip model
##~~~~~~~~~~~~~~~~~

if(zip.mod){

  matern.pri.feat.count <- c(500, .95, .1, .05) ## a, b, c, d
  matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.count[1:2],
                                           prior.sigma = matern.pri.feat.count[3:4])

  comps <- ~
    feat.count.int(1) +
    feat.count.field(cbind(x, y), model = matern.feat.count)

  zip.lik <- bru_obs(
    family = "zeroinflatedpoisson1",
    data = run.dat,
    formula = feat.count ~ feat.count.int + feat.count.field,
    E = total.count
  )

  fit.zip <- try(
  {
    cat('\n')
    for(i in 1){print(glue('{lr.n}: fitting zip model'))}
    cat('\n')

    bru(comps,
        zip.lik,
        # options
        options = list(bru_verbose = 4,
                       bru_max_iter = 1)
        )

  }, silent = T)

  if(exists("fit.zip")){
    try(
    {
      cat('\n')
      for(i in 1){print(glue('{lr.n}:  predicting from fitted zip model'))}
      cat('\n')

      pred.zip <- predict(
        fit.zip, run.dat,
        ~ {
          scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
          lambda <- exp( feat.count.int + feat.count.field )
          expect_param <- lambda * total.count
          expect <- scaling_prob * expect_param
          variance <- scaling_prob * expect_param *
            (1 + (1 - scaling_prob) * expect_param)
          list(
            lambda = lambda,
            expect = expect,
            variance = variance,
            obs_prob = (1 - scaling_prob) * (feat.count == 0) +
              scaling_prob * dpois(feat.count, expect_param)
          )
        },
        n.samples = n.samp
      )
    }, silent = T)
  }

  if(exists("fit.zip")){
    zip_pit <- fit.zip$cpo$pit * c(NA_real_, 1)[1 + (run.dat$feat.count > 0)]
  }else{
    zip_pit <- NA
  }

  if(exists("pred.zip")){
    expect_zip <- pred.zip$expect
    expect_zip$pred_var <- pred.zip$variance$mean + expect_zip$sd^2
    expect_zip$log_score <- -log(pred.zip$obs_prob$mean)
  }else{
    expect_zip$expect <- run.dat
    expect_zip$mean <- NA
    expect_zip$median <- NA
    expect_zip$pred_var <- NA
    expect_zip$log_score <- NA
  }
}

##~~~~~~~~~~~~~~~~~
## zap model
##~~~~~~~~~~~~~~~~~
if(zap.mod){

  matern.pri.feat.count <- c(500, .95, .1, .05) ## a, b, c, d
  matern.pri.feat.present <- c(1000, .95, 10, .05) ## a, b, c, d
  matern.feat.count <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                           constr = TRUE, # integrate-to-zero constraint
                                           prior.range = matern.pri.feat.count[1:2],
                                           prior.sigma = matern.pri.feat.count[3:4])
  ## TODO, tune this prior? it's hard with the 0s and 1s areas in space
  ## and the logit transform
  matern.feat.present <- inla.spde2.pcmatern(mesh=mesh, alpha = 2,
                                             constr = TRUE, # integrate-to-zero constraint
                                             prior.range = matern.pri.feat.present[1:2],
                                             prior.sigma = matern.pri.feat.present[3:4])


  comps <- ~
    ## for modeling the total
    ##  total.count.int(1) +
    ##  total.count.field(cbind(x, y), model = matern.feat.present) +
    # for modeling feature presence/absence
    feat.present.int(1) +
    feat.present.field(cbind(x, y), model = matern.feat.present) +
    # for modeling feature counts (when present is using a zip/zap)
    feat.count.int(1) +
    feat.count.field(cbind(x, y), model = matern.feat.count)

  zap.feat.lik <- bru_obs(
    family = "nzpoisson",
    data = run.dat[feat.count > 0, ],
    # feat.count.lambda  = exp(feat.count.eta)*total.count.lambda
    # feat.count.lambda  = exp(feat.count.eta)*exp(total.count.eta)
    # feat.count.lambda = exp(feat.count.eta + total.count.eta)
    formula = feat.count ~ feat.count.int + feat.count.field, #+ total.count.int + total.count.field,
    E = total.count
  )

  zap.pres.lik <- bru_obs(
    family = "binomial",
    data = run.dat,
    formula = feat.present ~ feat.present.int + feat.present.field
  )

  try(
  {
    cat('\n')
    for(i in 1){print(glue('{lr.n}: fitting zap model'))}
    cat('\n')

    fit.zap <- bru(comps,
                   zap.feat.lik,
                   zap.pres.lik,
                   # options
                   options = list(bru_verbose = 4,
                                  bru_max_iter = 1)
                   )
  }, silent = T)

  if(exists("fit.zap")){
    try(
    {
      cat('\n')
      for(i in 1){print(glue('{lr.n}: predicting from fitted zap model'))}
      cat('\n')

      pred.zap <- predict(
        fit.zap,
        run.dat,
        ~ {
          presence_prob <- plogis( feat.present.int + feat.present.field )
          lambda <- exp(feat.count.int + feat.count.field)
          expect_param <- presence_prob * lambda * total.count
          expect <- expect_param / (1 - exp(-lambda * total.count))
          variance <- expect * (1 - exp(-lambda * total.count) * expect)
          list(
            presence = presence_prob,
            lambda = lambda,
            expect = expect,
            variance = variance,
            obs_prob = (1 - presence_prob) * (feat.count == 0) +
              (feat.count > 0) * presence_prob * dpois(feat.count, expect_param) /
              (1 - dpois(0, expect_param))
          )
        },
        n.samples = n.samp
      )
    }, silent = T)
  }

  if(exists("fit.zap")){
    zap_pit <- rep(NA_real_, nrow(run.dat))
    zap_pit[run.dat$feat.count > 0] <- fit.zap$cpo$pit[-seq_len(nrow(run.dat))]
  }

  if(exists("pred.zap")){
    presence_zap <- pred.zap$presence
    expect_zap <- pred.zap$expect
    expect_zap$pred_var <- pred.zap$variance$mean + expect_zap$sd^2
    expect_zap$log_score <- -log(pred.zap$obs_prob$mean)
  }else{
    expect_zap$expect <- run.dat
    expect_zap$mean <- NA
    expect_zap$median <- NA
    expect_zap$pred_var <- NA
    expect_zap$log_score <- NA
  }
}

######################
## model comparison
######################
cat('\n')
for(i in 1){print(glue('{lr.n}:  model comparison'))}
cat('\n')

comp.df <- data.frame(
  count = rep(run.dat$feat.count, times = 3),
  pred_mean = c(
    expect_poi$mean,
    expect_zip$mean,
    expect_zap$mean
  ),
  pred_var = c(
    expect_poi$pred_var,
    expect_zip$pred_var,
    expect_zap$pred_var
  ),
  pred_median = c(
    expect_poi$median,
    expect_zip$median,
    expect_zap$median
  ),
  log_score = c(
    expect_poi$log_score,
    expect_zip$log_score,
    expect_zap$log_score
  ),
  pit = c(
    poi_pit,
    zip_pit,
    zap_pit
  ),
  Model = rep(c("Poisson", "ZIP", "ZAP"), each = nrow(run.dat))
)

resid.p <- ggplot(comp.df) +
  geom_point(aes(pred_mean, count - pred_mean, color = Model)) +
  ggtitle("Residuals")

pit.p <- ggplot(comp.df) +
  stat_ecdf(aes(pit, color = Model), na.rm = TRUE) +
  scale_x_continuous(expand = c(0, 0)) +
      ggtitle("PIT")

diag.p <- patchwork::wrap_plots(resid.p, pit.p, nrow = 1, guides = "collect")

comp.df <- comp.df %>%
  mutate(
    AE = abs(count - pred_median),
    SE = (count - pred_mean)^2,
    DS = (count - pred_mean)^2 / pred_var + log(pred_var),
    LG = log_score
  )

scores <- comp.df %>%
  group_by(Model) %>%
  summarise(
    MAE = mean(AE, na.rm = T),
    RMSE = sqrt(mean(SE, na.rm = T)),
    MDS = mean(DS, na.rm = T),
    MLG = mean(LG, na.rm = T)
  ) %>%
  left_join(
    data.frame(
      Model = c("Poisson", "ZIP", "ZAP"),
      Order = 1:3
    ),
    by = "Model"
  ) %>%
  arrange(Order) %>%
  select(-Order)
knitr::kable(scores)


#######################
## make plots
#######################

## TODO start here, make comparison plots
# poi (-, density, count, count resid)
# zip (-, density, count, count resid)
# zap (present, density, count, count resid)
# obs (present, density, count)

## TODO then, update save outputs section, and run

cat('\n')
for(i in 1){print(glue('{lr.n}:  making plots'))}
cat('\n')

## first make plots with constant color range per column

prob.zlim <- c(0, 1)
prob.cols <- (magma(256))

dens.zlim <- c(ifelse(exists("pred.poi"), NA, pred.poi$lambda$mean),
               ifelse(exists("pred.zip"), NA, pred.zip$lambda$mean),
               ifelse(exists("pred.zap"), NA, pred.zap$lambda$mean),
               run.dat[, feat.count / total.count]) |> range(na.rm = T)
dens.cols <- cividis(256)

cnt.zlim <- c(ifelse(exists("pred.poi"), NA, pred.poi$expect$mean),
              ifelse(exists("pred.zip"), NA, pred.zip$expect$mean),
              ifelse(exists("pred.zap"), NA, pred.zap$expect$mean),
              run.dat[, feat.count]) |> range(na.rm = T)
cnt.cols <- viridis(256)

res.zlim <- c(ifelse(exists("pred.poi"), NA, run.dat$feat.count - pred.poi$expect$mean),
              ifelse(exists("pred.zip"), NA, run.dat$feat.count - pred.zip$expect$mean),
              ifelse(exists("pred.zap"), NA, run.dat$feat.count - pred.zap$expect$mean),
              run.dat[, feat.count / total.count]) |> range(na.rm = T)
res.cols <- turbo(256)

png(file.path(o.d, glue('{lr.n}-model-prediction-comparisons-fixed-colors.png')),
    width = (qp.res.x / qp.res.y) * 4 / 3 * 13 + 8, height = 4 / 3 * 13, units = 'in', res = 300)
par(mfrow = c(4, 4),
    mai = c(.62, 0.82, .62, 1.22))

if(exists("pred.poi")){

  fields.style();quilt.plot(pred.poi$expect$x,
                            pred.poi$expect$y,
                            1 - dpois(0, pred.poi$expect$mean),
                            main = glue('{lr.n} - Poisson - Estimated Prob of Occurence'),
                            zlim = prob.zlim, nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.poi$lambda$x,
                            pred.poi$lambda$y,
                            pred.poi$lambda$mean,
                            main = glue('{lr.n} - Poisson - Estimated Density'),
                            zlim = dens.zlim, nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.poi$expect$x,
                            pred.poi$expect$y,
                            pred.poi$expect$mean,
                            main = glue('{lr.n} - Poisson - Estimated Count'),
                            zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.poi$expect$x,
                            pred.poi$expect$y,
                            run.dat$feat.count - pred.poi$expect$mean,
                            main = glue('{lr.n} - Poisson - Count Residuals'),
                            zlim = res.zlim, nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
}else{
  for(i in 1:4){ plot.new() }
}

if(exists("pred.zip")){
  fields.style();quilt.plot(pred.zip$expect$x,
                            pred.zip$expect$y,
                            1 - dpois(0, pred.zip$expect$mean),
                            main = glue('{lr.n} - ZIP - Estimated Prob of Occurence'),
                            zlim = prob.zlim, nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zip$lambda$x,
                            pred.zip$lambda$y,
                            pred.zip$lambda$mean,
                            main = glue('{lr.n} - ZIP - Estimated Density'),
                            zlim = dens.zlim, nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zip$expect$x,
                            pred.zip$expect$y,
                            pred.zip$expect$mean,
                            main = glue('{lr.n} - ZIP - Estimated Count'),
                            zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zip$expect$x,
                            pred.zip$expect$y,
                            run.dat$feat.count - pred.zip$expect$mean,
                            main = glue('{lr.n} - ZIP - Count Residuals'),
                            zlim = res.zlim, nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
}else{
  for(i in 1:4){ plot.new() }
}

if(exists("pred.zap")){
  fields.style();quilt.plot(pred.zap$expect$x,
                            pred.zap$expect$y,
                            1 - dpois(0, pred.zap$expect$mean),
                            main = glue('{lr.n} - ZAP - Estimated Prob of Occurence'),
                            zlim = prob.zlim, nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zap$lambda$x,
                            pred.zap$lambda$y,
                            pred.zap$lambda$mean,
                            main = glue('{lr.n} - ZAP - Estimated Density'),
                            zlim = dens.zlim, nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zap$expect$x,
                            pred.zap$expect$y,
                            pred.zap$expect$mean,
                            main = glue('{lr.n} - ZAP - Estimated Count'),
                            zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zap$expect$x,
                            pred.zap$expect$y,
                            run.dat$feat.count - pred.zap$expect$mean,
                            main = glue('{lr.n} - ZAP - Count Residuals'),
                            zlim = res.zlim, nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
}else{
  for(i in 1:4){ plot.new() }
}

fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat$feat.present,
                          main = glue('{lr.n} - Data - Prob of Occurence'),
                          zlim = prob.zlim, nlevel = 256, col = prob.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat[, feat.count / total.count],
                          main = glue('{lr.n} - Data -  Density'),
                          zlim = dens.zlim, nlevel = 256, col = dens.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat[, feat.count],
                          main = glue('{lr.n} - Data -  Counts'),
                          zlim = cnt.zlim, nlevel = 256, col = cnt.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)

fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat[, total.count],
                          main = 'Data -  NUMI',
                          nlevel = 256, col = res.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
dev.off()


png(file.path(o.d, glue('{lr.n}-model-prediction-comparisons-free-colors.png')),
    width = (qp.res.x / qp.res.y) * 4 / 3 * 13 + 8, height = 4 / 3 * 13, units = 'in', res = 300)
par(mfrow = c(4, 4),
    mai = c(.62, 0.82, .62, 1.22))

if(exists("pred.poi")){
  fields.style();quilt.plot(pred.poi$expect$x,
                            pred.poi$expect$y,
                            1 - dpois(0, pred.poi$expect$mean),
                            main = glue('{lr.n} - Poisson - Estimated Prob of Occurence'),
                            nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.poi$lambda$x,
                            pred.poi$lambda$y,
                            pred.poi$lambda$mean,
                            main = glue('{lr.n} - Poisson - Estimated Density'),
                            nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.poi$expect$x,
                            pred.poi$expect$y,
                            pred.poi$expect$mean,
                            main = glue('{lr.n} - Poisson - Estimated Count'),
                            nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.poi$expect$x,
                            pred.poi$expect$y,
                            run.dat$feat.count - pred.poi$expect$mean,
                            main = glue('{lr.n} - Poisson - Count Residuals'),
                            nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
}else{
  for(i in 1:4){ plot.new() }
}

if(exists("pred.zip")){
  fields.style();quilt.plot(pred.zip$expect$x,
                            pred.zip$expect$y,
                            1 - dpois(0, pred.zip$expect$mean),
                            main = glue('{lr.n} - ZIP - Estimated Prob of Occurence'),
                            nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zip$lambda$x,
                            pred.zip$lambda$y,
                            pred.zip$lambda$mean,
                            main = glue('{lr.n} - ZIP - Estimated Density'),
                            nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zip$expect$x,
                            pred.zip$expect$y,
                            pred.zip$expect$mean,
                            main = glue('{lr.n} - ZIP - Estimated Count'),
                            nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zip$expect$x,
                            pred.zip$expect$y,
                            run.dat$feat.count - pred.zip$expect$mean,
                            main = glue('{lr.n} - ZIP - Count Residuals'),
                            nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
}else{
  for(i in 1:4){ plot.new() }
}

if(exists("pred.zap")){
  fields.style();quilt.plot(pred.zap$expect$x,
                            pred.zap$expect$y,
                            1 - dpois(0, pred.zap$expect$mean),
                            main = glue('{lr.n} - ZAP - Estimated Prob of Occurence'),
                            nlevel = 256, col = prob.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zap$lambda$x,
                            pred.zap$lambda$y,
                            pred.zap$lambda$mean,
                            main = glue('{lr.n} - ZAP - Estimated Density'),
                            nlevel = 256, col = dens.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zap$expect$x,
                            pred.zap$expect$y,
                            pred.zap$expect$mean,
                            main = glue('{lr.n} - ZAP - Estimated Count'),
                            nlevel = 256, col = cnt.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  fields.style();quilt.plot(pred.zap$expect$x,
                            pred.zap$expect$y,
                            run.dat$feat.count - pred.zap$expect$mean,
                            main = glue('{lr.n} - ZAP - Count Residuals'),
                            nlevel = 256, col = res.cols,
                            nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
}else{
  for(i in 1:4){ plot.new() }
}

fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat$feat.present,
                          main = glue('{lr.n} - Data - Prob of Occurence'),
                          nlevel = 256, col = prob.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat[, feat.count / total.count],
                          main = glue('{lr.n} - Data -  Density'),
                          nlevel = 256, col = dens.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat[, feat.count],
                          main = glue('{lr.n} - Data -  Counts'),
                          nlevel = 256, col = cnt.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
fields.style();quilt.plot(run.dat$x,
                          run.dat$y,
                          run.dat[, total.count],
                          main = 'Data -  NUMI',
                          nlevel = 256, col = res.cols,
                          nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
dev.off()



## poi.p.d <-  ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = expect_poi, aes(x, y, fill = mean / total.count)) +
##   ggtitle("Poisson predicted density")

## poi.p.c <-  ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = expect_poi, aes(x, y, fill = mean)) +
##   ggtitle("Poisson predicted counts")

## zip.p.d <-  ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = expect_zip, aes(x, y, fill = mean / total.count)) +
##   ggtitle("ZIP predicted density")

## zip.p.c <-  ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = expect_zip, aes(x, y, fill = mean)) +
##   ggtitle("ZIP predicted counts")

## zap.p.p <- ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = presence_zap, aes(x, y, fill = mean)) +
##   geom_point(data = run.dat[feat.count > 0, ], aes(x, y),
##              color = "firebrick", size = 1, pch = 4, alpha = 1) +
##   ggtitle("Presence probability")
## zap.p.d <- ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = expect_zap, aes(x, y, fill = mean / total.count)) +
##   ggtitle("Expected feature counts per total count")
## # geom_sf(data = nests, color = "firebrick", size = 1, pch = 4, alpha = 0.2) +
## zap.p.c <- ggplot() +
##   geom_fm(data = mesh) +
##   geom_raster(data = expect_zap, aes(x, y, fill = mean)) +
##   #    geom_sf(data = nests, color = "firebrick", size = 1, pch = 4, alpha = 0.2) +
##   ggtitle("Expected feature counts per total count")

## patchwork::wrap_plots(zap.p.p, zap.p.d, zap.p.c, nrow = 1)






cat('\n')
for(i in 1){print(glue('{lr.n}:  saving outputs'))}
cat('\n')
## saving each fit objects takes up about 2gb!
## saveRDS(fit.pois.pois, file = file.path(o.d, "fitted-models",
##                                         glue('fit-pois-pois-{lr.n}.rds')))

# save model diagnostic plots
png(file.path(o.d, glue("{lr.n}-diagnostic-plots.png")), width = 10, height = 4, units = 'in', res = 100)
print(diag.p)
dev.off()

# save model comparison metrics
fwrite(scores, file = file.path(o.d, glue("{lr.n}-scores.csv")))

if(exists("pred.poi")){
  saveRDS(pred.poi, file = file.path(o.d, "prediction-objects",
                                     glue('{lr.n}-pred-poi.rds')))
}

if(exists("pred.zip")){
  saveRDS(pred.zip, file = file.path(o.d, "prediction-objects",
                                     glue('{lr.n}-pred-zip.rds')))
}

if(exists("pred.zap")){
  saveRDS(pred.zap, file = file.path(o.d, "prediction-objects",
                                     glue('{lr.n}-pred-zap.rds')))
}

# clean up
rm(run.dat, fit.poi, fit.zip, fit.zap,
   pred.poi, pred.zip, pred.zap,
   expect_poi, expect_zip, expect_zap);gc()


## TODO: CRPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ## get the CRPS as in:
## ## https://inlabru-org.github.io/inlabru/articles/prediction_scores.html
## # some large value, so that 1-F(K) is small
## pred <- generate(fit.pois.pois, run.dat,

##                  formula = ~ {

##                    # linear predictors for lig, feat, and remainder
##                    eta.l.feat <- feat.l.count.int + feat.l.count.field
##                    eta.r.feat <- feat.r.count.int + feat.r.count.field
##                    eta.remain <- remain.count.int + remain.count.field

##                    # rates for lig, feat, and remainder
##                    lambda.feat.l <- exp(feat.l.count.int + feat.l.count.field)
##                    lambda.feat.r <- exp(feat.r.count.int + feat.r.count.field)
##                    lambda.remain <- exp(remain.count.int + remain.count.field)
##                    # they are additive for the total on the count scale
##                    lambda.total <- lambda.feat.l + lambda.feat.r + lambda.remain

##                    # grab the observations
##                    feat.l.count <- run.dat[, feat.l.count]
##                    feat.r.count <- run.dat[, feat.r.count]
##                    total.count  <- run.dat[, total.count]

##                    # max Ks to evaluate
##                    max.K.l <- max(feat.l.count) + 4 * sqrt(max(feat.l.count))
##                    max.K.r <- max(feat.r.count) + 4 * sqrt(max(feat.r.count))
##                    max.K.t <- max(total.count)  + 4 * sqrt(max(total.count))
##                    k.l <- seq(0, max.K.l)
##                    k.r <- seq(0, max.K.r)
##                    k.t <- seq(0, max.K.t)

##                    crps.l.df <- do.call(
##                      cbind,
##                      lapply(
##                        seq_along(feat.l.count),
##                        function(i) {
##                          Fpred <- ppois(k.l, lambda = lambda.feat.l[i])
##                          data.frame(
##                            k = c(k.l, k.l),
##                            i = c(i, i),
##                            type = rep(c("F", "residual"), each = length(Fpred)),
##                            value = c(Fpred, Fpred - (feat.l.count[i] <= k.l))
##                          )
##                        }
##                      )
##                    )

##                    crps.r.df <- do.call(
##                      cbind,
##                      lapply(
##                        seq_along(feat.r.count),
##                        function(i) {
##                          Fpred <- ppois(k.r, lambda = lambda.feat.r[i])
##                          data.frame(
##                            k = c(k.r, k.r),
##                            i = c(i, i),
##                            type = rep(c("F", "residual"), each = length(Fpred)),
##                            value = c(Fpred, Fpred - (feat.r.count[i] <= k.r))
##                          )
##                        }
##                      )
##                    )

##                    crps.t.df <- do.call(
##                      rbind, # this was cbind in the inlabru code
##                      lapply(
##                        seq_along(total.count),
##                        function(i) {
##                          Fpred <- ppois(k.t, lambda = lambda.total[i])
##                          data.frame(
##                            k = c(k.t, k.t),
##                            i = c(i, i),
##                            type = rep(c("F", "residual"), each = length(Fpred)),
##                            value = c(Fpred, Fpred - (total.count[i] <= k.t))
##                          )
##                        }
##                      )
##                    )

##                    # return list of everything we want
##                    list(lambda.total = lambda.total,
##                         #total.obs.prob = NULL, # TODO
##                         lambda.remain = lambda.remain,
##                         lambda.feat.l = lambda.feat.l,
##                         lambda.feat.r = lambda.feat.r,
##                         lambda.feat.l.density.per.count = lambda.feat.l / lambda.total,
##                         lambda.feat.r.density.per.count = lambda.feat.r / lambda.total,
##                         crps.l = crps.l.df,
##                         crps.r = crps.r.df,
##                         crps.t = crps.t.df
##                         #feat.var = NULL, # TODO
##                         #feat.obs.prob = NULL # TODO
##                         )
##                  },
##                  n.samples = 5
##                  )

## ## restructure pred pieces by type, instead of pieces by draw
## pred.l <- list()
## pred.l[['lambda.total']] <- do.call('cbind',
##                                     lapply(pred,
##                                            function(x){
##                                              x$lambda.total}))
## pred.l[['lambda.remain']] <- do.call('cbind',
##                                      lapply(pred,
##                                             function(x){
##                                               x$lambda.remain}))
## pred.l[['lambda.feat.l']] <- do.call('cbind',
##                                      lapply(pred,
##                                             function(x){
##                                               x$lambda.feat.l}))
## pred.l[['lambda.feat.r']] <- do.call('cbind',
##                                      lapply(pred,
##                                             function(x){
##                                               x$lambda.feat.r}))
## pred.l[['lambda.feat.l.density.per.count']] <- do.call('cbind',
##                                                        lapply(pred,
##                                                               function(x){
##                                                                 x$lambda.feat.l.density.per.count}))
## pred.l[['lambda.feat.r.density.per.count']] <- do.call('cbind',
##                                                        lapply(pred,
##                                                               function(x){
##                                                                 x$lambda.feat.r.density.per.count}))
## pred.l[['crps.t']] <- lapply()



## F.est.l <-
##   (pred.dt %>%
##      filter(type == "F") %>%
##      group_by(i) %>%
##      summarise(F = sum(mean), groups = "drop") %>%
  ##      pull("F"))
  ## crps.score.l <-
  ##   (pred$crps.l %>%
  ##      filter(type == "residual") %>%
  ##      group_by(i) %>%
  ##      summarise(crps = sum(mean^2), groups = "drop") %>%
  ##      pull(crps))
  ## # Check that the cutoff point K has nearly probability mass 1 below it,
  ## # for all i:
  ## min(F_estimate)


  ## END OF CRPS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## OLDER CODE BELOW - aoz april 2025

  #pdf(file.path(o.d, 'data-vs-model-pois-pois.pdf'), width = 13, height = 7)

  ## png(file.path(o.d, glue('data-vs-model-pois-pois-{lr.n}.png')),
  ##     width = (qp.res.x / qp.res.y) * 13 + 8, height = 13, units = 'in', res = 300)
  ## par(mfrow = c(3, 3),
  ##     mai = c(.62, 0.82, .62, 1.22))
  ## fields.style();quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
  ##                           main = glue('Observed {lr.n} Counts'),
  ##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x, cex = 1.4)
  ## fields.style();quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count / total.count],
  ##                           main = glue('Observed {lr.n} Counts per Total Counts'),
  ##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  ## fields.style();quilt.plot(run.dat[, x], run.dat[, y], run.dat[, total.count],
  ##                           main = 'Total Counts',
  ##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  ## fields.style();quilt.plot(pred.pois.pois$feat[, x],
  ##                           pred.pois.pois$feat[, y],
  ##                           pred.pois.pois$feat[, mean],
  ##                           main = glue('Estimated {lr.n} Counts' ),
  ##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
  ## fields.style();quilt.plot(pred.pois.pois$feat[, x],
##                           pred.pois.pois$feat[, y],
##                           pred.pois.pois$feat.density.per.count[, mean],
##                           main = glue('Estimated {lr.n} Counts per Total Counts' ),
##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
## fields.style();quilt.plot(pred.pois.pois$total[, x],
##                           pred.pois.pois$total[, y],
##                           pred.pois.pois$total[, mean],
##                           main = glue('Estimated Total Counts'),
##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
## if(!('feat.count' %in% colnames(pred.pois.pois$feat))){
##   pred.pois.pois$feat <-
##     merge(pred.pois.pois$feat,
##           run.dat[, .(x, y, feat.count, total.count)])
## }
## fields.style();quilt.plot(pred.pois.pois$feat[, x],
##                           pred.pois.pois$feat[, y],
##                           pred.pois.pois$feat[, feat.count] - pred.pois.pois$feat[, mean],
##                           main = glue('Residual {lr.n} Counts' ),
##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
## fields.style();quilt.plot(pred.pois.pois$feat[, x],
##                           pred.pois.pois$feat[, y],
##                           pred.pois.pois$feat[, feat.count / total.count] -
##                             pred.pois.pois$feat.density.per.count[, mean],
##                           main = glue('Residual {lr.n} Counts per Total Counts' ),
##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
## fields.style();quilt.plot(pred.pois.pois$total[, x],
##                           pred.pois.pois$total[, y],
##                           pred.pois.pois$feat[, total.count] - pred.pois.pois$total[, mean],
##                           main = glue('Residual Total Counts' ),
##                           nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x)
## dev.off()

## # clean up
## rm(run.dat, fit.pois.pois, pred.pois.pois);gc()


## ########################
## ## feature only model ##
## ########################
## fit.feat.only.zip <- bru(
##   # all model components
##   zip.comps,
##   # feature only lik
##   feat.only.count.zip.l,
##   # options
##   options = list(bru_verbose = 4,
##                  bru_max_iter = 1
##                  )
## )
## pred.feat.only.zip <- predict(
##   fit.feat.only.zip, run.dat,
##   ~ {
##     scaling_prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
##     lambda <- exp( feat.count.int + feat.count.field )
##     expect_param <- lambda * 1
##     expect <- scaling_prob * expect_param
##     variance <- scaling_prob * expect_param *
##       (1 + (1 - scaling_prob) * expect_param)
##     list(
##       feat = expect,
##       feat.var = variance,
##       feat.obs.prob = (1 - scaling_prob) * (feat.present == 0) +
##         scaling_prob * dpois(feat.count, expect_param)
##     )
##   },
##   n.samples = 100
## )
## #pdf(file.path(o.d, 'data-vs-model-feat-only-zip.pdf'), width = 13, height = 7)
## png(file.path(o.d, 'data-vs-model-feat-only-zip.png'), width = 13, height = 7, units = 'in', res = 300)
## par(mfrow = c(3, 3))
## quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count],
##            main = glue('Observed {lr.n} Counts' ))
## quilt.plot(run.dat[, x], run.dat[, y], run.dat[, feat.count / total.count],
##            main = glue('Observed {lr.n} Counts per Total Counts'))
## quilt.plot(run.dat[, x], run.dat[, y], run.dat[, total.count],
##            main = 'Total Counts')
## quilt.plot(pred.feat.only.zip$feat[, x],
##            pred.feat.only.zip$feat[, y],
##            pred.feat.only.zip$feat[, mean],
##            main = glue('Estimated {lr.n} Counts' ))
## quilt.plot(pred.feat.only.zip$feat[, x],
##            pred.feat.only.zip$feat[, y],
##            pred.feat.only.zip$feat.var[, mean],
##            main = glue('Variance of {lr.n} Counts per Total Counts' ))
## dev.off()
## ####################
## ## pois zip model ##
## ####################
## fit.pois.zip <- bru(
##   # all model components
##   zip.comps,
##   # total counts lik
##   total.count.pois.l,
##   # feature lik
##   feat.count.zip.l,
##   # options
##   options = list(bru_verbose = 4,
##                  bru_max_iter = 1)
## )
## fit.pois.zap <- bru(
##   zip.comps,
##   total.count.pois.l, feat.count.zap.l, feat.present.zap.l,
##   options = list(bru_verbose = 4)
## )

## ## get list of bru standardised names
## pred.pois.pois.1 <- generate(fit.pois.pois, n.samples = 1)[[1]]
## pred.zip.1 <- generate(fit.zip, n.samples = 1)[[1]]
## #############
## ## predict ##
## #############



## pred.zip <- predict(
##   fit.zip, run.dat,
##   ~ {
##     total.scaling.prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
##     total.lambda <- exp(total.int + total.field)
##     total.expect.param <- total.lambda * 1
##     total.expect <- total.scaling.prob * total.expect.param
##     total.var <- total.scaling.prob * total.expect.param *
##       (1 + (1 - total.scaling.prob) * total.expect.param)

##     feat.scaling.prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1_2_)
##     feat.lambda <- exp(feat.int + feat.field + total.int + total.field)
##     feat.expect.param <- feat.lambda * 1
##     feat.expect <- feat.scaling.prob * feat.expect.param
##     feat.var <- feat.scaling.prob * feat.expect.param *
##       (1 + (1 - feat.scaling.prob) * feat.expect.param)

##     list(total.lambda = total.lambda,
##          total.expect = total.expect,
##          total.var = total.var,
##          total.obs.prob = (1 - total.scaling.prob) * (total == 0) +
##            total.scaling.prob * dpois(total, total.expect.param),
##          total.scaling.prob = total.scaling.prob,
##          feat.lambda = feat.lambda,
##          feat.expect = feat.expect,
##          feat.var = feat.var,
##          feat.obs.prob = (1 - feat.scaling.prob) * (feat.count == 0) +
##            feat.scaling.prob * dpois(feat.count, feat.expect.param),
##          feat.scaling.prob = feat.scaling.prob
##          )
##   },
##   n.samples = 50
## )


## total.expect.zip <- pred.zip$total.expect
## total.expect.zip$pred_var <- pred.zip$total.var$mean + total.expect.zip$sd^2
## total.expect.zip$log.score <- -log(pred.zip$total.obs.prob$mean)
## feat.expect.zip <- pred.zip$feat.expect
## feat.expect.zip$pred_var <- pred.zip$feat.var$mean + feat.expect.zip$sd^2
## feat.expect.zip$log.score <- -log(pred.zip$feat.obs.prob$mean)

## total.zip.p <- ggplot(data = total.expect.zip,
##                       aes(x, y, fill = mean)) +
##   geom_tile()
## feat.zip.p <- ggplot(data = feat.expect.zip,
##                      aes(x, y, fill = mean)) +
##   geom_tile()

## +
##   geom_sf(data = nests, color = "firebrick", size = 1, pch = 4, alpha = 0.2) +
##   ggtitle("Nest intensity per ~1.25 ha")

## stk.y <- inla.stack(data = list(amount = y, # for single model
##                                 alldata = cbind(y, NA)), # for joint model
##                     A = list(A.est, 1),
##                     effects = list(
##                       list(y.field = 1:spde$n.spde),
##                       list(y.intercept = rep(1, length(y)))),
##                     tag = "est.y"
##                     )


## stk.z <- inla.stack(data = list(occurrence = z, # for single model
##                                 alldata = cbind(NA, z)), # for joint model
##                     A = list(A.est, 1),
##                     effects = list(
##                       list(z.field = 1:spde$n.spde,
##                            zc.field = 1:spde$n.spde),
##                       list(z.intercept = rep(1, length(z)))),
##                     tag = "est.z"
##                     )

## ## make a quiet version of INLA to parse logs
## #inla.q <- quietly(inla)

## # spatial response only



## ## # spatial binomial only
## ## f.z <- occurrence ~ -1 + z.intercept + f(z.field, model = spde)
## ## out.z <- inla(f.z, family = "binomial",
## ##               data = inla.stack.data(stk.z),
## ##               control.predictor = list(A = inla.stack.A(stk.z), compute = TRUE),
## ##               control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
## ##               control.inla = list(strategy = "laplace"))

## ## print(glue("time to fit spatial binom only is: {toc()[3]}"))

## #######################
## # spatial joint model #
## # the real deal       #
## #######################

## full.start.time <- toc()[3] / 60
## stk.yz <- inla.stack(stk.y, stk.z)
## f.yz <- alldata ~ -1 + y.intercept + z.intercept +
##   f(y.field, model = spde) + f(z.field, model = spde) +
##   f(zc.field, copy = "y.field", fixed = FALSE)
## out.yz <- inla(f.yz, family = c("gamma", "binomial"),
##                data = inla.stack.data(stk.yz),
##                control.predictor = list(A = inla.stack.A(stk.yz), compute = TRUE),
##                control.compute = list(dic = TRUE, config = TRUE, return.marginals.predictor=TRUE),
##                control.inla = list(strategy = "laplace"),
##                verbose = F)

## full.start.elapsed.time <- (glue("time to fit full model at time {time.i} is: {round((toc()[3]/60 - full.start.time) ,2)} min"))
## # time.mat[2, time.i] <- full.start.elapsed.time # for profile timing
## print(full.start.elapsed.time)

## ## # non-spatial joint
## ## f.nospatial <- alldata ~ -1 + z.intercept + y.intercept
## ## out.nospatial <- inla(f.nospatial, family=c("gamma", "binomial"),
## ##                       data=inla.stack.data(stk.yz),
## ##                       control.predictor=list(A=inla.stack.A(stk.yz)),
## ##                       control.compute=list(dic=TRUE),
## ##                       control.inla=list(strategy="laplace"),
## ##                       verbose = F)

## ## # joint but not shared spatial component
## ## f.noshared <- alldata ~ -1 + y.intercept + z.intercept +
## ##   f(y.field, model=spde) + f(z.field, model=spde)
## ## out.noshared <- inla(f.noshared, family=c("gamma", "binomial"),
## ##                      data=inla.stack.data(stk.yz),
## ##                      control.predictor=list(A=inla.stack.A(stk.yz)),
## ##                      control.compute=list(dic=TRUE),
## ##                      control.inla=list(strategy="laplace"))

## ##save(file = file.path(out.top.dir, out.sub.dir, "dual-likelihood-test-fit.Rdata"), list = ls())
## ##load(file.path(out.top.dir, "2024-02-15/dual-likelihood-test-fit.Rdata"))

## ## perform model comparison between the three models using DIC
## # note: need to be careful with the DIC when using multiple likelihoods (joint)
## #       must sum local DIC for each obs
## # we want lower DIC which indicates better fit

## ## TODO rewrite all numerical outputs into one output row -> csv
## ## write.csv(rbind(nospatial=tapply(out.nospatial$dic$local.dic, out.nospatial$dic$family, sum),
## ##                 noshared=tapply(out.noshared$dic$local.dic, out.noshared$dic$family, sum),
## ##                 separate=c(out.y$dic$dic, out.z$dic$dic),
## ##                 joint=tapply(out.yz$dic$local.dic, out.yz$dic$family, sum)),
## ##           file.path(o.d, glue("{lr.n}-DIC-comp.csv")))

## ## calculate the expected value of the posterior marginal
## ## distributions on the response scale
## idy <- which(z == 1)
## exp.y <- sapply(out.yz$marginals.linear.predictor[idy],
##                 function(m){inla.emarginal(exp, m)})

## ## TODO is this which call right? it's what is in the INLA book pg 276
## idz <- which(!is.na(z))
## ## NOTE 28aug2024 some of the marginals for z (pre transform) have Inf height
## ## this breaks inla.emarginal, this the tryCatch. But it's not a big deal since we
## ## exp.z is only used to store mean(exp.z) as a summary metric
## exp.z <- tryCatch({sapply(out.yz$marginals.linear.predictor[length(lr.dat) + idz],
##                           function(m){inla.emarginal(inla.link.invlogit, m)})
## },
## error = function(err){
##   NA
## })


## ## # compare these values to the observed mean for the signal occurrence
## ## # and amount -- a nice simple sanity check to make sure everything is
## ## # as expected
## ## means.comp <- matrix(c(yPositive.obs = mean(lr.dat[which(lr.dat > 0)]),
## ##                        yPositive.pred = mean(exp.y),
## ##                        z.obs = mean(z, na.rm = T), z.pred = mean(exp.z)),
## ##                      nrow = 2, byrow = T)
## ## colnames(means.comp) <- c("obs", "pred")
## ## rownames(means.comp) <- c("yPositive", "z")
## ## write.csv(means.comp, file.path(o.d, glue("{lr.n}-obs-est-means-comp.csv")))

## ## plot the posterior marginal distributions for y.int, z.int,
## ## tau for separate and joint models, and beta_1 for joint
## ## as in figure 8.10 pg 277
## #print("starting first plot")

## png(file.path(o.d, glue("{lr.n}-{time.i}-post-marg-dist.png")),
##     width = 8, height = 8, units = "in", res = 300)
## par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=2:0)
## plot(inla.smarginal(out.yz$marginals.fixed[[1]]), type='l', ylab='Density',
##      ylim=c(0,max(out.yz$marginals.fixed[[1]][,2], out.y$marginals.fixed[[1]][,2])),
##      xlab=expression(b[0]^y))
## lines(inla.smarginal(out.y$marginals.fixed[[1]]), lty=2)
## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## plot(inla.smarginal(out.yz$marginals.fixed[[2]]), type='l', ylab='Density',
##      #ylim=c(0,max(out.yz$marginals.fixed[[2]][,2], out.z$marginals.fixed[[1]][,2])),
##      xlab=expression(b[0]^z))
## #lines(inla.smarginal(out.z$marginals.fixed[[1]]), lty=2)
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
## png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-effects-ordered.png")),
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
  ##                        observed.y.out=sum(out.yz$summary.random$y.field[in.d,4]>0 |
  ##                                             out.yz$summary.random$y.field[in.d,6]<0),
  ##                        observed.z.out=sum(out.yz$summary.random$z.field[in.d,4]>0 |
  ##                                             out.yz$summary.random$z.field[in.d,6]<0))
  ## #write.csv(spat.sig, file.path(o.d, glue("{lr.n}-spat-sig.csv")), row.names = F)

  ## out.yz.summ <- data.frame(lr.n = lr.n,
  ##                           time = time.i,
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
  ##                                   mean(lr.dat[which(lr.dat > 0)]),
  ##                                   mean(exp.y),
  ##                                   mean(z, na.rm = T),
  ##                                   mean(exp.z),
  ##                                   as.numeric(spat.sig[1, ]))
  ##                           )

  ## write.csv(out.yz.summ, file.path(o.d, glue("{lr.n}-{time.i}-output-summary.csv")), row.names = F)


  ## ## extract the spatial effects
  ## y.field <- inla.spde2.result(out.yz, "y.field", spde)
  ## z.field <- inla.spde2.result(out.yz, "z.field", spde)

  ## png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-effects-marginals.png")),
  ##     width = 8, height = 8, units = "in", res = 300)
  ## par(mfrow=c(2,2))
  ## plot.default(y.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
  ## plot.default(y.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
  ## plot.default(z.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
  ## plot.default(z.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')
  ## dev.off()

  ## #### plot spatial random fields

  ## ## make projection grid

  ## # predict resolution (number of pixels in x and y)
  ## xy.pred <- make.pred.grid(xy.obs = xy.obs,
  ##                           retain.shape = T,
  ##                           res = pred.grid.res,
  ##                           boundary.poly = boundary.coords)
  ## if(j.i  == 1){
  ##   write.csv(xy.pred, file.path(o.d, glue("xy-pred-locs-{time.i}.csv")), row.names = F)
  ## }

  ## A.pred <- inla.spde.make.A(
  ##   mesh = mesh,
  ##   loc = xy.pred)

  ## rf.grid <- list(xi.mean = fm_evaluate(mesh, field = out.y$summary.random$y.field$mean, loc = xy.pred),
  ##                 xi.sd = fm_evaluate(mesh, field = out.y$summary.random$y.field$sd, loc = xy.pred),
  ##                 obs = lr.dat,
  ##                 u.mean = fm_evaluate(mesh, field = out.yz$summary.random$z.field$mean, loc = xy.pred),
  ##                 u.sd = fm_evaluate(mesh, field = out.yz$summary.random$z.field$sd, loc = xy.pred),
  ##                 occ.obs = z)

  ## xy.in <- inout(xy.pred, boundary.coords)
  ## for (i in c(1:2, 4:5)) rf.grid[[i]][!xy.in] <- NA


  ## png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-fields.png")),
  ##       width = 12, height = 8, units = "in", res = 500)
  ##   par(mfrow=c(2,3), mar=c(0,0,0,0))
  ##   for (i in 1:6) {
  ##     if(grepl("obs", names(rf.grid)[i])){
  ##       quilt.plot(xy.obs, rf.grid[[i]],
  ##                  nx = 60, ny = 50,
  ##                  nlevel=100, col=viridis(100),
  ##                  legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
  ##                  xaxt = 'n', yaxt = 'n')
  ##     }else{
  ##       quilt.plot(xy.pred, rf.grid[[i]],
  ##                  nx = pred.grid.res * .9, ny = pred.grid.res * .9,
  ##                  nlevel=101, col=viridis(100),
  ##                  legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
  ##                  xaxt = 'n', yaxt = 'n')
  ##     }
  ##     lines(boundary.coords, col = "orange", lwd = 2)
  ##   }
  ##   dev.off()

  ##   #####################################
  ##   #--- Prediction of the responses ---#
  ##   #####################################

  ##   s1 <- inla.posterior.sample(n=1,result=out.yz)

  ##   names(s1[[1]])
  ##   grep("y.intercept", rownames(s1[[1]]$latent), fixed=TRUE)

  ##   ids <- lapply(c('y.intercept', 'z.intercept', 'y.field', 'z.field', 'zc.field'), function(x) grep(x, rownames(s1[[1]]$latent), fixed=TRUE))
  ##   pred.y.f <- function(s) exp(s$latent[ids[[1]], 1] + s$latent[ids[[3]], 1])
  ##   pred.z.f <- function(s) 1/(1 + exp(-(s$latent[ids[[2]], 1] + s$latent[ids[[4]], 1] + s$latent[ids[[5]], 1]))) ## NOTE! the latent zc.field is already beta.zc*y.field

  ##   s1000 <- inla.posterior.sample(1000, out.yz)

  ##   prd.y <- sapply(s1000, pred.y.f)
  ##   prd.z <- sapply(s1000, pred.z.f)
  ##   prd.y.cond.z <- matrix(rbinom(n = length(prd.z), size = 1, p = prd.z),
  ##                          ncol = ncol(prd.z))* prd.y

  ##   # project all draws to pred grid
  ##   prd.y.proj <- apply(prd.y, 2,
  ##                       function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
  ##   prd.z.proj <- apply(prd.z, 2,
  ##                       function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
  ##   prd.y.cond.z.proj <- apply(prd.y.cond.z, 2,
  ##                              function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})


  ##   # always write the means
  ##   write.csv(as.matrix(rowMeans(prd.y.proj), ncol = 1),
  ##             file.path(o.d, glue("{lr.n}-{time.i}-resp-est.csv")), row.names = F)
  ##   write.csv(as.matrix(rowMeans(prd.y.cond.z.proj), ncol = 1),
  ##             file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-est.csv")), row.names = F)
  ##   write.csv(as.matrix(apply(prd.y.cond.z.proj, 1, sd), ncol = 1),
  ##             file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-sd.csv")), row.names = F)

  ##   # only save some of the complete draws -- these are big files (.5Gb)
  ##   # and only save the "final" field draws
  ##   if(lr.n %in% lr.to.save.draws){
  ##     write.csv(prd.y.cond.z.proj,
  ##               file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-draws.csv")), row.names = F)
  ##   }

  ##   # prepare to plot obs, est, sd for y, z, y|z

  ##   prd <- list(resp.obs = lr.dat)
  ##   prd$y.mean <- fm_evaluate(mesh, field=rowMeans(prd.y), loc = xy.pred)
  ##   prd$y.sd <- fm_evaluate(mesh, field=apply(prd.y, 1, sd), loc = xy.pred)

  ##   prd$occ.obs <- z
  ##   prd$z.mean <- fm_evaluate(mesh, field=rowMeans(prd.z), loc = xy.pred)
  ##   prd$z.sd <- fm_evaluate(mesh, field=apply(prd.z, 1, sd), loc = xy.pred)

  ##   prd$cond.resp.obs <- lr.dat
  ##   prd$y.cond.z.mean <- fm_evaluate(mesh, field=rowMeans(prd.y.cond.z), loc = xy.pred)
  ##   prd$y.cond.z.sd <- fm_evaluate(mesh, field=apply(prd.y.cond.z, 1, sd), loc = xy.pred)

  ##   # knock out places where we don't want to plot/predict
  ##   for (j in c(2:3, 5:6, 8:9)) prd[[j]][!xy.in] <- NA

  ##   titles <- c("Y Obs", "Y Est", "Y SD",
  ##               "Z Obs", "Z Est", "Z SD",
  ##               "Y|Z Obs", "Y|Z Est", "Y|Z SD")

  ##   # *** Code for Figure 8.14
  ##   png(file.path(o.d, glue("{lr.n}-{time.i}-response-occurrence-fields.png")),
  ##       width = 12, height = 12, units = "in", res = 500)
  ##   par(mfcol=c(3,3), mar=c(0,0,0,1))
  ##   for (i in 1:length(prd)) {

  ##     # set xlim and ylim for consistency
  ##     xl <- range(xy.obs[, 1]) + c(-.05, .05) * (max(xy.obs[, 1]) - min(xy.obs[, 1]))
  ##     yl <- range(xy.obs[, 2]) + c(-.05, .05) * (max(xy.obs[, 2]) - min(xy.obs[, 2]))

  ##     ## first, set z/color ranges and flag large predicted values in estimates
  ##     if(i %in% c(1, 2, 7, 8)){

  ##       zl <- c(0, max(prd$resp.obs) * 1.01) ## mean and data obs plots
  ##       zcols <- c(viridis(100), "#b33609")

  ##       # this should be a shrinkage/smoothing model. outputs above max observed are truncated and
  ##       # plotted as dark red
  ##       prd[[i]][prd[[i]] > max(prd$resp.obs) * 1.01] <- max(prd$resp.obs) * 1.01

  ##     }else if(i %in% c(4, 5)){## prob of non-zero response mean/est and obs

  ##       zl <- 0:1
  ##       if(i == 4) zcols <- viridis(2)
  ##       if(i == 5) zcols <- viridis(100)

  ##     }else{ ## sd plots, no scaling

  ##       zl <- c(0, max(prd[[i]], na.rm = T))
  ##       zcols <- turbo(100)

  ##     }

  ##     ## second, plot the model outputs differently than the raw data

  ##     if(grepl("obs", names(prd)[i])){
  ##       bubblePlot(xy.obs[, 1], xy.obs[, 2], prd[[i]],
  ##                  col=zcols,
  ##                  legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
  ##                  legend.args=list(text="", line=0),
  ##                  xaxs = 'i', yaxs = 'i', # o.w. R adds 4% to pad axes and
  ##                  # then plot this won't align with quilt.plot
  ##                  xaxt = "n", yaxt = "n", bty = "n",
  ##                  xlim = xl, ylim = yl, zlim = zl)
  ##     }else{
  ##       quilt.plot(xy.pred, prd[[i]],
  ##                  nx = pred.grid.res * .9, ny = pred.grid.res * .9,
  ##                  nlevel=length(zcols), col=zcols,
  ##                  legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
  ##                  legend.args=list(text="", line=0),
  ##                  xaxt = 'n', yaxt = 'n', bty = "n",
  ##                  xlim = xl, ylim = yl, zlim = zl)
  ##     }

  ##     ## add title and boundary lines
  ##     text(max(xy.obs[, 1]), max(xy.obs[, 2]), adj = c(1, 0.5), titles[i], cex = 2)
  ##     lines(boundary.coords, lwd = 2)
  ##   }
  ##   dev.off()


  ##   ## png(file.path(o.d, "test.png"), width = 5, height = 5, units = "in", res = 300)
  ##   ## DefaultAssay(brain)<- "NeighborhoodToCell_ALRA"
  ##   ## SpatialFeaturePlot(brain, features = lr.n, slot = "data") +
  ##   ##   scale_fill_viridis_c(limits = c(0.2, 1))
  ##   ## dev.off()

  ##               final.elapsed.time <- (glue("final elapsed time at time {time.i} : {round(toc()[3]/60,2)} min"))
  ##               # time.mat[3, time.i] <- final.elapsed.time # for profile timing
  ##               print(final.elapsed.time)

  ## ######################################################################
  ## ######################################################################
  ## ######################################################################

  ## ## first attempted dual ZIP model is saved below


  ## ## ############################
  ## ## ## define the objects needed for a zero-inflated poisson model
  ## ## ############################

  ## ## ## define all the components in the model
  ## ## zip.comps <- ~
  ## ##   total.int(1) +
  ## ##   feat.int(1) +
  ## ##   total.field(cbind(x, y), model = matern) +
  ## ##   feat.field(cbind(x, y), model = matern)
  ## ## #  total.present.int(1) + # for ZAP model
  ## ## #  feat.present.int(1) + # for ZAP model
  ## ## #  total.present.field(geometry, model = matern) + # for ZAP model
  ## ## #  feat.present.field(geometry, model = matern) # for ZAP model

  ## ## ## define the likelihoods
  ## ## # we use a zero-inflated poisson model for total counts
  ## ## # which assumes the binomial zero process is correlated with the non-zero process:
  ## ## # P(y|...) = p*I(y=0) + (1-p)*Poisson(y)
  ## ## # (as opposed to a zero-adjusted poisson, which allows for the binomial zero process to significantly differ from the non-zero poisson process:
  ## ## # P(y|...) = p*I(y=0) + (1-p)*Poisson(y|y>0)
## ## # see https://inlabru-org.github.io/inlabru/articles/zip_zap_models.html
## ## # for more on differences and how to code
## ## zip.total.l <- inlabru::bru_obs(
## ##   family = "zeroinflatedpoisson1",
## ##   data = run.dat,
## ##   formula = total.count ~ total.int + total.field,
## ##   E = 1
## ## )
## ## zip.feat.l <- bru_obs(
## ##   family = "zeroinflatedpoisson1",
## ##   data = run.dat,
## ##   formula = feat.count ~ {
## ##     eta.feat.rate.per.count <- feat.int + feat.field
## ##     eta.E.count   <- total.int + total.field
## ##     # feat.rate = feat.rate.per.count*count
## ##     # eta.feat = log(feat.rate) = log(feat.rate.per.count) + log(count)
## ##     eta.feat <- eta.feat.rate.per.count + eta.E.count
## ##     eta.feat
## ##   },
## ##   E = 1
## ## )
## ## ## run inlabru
## ## fit.zip <- bru(
## ##   zip.comps,
## ##   zip.total.l, zip.feat.l,
## ##   options = list(bru_verbose = 4)
## ## )
## ## ## get list of bru standardised names
## ## pred.zip.1 <- generate(fit.zip, n.samples = 1)[[1]]

## ## ## predict

## ## pred.zip <- predict(
## ##   fit.zip, run.dat,
## ##   ~ {
## ##     total.scaling.prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1)
## ##     total.lambda <- exp(total.int + total.field)
## ##     total.expect.param <- total.lambda * 1
## ##     total.expect <- total.scaling.prob * total.expect.param
## ##     total.var <- total.scaling.prob * total.expect.param *
## ##       (1 + (1 - total.scaling.prob) * total.expect.param)

## ##     feat.scaling.prob <- (1 - zero_probability_parameter_for_zero_inflated_poisson_1_2_)
## ##     feat.lambda <- exp(feat.int + feat.field + total.int + total.field)
## ##     feat.expect.param <- feat.lambda * 1
## ##     feat.expect <- feat.scaling.prob * feat.expect.param
## ##     feat.var <- feat.scaling.prob * feat.expect.param *
## ##       (1 + (1 - feat.scaling.prob) * feat.expect.param)

## ##     list(total.lambda = total.lambda,
## ##          total.expect = total.expect,
## ##          total.var = total.var,
## ##          total.obs.prob = (1 - total.scaling.prob) * (total == 0) +
## ##            total.scaling.prob * dpois(total, total.expect.param),
## ##          total.scaling.prob = total.scaling.prob,
## ##          feat.lambda = feat.lambda,
## ##          feat.expect = feat.expect,
## ##          feat.var = feat.var,
## ##          feat.obs.prob = (1 - feat.scaling.prob) * (feat.count == 0) +
## ##            feat.scaling.prob * dpois(feat.count, feat.expect.param),
## ##          feat.scaling.prob = feat.scaling.prob
## ##          )
## ##   },
## ##   n.samples = 50
## ## )


## ## total.expect.zip <- pred.zip$total.expect
## ## total.expect.zip$pred_var <- pred.zip$total.var$mean + total.expect.zip$sd^2
## ## total.expect.zip$log.score <- -log(pred.zip$total.obs.prob$mean)
## ## feat.expect.zip <- pred.zip$feat.expect
## ## feat.expect.zip$pred_var <- pred.zip$feat.var$mean + feat.expect.zip$sd^2
## ## feat.expect.zip$log.score <- -log(pred.zip$feat.obs.prob$mean)

## ## total.zip.p <- ggplot(data = total.expect.zip,
## ##                       aes(x, y, fill = mean)) +
## ##   geom_tile()
## ## feat.zip.p <- ggplot(data = feat.expect.zip,
## ##                      aes(x, y, fill = mean)) +
## ##   geom_tile()

## ## +
## ##   geom_sf(data = nests, color = "firebrick", size = 1, pch = 4, alpha = 0.2) +
## ##   ggtitle("Nest intensity per ~1.25 ha")

## ## stk.y <- inla.stack(data = list(amount = y, # for single model
## ##                                 alldata = cbind(y, NA)), # for joint model
## ##                     A = list(A.est, 1),
## ##                     effects = list(
## ##                       list(y.field = 1:spde$n.spde),
## ##                       list(y.intercept = rep(1, length(y)))),
## ##                     tag = "est.y"
## ##                     )


## ## stk.z <- inla.stack(data = list(occurrence = z, # for single model
## ##                                 alldata = cbind(NA, z)), # for joint model
## ##                     A = list(A.est, 1),
## ##                     effects = list(
## ##                       list(z.field = 1:spde$n.spde,
## ##                            zc.field = 1:spde$n.spde),
## ##                       list(z.intercept = rep(1, length(z)))),
## ##                     tag = "est.z"
## ##                     )

## ## ## make a quiet version of INLA to parse logs
## ## #inla.q <- quietly(inla)

## ## # spatial response only



## ## ## # spatial binomial only
## ## ## f.z <- occurrence ~ -1 + z.intercept + f(z.field, model = spde)
## ## ## out.z <- inla(f.z, family = "binomial",
## ## ##               data = inla.stack.data(stk.z),
## ## ##               control.predictor = list(A = inla.stack.A(stk.z), compute = TRUE),
## ## ##               control.compute = list(dic = TRUE, return.marginals.predictor=TRUE),
## ## ##               control.inla = list(strategy = "laplace"))

## ## ## print(glue("time to fit spatial binom only is: {toc()[3]}"))

## ## #######################
## ## # spatial joint model #
## ## # the real deal       #
## ## #######################

## ## full.start.time <- toc()[3] / 60
## ## stk.yz <- inla.stack(stk.y, stk.z)
## ## f.yz <- alldata ~ -1 + y.intercept + z.intercept +
## ##   f(y.field, model = spde) + f(z.field, model = spde) +
## ##   f(zc.field, copy = "y.field", fixed = FALSE)
## ## out.yz <- inla(f.yz, family = c("gamma", "binomial"),
## ##                data = inla.stack.data(stk.yz),
## ##                control.predictor = list(A = inla.stack.A(stk.yz), compute = TRUE),
## ##                control.compute = list(dic = TRUE, config = TRUE, return.marginals.predictor=TRUE),
## ##                control.inla = list(strategy = "laplace"),
## ##                verbose = F)

## ## full.start.elapsed.time <- (glue("time to fit full model at time {time.i} is: {round((toc()[3]/60 - full.start.time) ,2)} min"))
## ## # time.mat[2, time.i] <- full.start.elapsed.time # for profile timing
## ## print(full.start.elapsed.time)

## ## ## # non-spatial joint
## ## ## f.nospatial <- alldata ~ -1 + z.intercept + y.intercept
## ## ## out.nospatial <- inla(f.nospatial, family=c("gamma", "binomial"),
## ## ##                       data=inla.stack.data(stk.yz),
## ## ##                       control.predictor=list(A=inla.stack.A(stk.yz)),
## ## ##                       control.compute=list(dic=TRUE),
## ## ##                       control.inla=list(strategy="laplace"),
## ## ##                       verbose = F)

## ## ## # joint but not shared spatial component
## ## ## f.noshared <- alldata ~ -1 + y.intercept + z.intercept +
## ## ##   f(y.field, model=spde) + f(z.field, model=spde)
## ## ## out.noshared <- inla(f.noshared, family=c("gamma", "binomial"),
## ## ##                      data=inla.stack.data(stk.yz),
## ## ##                      control.predictor=list(A=inla.stack.A(stk.yz)),
## ## ##                      control.compute=list(dic=TRUE),
## ## ##                      control.inla=list(strategy="laplace"))

## ## ##save(file = file.path(out.top.dir, out.sub.dir, "dual-likelihood-test-fit.Rdata"), list = ls())
## ## ##load(file.path(out.top.dir, "2024-02-15/dual-likelihood-test-fit.Rdata"))

## ## ## perform model comparison between the three models using DIC
## ## # note: need to be careful with the DIC when using multiple likelihoods (joint)
## ## #       must sum local DIC for each obs
## ## # we want lower DIC which indicates better fit

## ## ## TODO rewrite all numerical outputs into one output row -> csv
## ## ## write.csv(rbind(nospatial=tapply(out.nospatial$dic$local.dic, out.nospatial$dic$family, sum),
## ## ##                 noshared=tapply(out.noshared$dic$local.dic, out.noshared$dic$family, sum),
## ## ##                 separate=c(out.y$dic$dic, out.z$dic$dic),
## ## ##                 joint=tapply(out.yz$dic$local.dic, out.yz$dic$family, sum)),
## ## ##           file.path(o.d, glue("{lr.n}-DIC-comp.csv")))

## ## ## calculate the expected value of the posterior marginal
## ## ## distributions on the response scale
## ## idy <- which(z == 1)
## ## exp.y <- sapply(out.yz$marginals.linear.predictor[idy],
## ##                 function(m){inla.emarginal(exp, m)})

## ## ## TODO is this which call right? it's what is in the INLA book pg 276
## ## idz <- which(!is.na(z))
## ## ## NOTE 28aug2024 some of the marginals for z (pre transform) have Inf height
## ## ## this breaks inla.emarginal, this the tryCatch. But it's not a big deal since we
## ## ## exp.z is only used to store mean(exp.z) as a summary metric
## ## exp.z <- tryCatch({sapply(out.yz$marginals.linear.predictor[length(lr.dat) + idz],
## ##                           function(m){inla.emarginal(inla.link.invlogit, m)})
## ## },
## ## error = function(err){
## ##   NA
## ## })


## ## ## # compare these values to the observed mean for the signal occurrence
## ## ## # and amount -- a nice simple sanity check to make sure everything is
## ## ## # as expected
## ## ## means.comp <- matrix(c(yPositive.obs = mean(lr.dat[which(lr.dat > 0)]),
## ## ##                        yPositive.pred = mean(exp.y),
## ## ##                        z.obs = mean(z, na.rm = T), z.pred = mean(exp.z)),
## ## ##                      nrow = 2, byrow = T)
## ## ## colnames(means.comp) <- c("obs", "pred")
## ## ## rownames(means.comp) <- c("yPositive", "z")
## ## ## write.csv(means.comp, file.path(o.d, glue("{lr.n}-obs-est-means-comp.csv")))

## ## ## plot the posterior marginal distributions for y.int, z.int,
## ## ## tau for separate and joint models, and beta_1 for joint
## ## ## as in figure 8.10 pg 277
## ## #print("starting first plot")

## ## png(file.path(o.d, glue("{lr.n}-{time.i}-post-marg-dist.png")),
## ##     width = 8, height = 8, units = "in", res = 300)
## ## par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=2:0)
## ## plot(inla.smarginal(out.yz$marginals.fixed[[1]]), type='l', ylab='Density',
## ##      ylim=c(0,max(out.yz$marginals.fixed[[1]][,2], out.y$marginals.fixed[[1]][,2])),
## ##      xlab=expression(b[0]^y))
## ## lines(inla.smarginal(out.y$marginals.fixed[[1]]), lty=2)
## ## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## ## plot(inla.smarginal(out.yz$marginals.fixed[[2]]), type='l', ylab='Density',
## ##      #ylim=c(0,max(out.yz$marginals.fixed[[2]][,2], out.z$marginals.fixed[[1]][,2])),
## ##      xlab=expression(b[0]^z))
## ## #lines(inla.smarginal(out.z$marginals.fixed[[1]]), lty=2)
## ## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## ## plot(inla.smarginal(out.yz$marginals.hyperpar[[1]]), type='l', ylab='Density',
## ##      ylim=c(0,max(out.yz$marginals.hy[[1]][,2], out.y$marginals.hy[[1]][,2])),
## ##      xlab=expression(tau))
## ## lines(inla.smarginal(out.y$marginals.hyperpar[[1]]),lty=2)
## ## legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

## ## plot(inla.smarginal(out.yz$marginals.hyperpar[[6]]), type='l', ylab='Density',
## ##      ylim=c(0,max(out.yz$marginals.hy[[6]][,2])),  xlab=expression(beta[1]))
## ## dev.off()

## ## ## to check if the occurrence and amount of signaling are spatially
## ## ## dependent, we can test the significance of the spatial random
## ## ## effects zeta(s) and u(s)
## ## png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-effects-ordered.png")),
## ##     width = 8, height = 8, units = "in", res = 300)
## ## par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,1,0))
## ## ordxy <- intersect(order(out.yz$summary.random$y.field$mean), in.d) # in.d(omain)
## ## plot(out.yz$summary.random$y.field$mean[ordxy], type="l", ylab=expression(xi[i]), ylim=range(out.yz$summary.random$y.field[, 4:6]))
## ## for (i in c(4,6)) lines(out.yz$summary.random$y.field[ordxy,i], lty=2)
## ## abline(h=0, lty=3)
## ## ordxz <- intersect(order(out.yz$summary.random$z.field$mean), in.d)
## ## plot(out.yz$summary.random$z.field$mean[ordxz], type="l", ylab=expression(u[i]), ylim=range(out.yz$summary.random$z.field[, 4:6]))
## ## for (i in c(4,6)) lines(out.yz$summary.random$z.field[ordxz,i], lty=2)
## ## abline(h=0, lty=3)
## ## dev.off()

## ## ## #--- Significance of the spatial effects ---#
## ## spat.sig <- data.frame(n.nodes.in.d=length(in.d),
## ##                        expected.out=0.05*length(in.d),
## ##                        observed.y.out=sum(out.yz$summary.random$y.field[in.d,4]>0 |
## ##                                             out.yz$summary.random$y.field[in.d,6]<0),
## ##                        observed.z.out=sum(out.yz$summary.random$z.field[in.d,4]>0 |
## ##                                             out.yz$summary.random$z.field[in.d,6]<0))
## ## #write.csv(spat.sig, file.path(o.d, glue("{lr.n}-spat-sig.csv")), row.names = F)

## ## out.yz.summ <- data.frame(lr.n = lr.n,
## ##                           time = time.i,
## ##                           var = c("local.dic",
## ##                                   "family.dic",
## ##                                   "yPos-obs",
## ##                                   "yPos-pred",
## ##                                   "z-obs",
## ##                                   "z-pred",
## ##                                   colnames(spat.sig)
## ##                                   ),
## ##                           val = c(tapply(out.yz$dic$local.dic, out.yz$dic$family, sum),
## ##                                   # the tapply returns both local and family dic
## ##                                   mean(lr.dat[which(lr.dat > 0)]),
## ##                                   mean(exp.y),
## ##                                   mean(z, na.rm = T),
## ##                                   mean(exp.z),
## ##                                   as.numeric(spat.sig[1, ]))
## ##                           )

## ## write.csv(out.yz.summ, file.path(o.d, glue("{lr.n}-{time.i}-output-summary.csv")), row.names = F)


## ## ## extract the spatial effects
## ## y.field <- inla.spde2.result(out.yz, "y.field", spde)
## ## z.field <- inla.spde2.result(out.yz, "z.field", spde)

## ## png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-effects-marginals.png")),
## ##     width = 8, height = 8, units = "in", res = 300)
## ## par(mfrow=c(2,2))
## ## plot.default(y.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
## ## plot.default(y.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
## ## plot.default(z.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
## ## plot.default(z.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')
## ## dev.off()

## ## #### plot spatial random fields

## ## ## make projection grid

## ## # predict resolution (number of pixels in x and y)
## ## xy.pred <- make.pred.grid(xy.obs = xy.obs,
## ##                           retain.shape = T,
## ##                           res = pred.grid.res,
## ##                           boundary.poly = boundary.coords)
## ## if(j.i  == 1){
## ##   write.csv(xy.pred, file.path(o.d, glue("xy-pred-locs-{time.i}.csv")), row.names = F)
## ## }

## ## A.pred <- inla.spde.make.A(
## ##   mesh = mesh,
## ##   loc = xy.pred)

## ## rf.grid <- list(xi.mean = fm_evaluate(mesh, field = out.y$summary.random$y.field$mean, loc = xy.pred),
## ##                 xi.sd = fm_evaluate(mesh, field = out.y$summary.random$y.field$sd, loc = xy.pred),
## ##                 obs = lr.dat,
## ##                 u.mean = fm_evaluate(mesh, field = out.yz$summary.random$z.field$mean, loc = xy.pred),
## ##                 u.sd = fm_evaluate(mesh, field = out.yz$summary.random$z.field$sd, loc = xy.pred),
## ##                 occ.obs = z)

## ## xy.in <- inout(xy.pred, boundary.coords)
## ## for (i in c(1:2, 4:5)) rf.grid[[i]][!xy.in] <- NA


## ## png(file.path(o.d, glue("{lr.n}-{time.i}-spatial-fields.png")),
## ##     width = 12, height = 8, units = "in", res = 500)
## ## par(mfrow=c(2,3), mar=c(0,0,0,0))
## ## for (i in 1:6) {
## ##   if(grepl("obs", names(rf.grid)[i])){
## ##     quilt.plot(xy.obs, rf.grid[[i]],
## ##                nx = 60, ny = 50,
## ##                nlevel=100, col=viridis(100),
## ##                legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
## ##                xaxt = 'n', yaxt = 'n')
## ##   }else{
## ##     quilt.plot(xy.pred, rf.grid[[i]],
## ##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
## ##                nlevel=101, col=viridis(100),
## ##                legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0),
## ##                xaxt = 'n', yaxt = 'n')
## ##   }
## ##   lines(boundary.coords, col = "orange", lwd = 2)
## ## }
## ## dev.off()

## ## #####################################
## ## #--- Prediction of the responses ---#
## ## #####################################

## ## s1 <- inla.posterior.sample(n=1,result=out.yz)

## ## names(s1[[1]])
## ## grep("y.intercept", rownames(s1[[1]]$latent), fixed=TRUE)

## ## ids <- lapply(c('y.intercept', 'z.intercept', 'y.field', 'z.field', 'zc.field'), function(x) grep(x, rownames(s1[[1]]$latent), fixed=TRUE))
## ## pred.y.f <- function(s) exp(s$latent[ids[[1]], 1] + s$latent[ids[[3]], 1])
## ## pred.z.f <- function(s) 1/(1 + exp(-(s$latent[ids[[2]], 1] + s$latent[ids[[4]], 1] + s$latent[ids[[5]], 1]))) ## NOTE! the latent zc.field is already beta.zc*y.field

## ## s1000 <- inla.posterior.sample(1000, out.yz)

## ## prd.y <- sapply(s1000, pred.y.f)
## ## prd.z <- sapply(s1000, pred.z.f)
## ## prd.y.cond.z <- matrix(rbinom(n = length(prd.z), size = 1, p = prd.z),
## ##                        ncol = ncol(prd.z))* prd.y

## ## # project all draws to pred grid
## ## prd.y.proj <- apply(prd.y, 2,
## ##                     function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
## ## prd.z.proj <- apply(prd.z, 2,
## ##                     function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})
## ## prd.y.cond.z.proj <- apply(prd.y.cond.z, 2,
## ##                            function(s){fm_evaluate(mesh, field=s, loc = xy.pred)})


## ## # always write the means
## ## write.csv(as.matrix(rowMeans(prd.y.proj), ncol = 1),
## ##           file.path(o.d, glue("{lr.n}-{time.i}-resp-est.csv")), row.names = F)
## ## write.csv(as.matrix(rowMeans(prd.y.cond.z.proj), ncol = 1),
## ##           file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-est.csv")), row.names = F)
## ## write.csv(as.matrix(apply(prd.y.cond.z.proj, 1, sd), ncol = 1),
## ##           file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-sd.csv")), row.names = F)

## ## # only save some of the complete draws -- these are big files (.5Gb)
## ## # and only save the "final" field draws
## ## if(lr.n %in% lr.to.save.draws){
## ##   write.csv(prd.y.cond.z.proj,
## ##             file.path(o.d, glue("{lr.n}-{time.i}-cond-resp-draws.csv")), row.names = F)
## ## }

## ## # prepare to plot obs, est, sd for y, z, y|z

## ## prd <- list(resp.obs = lr.dat)
## ## prd$y.mean <- fm_evaluate(mesh, field=rowMeans(prd.y), loc = xy.pred)
## ## prd$y.sd <- fm_evaluate(mesh, field=apply(prd.y, 1, sd), loc = xy.pred)

## ## prd$occ.obs <- z
## ## prd$z.mean <- fm_evaluate(mesh, field=rowMeans(prd.z), loc = xy.pred)
## ## prd$z.sd <- fm_evaluate(mesh, field=apply(prd.z, 1, sd), loc = xy.pred)

## ## prd$cond.resp.obs <- lr.dat
## ## prd$y.cond.z.mean <- fm_evaluate(mesh, field=rowMeans(prd.y.cond.z), loc = xy.pred)
## ## prd$y.cond.z.sd <- fm_evaluate(mesh, field=apply(prd.y.cond.z, 1, sd), loc = xy.pred)

## ## # knock out places where we don't want to plot/predict
## ## for (j in c(2:3, 5:6, 8:9)) prd[[j]][!xy.in] <- NA

## ## titles <- c("Y Obs", "Y Est", "Y SD",
## ##             "Z Obs", "Z Est", "Z SD",
## ##             "Y|Z Obs", "Y|Z Est", "Y|Z SD")

## ## # *** Code for Figure 8.14
## ## png(file.path(o.d, glue("{lr.n}-{time.i}-response-occurrence-fields.png")),
## ##     width = 12, height = 12, units = "in", res = 500)
## ## par(mfcol=c(3,3), mar=c(0,0,0,1))
## ## for (i in 1:length(prd)) {

## ##   # set xlim and ylim for consistency
## ##   xl <- range(xy.obs[, 1]) + c(-.05, .05) * (max(xy.obs[, 1]) - min(xy.obs[, 1]))
## ##   yl <- range(xy.obs[, 2]) + c(-.05, .05) * (max(xy.obs[, 2]) - min(xy.obs[, 2]))

## ##   ## first, set z/color ranges and flag large predicted values in estimates
## ##   if(i %in% c(1, 2, 7, 8)){

## ##     zl <- c(0, max(prd$resp.obs) * 1.01) ## mean and data obs plots
## ##     zcols <- c(viridis(100), "#b33609")

## ##     # this should be a shrinkage/smoothing model. outputs above max observed are truncated and
## ##     # plotted as dark red
## ##     prd[[i]][prd[[i]] > max(prd$resp.obs) * 1.01] <- max(prd$resp.obs) * 1.01

## ##   }else if(i %in% c(4, 5)){## prob of non-zero response mean/est and obs

## ##     zl <- 0:1
## ##     if(i == 4) zcols <- viridis(2)
## ##     if(i == 5) zcols <- viridis(100)

## ##   }else{ ## sd plots, no scaling

## ##     zl <- c(0, max(prd[[i]], na.rm = T))
## ##     zcols <- turbo(100)

## ##   }

## ##   ## second, plot the model outputs differently than the raw data

## ##   if(grepl("obs", names(prd)[i])){
## ##     bubblePlot(xy.obs[, 1], xy.obs[, 2], prd[[i]],
## ##                col=zcols,
## ##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
## ##                legend.args=list(text="", line=0),
## ##                xaxs = 'i', yaxs = 'i', # o.w. R adds 4% to pad axes and
## ##                # then plot this won't align with quilt.plot
## ##                xaxt = "n", yaxt = "n", bty = "n",
## ##                xlim = xl, ylim = yl, zlim = zl)
## ##   }else{
## ##     quilt.plot(xy.pred, prd[[i]],
## ##                nx = pred.grid.res * .9, ny = pred.grid.res * .9,
## ##                nlevel=length(zcols), col=zcols,
## ##                legend.mar=7, #legend.args=list(text=names(prd)[i], line=0),
## ##                legend.args=list(text="", line=0),
## ##                xaxt = 'n', yaxt = 'n', bty = "n",
## ##                xlim = xl, ylim = yl, zlim = zl)
## ##   }

## ##   ## add title and boundary lines
## ##   text(max(xy.obs[, 1]), max(xy.obs[, 2]), adj = c(1, 0.5), titles[i], cex = 2)
## ##   lines(boundary.coords, lwd = 2)
## ## }
## ## dev.off()


## ## ## png(file.path(o.d, "test.png"), width = 5, height = 5, units = "in", res = 300)
## ## ## DefaultAssay(brain)<- "NeighborhoodToCell_ALRA"
## ## ## SpatialFeaturePlot(brain, features = lr.n, slot = "data") +
## ## ##   scale_fill_viridis_c(limits = c(0.2, 1))
## ## ## dev.off()

## ## final.elapsed.time <- (glue("final elapsed time at time {time.i} : {round(toc()[3]/60,2)} min"))
## ## # time.mat[3, time.i] <- final.elapsed.time # for profile timing
## ## print(final.elapsed.time)
