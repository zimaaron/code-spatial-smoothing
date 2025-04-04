
## combine all the mean estimates from an output dir
setwd("~/Dropbox/genetics/a-s-omics/")
i.d <- file.path(getwd(), "/data-inputs/mouse-embryo-raw")
o.d <- "~/Dropbox/genetics/a-s-omics//data-outputs/mouse-embryo-raw/2025-01-27" # initial output location

require(data.table)
require(deSolve)
require(ReacTran)
require(glue)
require(fields)

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

#########################################
## perform diffusion and sequestration ##
#########################################

## make buffer around spatial domain to minimize edge effects

# load one feature and grab xy coords
xy <- readRDS(p.f.p[1])$total[, .(x, y)]
# append on border of 20 extra cells
buff.bin.ct = 20
bin.size = 50
xy.b.buff <- data.table(x = rep(sort(unique(xy[, x])), each = buff.bin.ct),
                        y = rep((min(xy[, y]) - (1:buff.bin.ct) * bin.size),
                                length(unique(xy[, x])))) # bottom
xy.t.buff <- data.table(x = rep(sort(unique(xy[, x])), each = buff.bin.ct),
                        y = rep((max(xy[, y]) + (1:buff.bin.ct) * bin.size),
                                length(unique(xy[, x])))) # top
xy.w.buff <- rbind(xy, xy.b.buff, xy.t.buff)
xy.l.buff <- data.table(x = rep((min(xy.w.buff[, x]) - (1:buff.bin.ct) * bin.size),
                                length(unique(xy.w.buff[, y]))),
                        y = rep(sort(unique(xy.w.buff[, y])), each = buff.bin.ct)) # left
xy.r.buff <- data.table(x = rep((max(xy.w.buff[, x]) + (1:buff.bin.ct) * bin.size),
                                length(unique(xy.w.buff[, y]))),
                        y = rep(sort(unique(xy.w.buff[, y])), each = buff.bin.ct)) # right
xy.w.buff <- rbind(xy.w.buff, xy.l.buff, xy.r.buff)
xy.w.buff.map <- copy(xy.w.buff)
xy.w.buff.map[, domain := 0]
xy.w.buff.map[1:nrow(xy), domain := 1]
# sort into L->R, from T->B, as R fills in a matrix from a vec
# then map how the ordering of the original xy maps into the xy.w.buff
# in this vector to matrix ordering
xy.w.buff.map[, neg.y := -y]
setkey(xy.w.buff.map, neg.y, x)
xy.w.buff.map[, xy.idx := 1:.N]
xy.idx <- xy.w.buff.map[domain == 1, .(x, y, xy.idx)]
xy.y.map <- merge(xy, xy.idx)
y.xy.idx <- xy.y.map[, xy.idx]

# due to edge effects, we also plot a bit zoomed in
uniq.x <- xy[, 1] |> unlist() |> unique() |> sort()
uniq.y <- xy[, 2] |> unlist() |> unique() |> sort()
x.plot.lim <- uniq.x[c(5, length(uniq.x) - 4)]
y.plot.lim <- uniq.y[c(5, length(uniq.y) - 4)]

## function that performs diffusion and sequestration calculations
Diff2D <- function (t, # time
                    y, # ligand concentration
                    l.prod, # ligand production
                    r.conc, # receptor current concentration vector
                    produce.l, # T/F, should ligand be produced?
                    sequest.l, # T/F, should sequestered ligand be subtracted?
                    in.lig.dt, # initial values
                    pp.lig.dt, # post-production values
                    ps.lig.dt, # post-sequestration values
                    r.prod, # receptor concentration production
                    produce.r, # T/F, should receptor be produced ?
                    sequest.r, # T/F, should sequestered receptor be subtracted?
                    # each time step? if F, then the
                    # initial amount decays over time until it's gone
                    in.rec.dt, # initial values
                    pp.rec.dt, # post-production values
                    ps.rec.dt, # post-sequestration values
                    seques.dt, # pre-initialized data.table to store sequestration
                    lderiv.dt, # pre-initialized data.table to store lig derivs
                    seq.before.dif, #  T/F should sequestration happen before diffusion?
                    parms,
                    Nx, Ny, # dimensions of the concentration matrix
                    Dx = 1, Dy = 1, # diffusion coef in x/y directions
                    dx, dy, # distince between adjacent cell interfaces in x/y
                    zero.boundary = F # T/F for forcing 0 concentration at boundary
                    ) {

  # IMPORTANT NOTE :
  # since I don't know how to save the intermediate outputs from this
  # function when running it in ode.2d
  # we leverage data.table pointers to do this

  # one option for sequestration calculations
  elemwise.min <- function(x, y){
    both <- cbind(as.vector(x), as.vector(y))
    elem.min <- apply(both, 1, min)
    matrix(elem.min, nxlim = x.plot.lim, ylim = y.plot.lim, col = ncol(x), nrow = nrow(x))
  }

  # if seq.before.dif == T
  # 1) produce lig/rec
  # 2) sequesteration
  # 3) ligand diffusion
  # this is simpler in the ode framework, since diffusion happens "in between" time steps

  # if seq.before.dif == F
  # 1) produce lig/rec
  # 2) ligand diffusion
  # 3) sequestration
  # this is more complicated since diffusion now happens "in the middle of" each time step

  # cast lig and rec concentration into 2d domain
  CONC.l <- matrix(nrow = Nx, nxlim = x.plot.lim, ylim = y.plot.lim, col = Ny, y)
  CONC.r <- matrix(nrow = Nx, nxlim = x.plot.lim, ylim = y.plot.lim, col = Ny, r.conc)

  # store values of l and r as they come into the function
  in.lig.dt[time == t, lig := as.vector(CONC.l)]
  in.rec.dt[time == t, rec := as.vector(CONC.r)]

  if(seq.before.dif){

    ##################
    # seq.before.dif #
    ##################

    ## 1) lig and rec production
    if(produce.l){
      CONC.l <- CONC.l + l.prod
    }
    if(produce.r){
      CONC.r <- CONC.r + r.prod
    }

    # store post-production values of l and r
    pp.lig.dt[time == t, lig := as.vector(CONC.l)]
    pp.rec.dt[time == t, rec := as.vector(CONC.r)]

    ## 2) perform sequestration

    ## TODO, replace min with any other binding function
    seq.mat <- elemwise.min(CONC.l, CONC.r)

    if(sequest.l){
      CONC.l <- CONC.l - seq.mat
    }
    if(sequest.r){
      CONC.r <- CONC.r - seq.mat
    }

    # save post-sequestration concentrations and the sequestration values
    ps.lig.dt[time == t, lig := as.vector(CONC.l)]
    ps.rec.dt[time == t, rec := as.vector(CONC.r)]
    seques.dt[time == t, seq := as.vector(seq.mat)]

  }else{ # if/else seq.before.diff

    ##################
    # dif.before.seq #
    ##################

    # since diffusion is the ODE, this must happen "in between times"
    # so we juggle time 0 a bit to allow for this:

    # on time 0, we produce and diffuse
    # on all other times, we sequester, produce, diffuse

    if(t > 0){

      ## 1) perform sequestration

      ## TODO, replace min with any other binding function
      seq.mat <- elemwise.min(CONC.l, CONC.r)
    }else{

      # at time 0, just make a 0 sequestering matrix
      seq.mat <- matrix(0, nxlim = x.plot.lim, ylim = y.plot.lim, col = ncol(CONC.l), nrow = nrow(CONC.l))
    }

    if(sequest.l){
      CONC.l <- CONC.l - seq.mat
    }
    if(sequest.r){
      CONC.r <- CONC.r - seq.mat
    }

    ## save receptor and sequestration to global environ
    # since I don't know how to save the intermediate outputs from this
    # function when running it in ode.2d
    # i leverage that data.table uses pointers to do this

    # save post-sequestration concentrations and the sequestration values
    ps.lig.dt[time == t, lig := as.vector(CONC.l)]
    ps.rec.dt[time == t, rec := as.vector(CONC.r)]
    seques.dt[time == t, seq := as.vector(seq.mat)]


    ## 2) lig and rec production
    if(produce.l){
      CONC.l <- CONC.l + l.prod
    }
    if(produce.r){
      CONC.r <- CONC.r + r.prod
    }

    # store post-production values of l and r
    pp.lig.dt[time == t, lig := as.vector(CONC.l)]
    pp.rec.dt[time == t, rec := as.vector(CONC.r)]

  }

  ## 3) ligand diffusion/transport
  if(zero.boundary){
    Tran <- tran.2D(C = CONC.l,
                    D.x = Dx, D.y = Dy, # diffusion coefficient in x/y directions
                    dx = dx, dy = dy, # distance between adjacent cell interfaces in x/y dirdx                    C.x.up = 0, C.x.down = 0, C.y.up = 0, C.y.down = 0 # boundary concentrations
                    )
  }else{
    Tran <- tran.2D(C = CONC.l,
                    D.x = Dx, D.y = Dy, # diffusion coefficient in x/y directions
                    dx = 50, dy = dy # distance between adjacent cell interfaces in x/y dirs
                    )
  }

  # transport + reaction
  lderiv.dt[time == t, dl := as.vector(Tran$dC)]
  dCONC <- as.vector(Tran$dC) # + r*CONC)

  # examples of other ways to alter the ODE system:
  # Bioirrigation in a central spot
  # mid <- N/2
  # dCONC[mid, mid] <- dCONC[mid, mid] + irr*(Bc - CONC[mid, mid])

  return (list(dCONC, CONC.r))
}

## specify diffusivity and distance between cells, in micrometers, um, 1e-6 meters

# 50 spots binned, 500nm between spots, 500e-9*50 = 25e-6 == 25um
dx <- dy <- 25

# some numbers from MSR
diff.dict <- data.table(feat = lr.feats, diffusivity = 2) # default, in um^2/sec
diff.dict[feat == "IGF1", diffusivity := 1.59e-6 * 1000 ^ 2] # cm^2 /  sec -> um^2/sec
diff.dict[feat == "FGF2", diffusivity := 1.32e-4 * 1000 ^ 2/ 60] # cm^2 / min -> um^2/sec

## set up diffusion experiment  parameters
exp.specs <- data.table(expand.grid(time.step = c(1, 5, 10), # in second
                                    # if T, then R is constant every timestep
                                    # ow, then R produces and is sequestered at each time
                                    static.recept = c(T, F),
                                    # T/F should sequestration happen before diffusion?
                                    seq.then.dif = c(T, F)))

## loop through and run diffusion experiments
for(ii in 1:exp.specs[, .N]){

  time.step <- exp.specs[ii, time.step]
  static.r <- exp.specs[ii, static.recept]
  seq.then.dif <- exp.specs[ii, seq.then.dif]

  experiment.name <- glue("{ifelse(seq.then.dif,'seq-then-diff','diff-then-seq')}-{ifelse(exp.specs[ii,static.recept],'static-recept','dynamic-receptor')}-timestep-{time.step}")
  s.d <- file.path(o.d, experiment.name) # save directory for diffusion outputs
  dir.create(s.d)

  # loop through the mechanism preds and process all for eac experiment
  for(pred.fp in p.f.p){

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

    diffusivity = diff.dict[feat == l.n, diffusivity]

    cat(glue("on experiment {ii} of {exp.specs[,.N]} || mech {which(p.f.p==pred.fp)} of {length(p.f.p)} || {m.n}\n"))
    cat("\n")

    ## set up args and run diffusion
    ## final timestep will be final value minus 1
    ode.times <- seq(0, 101 * time.step, time.step)

    qp.res.x <- pred.t[, length(unique(x))]
    qp.res.y <- pred.t[, length(unique(y))]

    lig.prod <- rec.prod <- CONC.r <- numeric(nrow(xy.w.buff))
    lig.prod[y.xy.idx] <- pred.l$mean
    rec.prod[y.xy.idx] <- pred.r$mean
    if(static.r){
      CONC.r <- rec.prod
    }else{
      CONC.r <- rep(0, length(rec.prod))
    }
    init.rec.dt <- data.table(time = rep(ode.times, each = length(lig.prod)),
                              xy.idx = rep(1:length(lig.prod), times = length(ode.times)),
                              rec = -1)
    post.pro.rec.dt <- copy(init.rec.dt) # post production
    post.seq.rec.dt <- copy(init.rec.dt) # post seq
    seq.dt <- data.table(time = rep(ode.times, each = length(lig.prod)),
                         xy.idx = rep(1:length(lig.prod), times = length(ode.times)),
                         seq = -1)
    der.dt <-  data.table(time = rep(ode.times, each = length(lig.prod)),
                          xy.idx = rep(1:length(lig.prod), times = length(ode.times)),
                          dl = -1) # ligand derivative

    init.lig.dt <- data.table(time = rep(ode.times, each = length(lig.prod)),
                              xy.idx = rep(1:length(lig.prod), times = length(ode.times)),
                              lig = -1)
    post.pro.lig.dt <- copy(init.lig.dt) # post production
    post.seq.lig.dt <- copy(init.lig.dt) # post seq

    # run the diffusion ODE in time
    ode.out <- ode.2D(func = Diff2D, y = rep(0, length(rec.prod)), times = ode.times,
                      # ligand args
                      l.prod = lig.prod,
                      produce.l = T,
                      sequest.l = T,
                      in.lig.dt = init.lig.dt, # dt to store ligand init concentration
                      pp.lig.dt = post.pro.lig.dt, # dt to store ligand post-prod concentration
                      ps.lig.dt = post.seq.lig.dt, # dt to store ligand post-seq concentration
                      # recep args
                      r.conc = CONC.r,
                      r.prod = rec.prod,
                      produce.r = !static.r,
                      sequest.r = !static.r,
                      in.rec.dt = init.rec.dt, # dt to store receptor init concentration
                      pp.rec.dt = post.pro.rec.dt, # dt to store receptor post-prod concentration
                      ps.rec.dt = post.seq.rec.dt, # dt to store receptor post-seq concentration
                      # other values to store at each iteration
                      seques.dt = seq.dt, # dt to store sequestration
                      lderiv.dt = der.dt, # dt to store ligand diffusion derivative
                      # domain and diffusion args
                      seq.before.dif = seq.then.dif,
                      Nx = length(unique(xy.w.buff[, x])), Ny = length(unique(xy.w.buff[, y])),
                      Dx = diffusivity, Dy = diffusivity,
                      dx = dx, dy = dy,
                      zero.boundary = F,
                      parms = NULL, lrw = 10000000,
                      dimens = c(length(unique(xy.w.buff[, x])), length(unique(xy.w.buff[, y]))))
    ## plot diffusion/sequestration results

    # get the initial and final state (at second to last time) for first page
    l.conc.init    <- post.pro.lig.dt[time == ode.times[1], lig] # L starts as 0, get post-prod val
    r.conc.init    <- post.pro.rec.dt[time == ode.times[1], rec] # R starts as 0, get post-prod val
    seq.rate.init  <- seq.dt[time == ode.times[1], seq]
    if(!seq.then.dif){
      # sequestration happens in next step
      seq.rate.init  <- seq.dt[time == ode.times[1 + 1], seq]
    }
    seq.total.init <- seq.dt[time <= ode.times[1], sum(seq), by = xy.idx][, V1]

    l.conc.final    <- post.pro.lig.dt[time == tail(ode.times, 2)[1], lig]
    r.conc.final    <- post.pro.rec.dt[time == tail(ode.times, 2)[1], rec]
    seq.rate.final  <- seq.dt[time == tail(ode.times, 2)[1], seq]
    seq.total.final <- seq.dt[time <= tail(ode.times, 2)[1], sum(seq), by = xy.idx][, V1]

    # initialize pdf to store multiple plot pages
    pdf(file = file.path(s.d, glue('{experiment.name}-{l.n}-{r.n}-01-timestep.pdf')), onefile = T,
        width = (qp.res.x / qp.res.y) * 13 + 17.5, height = 13)

    # setup plot panel
    par(mfrow = c(3, 5), mai = c(.62, 0.82, .62, 1.22))

    ## set color palettes

    # diverging, for derivative
    brbg <- colorspace::diverging_hcl(64, h = c(180, 50),
                                      c = 80, l = c(20, 95),
                                      power = c(0.7, 1.3))
    viri <- viridis(64) # lig
    plas <- plasma(64)  # rec
    mako <- mako(64)    # seq

    ## plot first page summary

    lig.lim <- range(c(l.conc.init, l.conc.final))
    rec.lim <- range(c(r.conc.init, r.conc.final))
    seq.lim <- range(c(seq.rate.init, seq.rate.final))

    if(min(lig.lim) < 0) lig.lim[1] <- 0
    if(min(rec.lim) < 0) rec.lim[1] <- 0
    if(min(seq.lim) < 0) seq.lim[1] <- 0

    # row 1
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.init,
                              main = glue("{l.n} count at time {ode.times[1]}"),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                              zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.final,
                              main = glue("{l.n} count at time {tail(ode.times, 2)[1]}"),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                              zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    # row 2
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.init,
                              main = glue("{r.n} count at time {ode.times[1]}", ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                              zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.final,
                              main = glue("{r.n} count at time {tail(ode.times, 2)[1]}", ),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                              zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    # row 3
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], seq.rate.init,
                              main = glue("{m.n} sequestration count at time {ode.times[1]}"),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                              zlim = seq.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = mako)
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], seq.rate.final,
                              main = glue("{m.n} sequestration count at time {tail(ode.times, 2)[1]}"),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                              zlim = seq.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = mako)
    fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], seq.total.init,
                              main = glue("{m.n} cumulative sequestration at time {tail(ode.times, 2)[1]}"),
                              nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x, xlim = x.plot.lim, ylim = y.plot.lim, col = mako)
    plot(0,type='n',axes=FALSE,ann=FALSE)
    plot(0,type='n',axes=FALSE,ann=FALSE)

    for(tt in ode.times[-length(ode.times)]){
      # regardless of time step, only plot on integer times (to avoid
      # thosand page pdfs that mostly look the same ...)

      # grab values to plot for this timestep

      l.conc.in <- init.lig.dt[time == tt, lig]
      l.conc.pp <- post.pro.lig.dt[time == tt, lig]
      l.conc.ps <- post.seq.lig.dt[time == tt, lig]
      l.conc.dl <- der.dt[time == tt, dl]

      r.conc.in <- init.rec.dt[time == tt, rec]
      r.conc.pp <- post.pro.rec.dt[time == tt, rec]
      r.conc.ps <- post.seq.rec.dt[time == tt, rec]

      seq.rate  <- seq.dt[time == tt, seq]
      seq.total <- seq.dt[time <= tt, sum(seq), by = xy.idx][, V1]

      # get some values from the next time step, if possible

      t.idx <- which(ode.times == tt)
      tt1 <- ode.times[t.idx + 1]

      if(t.idx <= (length(ode.times) - 1)){
        l.conc.pd  <- init.lig.dt[time == tt1, lig] # post diffusion
        l.conc.ps1 <- post.seq.lig.dt[time == tt1, lig] # next sequestration
        r.conc.pd  <- init.rec.dt[time == tt1, rec] # post diffusion
        r.conc.ps1 <- post.seq.rec.dt[time == tt1, rec] # next sequestration
      }else{
        l.conc.pd  <- rep(NA, length(l.conc.in))
        l.conc.ps1 <- rep(NA, length(l.conc.in))
        r.conc.pd  <- rep(NA, length(r.conc.in))
        r.conc.ps1 <- rep(NA, length(r.conc.in))
      }

      if(seq.then.dif){
        # since sequestration happens in time i+1, we relabel some things
        # | time i-1 | ----- time i ------- | -------- time i+1--------- | ...
        # |   function at time i    |   function at time i+1   | function at time i+2 | ...
        # | seq.i-1  * prod.i dif.i | seq.i * prod.i+1 dif.i+1 | seq.i+1 * ...
        if(tt > 0){
          # post-sequestration is really where the next time starts
          l.conc.in <- l.conc.ps
          r.conc.in <- r.conc.ps
        }
        # and grab sequestration from next time for the current one
        l.conc.ps <- l.conc.ps1
        r.conc.ps <- r.conc.ps1
      }

      # get color limits

      lig.lim <- range(c(l.conc.in, l.conc.pp, l.conc.ps, l.conc.pd, l.conc.ps1), na.rm = T)
      rec.lim <- range(c(r.conc.in, r.conc.pp, r.conc.ps, r.conc.pd), na.rm = T)
      #seq.lim <- range(c(seq.rate.init, seq.rate.final, seq.rate), na.rm = T)

      ## plot top row: lig

      fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.in,
                                main = glue("Init {l.n} count at time {tt}"),
                                nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
      fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.pp,
                                main = glue("Post Prod {l.n} count at time {tt}", ),
                                nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
      if(seq.then.dif){
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.ps,
                                  main = glue("Post Seq {l.n} count at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.dl,
                                  main = glue("d({l.n})/dt at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  xlim = x.plot.lim, ylim = y.plot.lim, col = brbg)
        # post-diff at time i is the same as init at time i+1
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.pd,
                                  main = glue("Post Diff {l.n} count at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
      }else{
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.dl,
                                  main = glue("d({l.n})/dt at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  xlim = x.plot.lim, ylim = y.plot.lim, col = brbg)
        # post-diff at time i is the same as init at time i+1
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.pd,
                                  main = glue("Post Diff {l.n} count at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], l.conc.ps1,
                                  main = glue("Post Seq {l.n} count at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = lig.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = viri)
      }

      ## middle row: rec
      fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.in,
                                main = glue("Init {r.n} count at time {tt}"),
                                nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
      fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.pp,
                                main = glue("Post Prod {r.n} count at time {tt}", ),
                                nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
      if(seq.then.dif){
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.ps,
                                  main = glue("Post Seq {r.n} count at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
        plot(0,type='n',axes=FALSE,ann=FALSE) # no derivative
        plot(0,type='n',axes=FALSE,ann=FALSE) # and no change post-diffusion
      }else{
        plot(0,type='n',axes=FALSE,ann=FALSE) # no derivative
        # this is duplicated, but it allows alignment of the fields that go into sequestration
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.pp,
                                  main = glue("Post Prod {r.n} count at time {tt}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
        fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], r.conc.ps1,
                                  main = glue("Post Seq {r.n} count at time {tt+1}", ),
                                  nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                  zlim = rec.lim, xlim = x.plot.lim, ylim = y.plot.lim, col = plas)
      }

      ## bottom row: seq
      plot(0,type='n',axes=FALSE,ann=FALSE)
      if(!seq.then.dif){
        plot(0,type='n',axes=FALSE,ann=FALSE)
        plot(0,type='n',axes=FALSE,ann=FALSE)
      }
      fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], seq.rate,
                                main = glue("{m.n} sequestration count at time {tt}"),
                                nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                xlim = x.plot.lim, ylim = y.plot.lim, col = mako)
      fields.style();quilt.plot(x = xy.w.buff.map[, x], y = xy.w.buff.map[, y], seq.total,
                                main = glue("{m.n} cumulative sequestration by time {tt}"),
                                nx = qp.res.x, ny = qp.res.y, asp = qp.res.y / qp.res.x,
                                xlim = x.plot.lim, ylim = y.plot.lim, col = mako)
      if(seq.then.dif){
        plot(0,type='n',axes=FALSE,ann=FALSE)
        plot(0,type='n',axes=FALSE,ann=FALSE)
      }

    }
    dev.off()
  }
}






# save output
# write.csv(a.r, file.path(o.d, "0000-all-smooth-ests.csv"), row.names = F)
