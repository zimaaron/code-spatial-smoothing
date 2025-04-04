# 2024-dec-6
# aoz

# This script loads in the raw raw raw counts, subsets, aggregates to 10, 25, and 50 unit squares, and saves the output for later modeling.

require(glue)
require(ggplot2)
require(reshape2)
require(Matrix)
require(Seurat)
require(dplyr)
require(data.table)
require(splancs)
require(fmesher)

setwd("~/Dropbox/genetics/a-s-omics/data-inputs/mouse-embryo-raw/")

run.entire.preprocess <- FALSE

##############################################################
#### Load data (raw counts, unbinned) -- VERY LARGE FILE ~
##############################################################
## this is a data frame with columns: gene, x, y, and value
if(run.entire.preprocess){
  full.raw <- read.delim("./full-raw-t16-data.tsv", comment.char="#") |> as.data.table()
  setnames(full.raw, 'geneID', 'feat')


  ## # Number of spots?
  ## nrow(E16.5_E1S3_GEM_bin1) # 552290433 measurements
  ## # Range of X and Y?
  ## range(E16.5_E1S3_GEM_bin1$x) # 0 44099
  ## range(E16.5_E1S3_GEM_bin1$y) # 0 26459

  # Limit to just ligand and receptor genes
  lr.names <- readRDS('ligand-receptor-names.RDs')
  setkey(full.raw, feat)
  sub <- subset(full.raw, feat %in% lr.names)
  nrow(sub) # 28044528
}

# define boundary (found visually, see code below)
boundary <- matrix(ncol = 2, data = c(20000, 9500,
                                      19000, 9500,
                                      17000, 9600,
                                      13250, 10000,
                                      10500, 11000,
                                      8200, 13000,
                                      6250, 16000,
                                      6000, 18500,
                                      6800, 20800,
                                      7300, 22300,
                                      8700, 23400,
                                      10500, 24400,
                                      16000, 25500,
                                      19120, 22800,
                                      18900, 19200,
                                      20400, 19100,
                                      22000, 22000,
                                      24200, 22200,
                                      27800, 21800,
                                      29300, 22800,
                                      30600, 22600,
                                      31700, 21500,
                                      32200, 20800,
                                      32800, 19300,
                                      32800, 17500,
                                      32300, 15300,
                                      31400, 13700,
                                      29600, 11700,
                                      28500, 10800,
                                      27000, 9800,
                                      26500, 9600,
                                      25500, 9500,
                                      20000, 9500), byrow = T)

if(run.entire.preprocess){
  feat.in.d <- which(splancs::inout(as.points(sub[, .(x, y)]), boundary))
  sub.bounded <- sub[feat.in.d, ]
  dim(sub.bounded) # 23660265, 15% reduction
  fwrite(sub.bounded, file.path('partial-raw-t16-data.csv'))

  total.cts <- full.raw[, list(total.count = sum(MIDCount)), by = .(x, y)]
  total.in.d <- which(splancs::inout(as.points(total.cts[, .(x, y)]), boundary))
  total.cts <- total.cts[total.in.d, ]
  fwrite(total.cts, file.path('total-raw-t16-data.csv'))
}

## binning the data
for(agg.size in c(25, 50, 75, 100)){

  ## load in clean data
  sub.bounded <- fread('./partial-raw-t16-data.csv')
  total.cts <- fread(file.path('total-raw-t16-data.csv'))

  ## define aggregation bins

  x.breaks = seq(6000, 32800, by = agg.size)
  x.break.val = data.table(xbin = 1:(length(x.breaks) - 1),
                           x = x.breaks[-length(x.breaks)] + agg.size / 2)
  y.breaks = seq(9500, 25500, by = agg.size)
  y.break.val = data.table(ybin = 1:(length(y.breaks) - 1),
                           y = y.breaks[-length(y.breaks)] + agg.size / 2)

  ## aggregate feature

  print(glue('aggregating feature to agg.size {agg.size}'))

  sub.bounded[, xbin := cut(x, breaks = x.breaks, labels = FALSE)]
  sub.bounded[, ybin := cut(y, breaks = y.breaks, labels = FALSE)]
  all.dat <- sub.bounded[, list(feat.count = sum(MIDCount)),
                         by = .(xbin, ybin, feat)]
  # merge on x and y vals that are the center of the bins
  all.dat <- merge(all.dat, x.break.val,
                   by = 'xbin', all.x = T, all.y = F)
  all.dat <- merge(all.dat, y.break.val,
                   by = 'ybin', all.x = T, all.y = F)
  all.dat <- na.omit(all.dat)
  # again, drop point in case centers of bins are now outside boundary
  feat.in.d <- which(splancs::inout(as.points(all.dat[, .(x, y)]), boundary))
  all.dat <- all.dat[feat.in.d, ]

  ## aggregate total

  print(glue('aggregating total to agg.size {agg.size}'))

  total.cts[, xbin := cut(x, breaks = x.breaks, labels = FALSE)]
  total.cts[, ybin := cut(y, breaks = y.breaks, labels = FALSE)]
  total.cts <- total.cts[, list(total.count = sum(total.count)),
                         by = .(xbin, ybin)]
  # merge on x and y vals that are the center of the bins
  total.cts <- merge(total.cts, x.break.val,
                     by = 'xbin', all.x = T, all.y = F)
  total.cts <- merge(total.cts, y.break.val,
                     by = 'ybin', all.x = T, all.y = F)
  total.cts <- na.omit(total.cts)
  # again, drop point in case centers of bins are now outside boundary
  total.in.d <- which(splancs::inout(as.points(total.cts[, .(x, y)]), boundary))
  total.cts <- total.cts[total.in.d, ]

  # calculate the totals per (x,y bins)
  all.dat <- all.dat[, .(x, y, feat, feat.count)]
  all.dat[, feat.observed := T]
  total.cts <- total.cts[, .(x, y, total.count)]
  total.cts[, total.observed := T]

  # add in zeros in unobserved bins since this raw dataset only
  # include locations where at least one feature measured non-zero

  # get all possible x,y locs
  all.xy <- expand.grid(x = x.break.val[, x],
                        y = y.break.val[, y]) |> data.table()
  interior.xy <- which(splancs::inout(as.points(all.xy[, .(x, y)]), boundary))
  all.xy <- all.xy[interior.xy, ]
  setkey(all.dat, feat, x, y)
  setkey(all.xy, x, y)

  # per gene, append zeros of (x,y) is missing
  zeros.dat.list <- list()
  list.idx <- 1
  # turn off warning due to data.table warning bug in `:=` call after setdiff
  options(warn = -1)
  for(gg in all.dat[, unique(feat)]){
    if(list.idx == 1 | (list.idx %% 100) == 0){
      print(glue('-- on gene {list.idx}'))
    }
    tmp <- subset(all.dat, feat == gg)
    tmp.xy <- tmp[, .(x, y)]
    missing.dat <- setdiff(all.xy, tmp.xy)
    missing.dat[, `:=` (feat = gg,
                        feat.count = 0)]
    if(missing.dat[, .N] > 0){
      zeros.dat.list[[list.idx]] <- missing.dat
      list.idx <- list.idx + 1
    }
  }
  options(warn = 0) # turn warnings back on
  zeros.dat <- rbindlist(zeros.dat.list)
  all.dat <- rbind(all.dat, zeros.dat[, feat.observed := FALSE])
  # sanity check, should be true
  # nrow(all.dat) == nrow(all.xy) * all.dat[, length(unique(feat))]
  setkey(all.dat, feat, x, y)

  # merge on totals
  all.dat <- merge(all.dat, total.cts,
                   by = c('x', 'y'),
                   all.x = T)
  # identify rows that had no data observed, fill in total and total.observed
  all.dat[which(is.na(total.observed == TRUE)), `:=` (total.count = 0,
                                                      total.observed = F)]
  all.dat[feat.count > 0, feat.present := TRUE]
  all.dat[feat.count == 0, feat.present := FALSE]
  all.dat[total.count > 0, total.present := TRUE]
  all.dat[total.count == 0, total.present := FALSE]

  ## make a few other objects while we're here:
  # 1) xy.obs
  # 2) boundary
  # 3) mesh at a few different resolutions

  xy.obs <- all.xy

  ## make mesh
  coords <- as.matrix(all.xy)

  # concaveman returns coords in CLOCKWISE fashion
  # INLA interprets interior to be defined relative to COUTERCLOCKWISE order!
  boundary.coords <- boundary
  domain <- boundary.coords[nrow(boundary.coords):1, ] |> fmesher::fm_segm()

  # make mesh
  # make mesh
  max.edge = diff(range(coords[,1]))/(3*5)
  bound.outer = diff(range(coords[,1]))/3
  mesh <- fm_mesh_2d(
    # coords,
    ## using coords creates more mesh at high resolution beyond the domain of the data
    boundary = list(fm_segm(c(domain))),
    ## using the boundary ensures that the coarser resolution starts
    ## right outside the data locations. this is likely more prone to
    ## edge effects
    max.edge = c(.1, 2) * max.edge,
    offset = c(max.edge, bound.outer / 1.5),
    cutoff = max.edge / 10,
    min.angle = 21)

  #  mesh.p <- ggplot() + inlabru::gg(mesh) + ggtitle(glue("{mesh$n} mesh nodes for {nrow(coords)} data locations"))
  #  ggsave(mesh.p, filename = file.path(o.d, "mesh-plot.png"))

  # check which mesh vertices are inside the domain
  ds1 <- splancs::as.points(boundary.coords)
  #  ds1.poly <- ds1[chull(ds1),]
  in.dom <- which(splancs::inout(mesh$loc[, 1:2], ds1))

  save("xy.obs", "coords", "domain", "in.dom", "mesh", "boundary.coords", "all.dat",
       file = file.path(glue("pre-fit-obj-binned-{agg.size}.Rdata")))
  fwrite(na.omit(all.dat), glue('ligand-receptor-binned-{agg.size}.csv'))
}

#########################################################
#### code to pick LR pairs with relatively high counts ##
#########################################################

## to get counts matrix:
# l.o@assays$RNA@layers$counts
## to get counts for a feature name:
# l.o@assays$RNA@layers$counts[which(l.n==<ligand>), ]

load(file.path(i.d, "ligand-only.Robj"))
l.o <- ligands.only
l.n <- l.o@assays$RNA@features |> rownames()
load(file.path(i.d, "recept-only.Robj"))
r.o <- receptors.only
r.n <- r.o@assays$RNA@features |> rownames()

# check the percent of non-zero counts per feature
l.nonzero <- rowMeans(l.o@assays$RNA@layers$counts > 0)
r.nonzero <- rowMeans(r.o@assays$RNA@layers$counts > 0)
#par(mfrow = c(1, 2));hist(l.nonzero);hist(r.nonzero)

## look for pairs in the LR hierarchy with both L and R having a (relatively) high percent of non-zero counts
load(file.path(i.d, "fantom-hierarchy.Robj"))
hier <- as.data.table(output)
lr.pairs <- 0
top.ct <- 5
n.pairs <- 100 # how many to look for?
while(lr.pairs < n.pairs & top.ct < 250){
  top.ct <- top.ct + 1

  # get top.ct densest Ls and Rs
  l.dense <- l.n[order(l.nonzero, decreasing = T)] |> head(, n = top.ct)
  r.dense <- r.n[order(r.nonzero, decreasing = T)] |> head(, n = top.ct)

  # see how many LR pairs we get using the top.ct densest of each
  lr.pairs <-   hier[LIGAND %in% toupper(l.dense) & RECEPTOR %in% toupper(r.dense), .N]
}


## ## find LR pairs where both are 'dense'
# use capwords to align upper/lower from the hiearchy with the raw data names
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
  {s <- substring(s, 2); if(strict) tolower(s) else s},
  sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}
lr.pairs <- hier[LIGAND %in% toupper(l.dense) & RECEPTOR %in% toupper(r.dense),
                 .(LIGAND.CLASS, LIGAND.FAMILY, capwords(LIGAND, strict = T),
                   RECEPTOR.CLASS, RECEPTOR.FAMILY, capwords(RECEPTOR, strict = T), MECHANISM)]
setnames(lr.pairs, c('V3', 'V6'), c('LIGAND', 'RECEPTOR'))
l.non0 <- data.table(LIGAND = l.n,
                     l.non0 = l.nonzero)
r.non0 <- data.table(RECEPTOR = r.n,
                     r.non0 = r.nonzero)
lr.pairs <- merge(lr.pairs, l.non0, by = 'LIGAND')
lr.pairs <- merge(lr.pairs, r.non0, by = 'RECEPTOR')
lr.pairs[, lr.non0.s := l.non0 + r.non0]
lr.pairs[, lr.non0.p := l.non0 * r.non0]

## to start, pick a few LR pair from each of the following classes:
## Growth Factor
lr.pairs[LIGAND.CLASS == 'Growth Factor', ][order(lr.non0.s, decreasing = T), ][1:3, ]
lr.pairs[LIGAND.CLASS == 'Growth Factor', ][order(lr.non0.p, decreasing = T), ][1:3, ]
## top 3 by either metric:   MDK_ITGB1, IGF2_IGF2R, IGF2_IGF1R

## Matrix
lr.pairs[LIGAND.CLASS == 'Matrix', ][order(lr.non0.s, decreasing = T), ][1:3, ]
lr.pairs[LIGAND.CLASS == 'Matrix', ][order(lr.non0.p, decreasing = T), ][1:3, ]
## top 3 by either metric:  COL1A2_CD44, COL1A1_CD44, FN1_CD44

## Spatial Guidance
lr.pairs[LIGAND.CLASS == 'Spatial Guidance', ][order(lr.non0.s, decreasing = T), ][1:3, ]
lr.pairs[LIGAND.CLASS == 'Spatial Guidance', ][order(lr.non0.p, decreasing = T), ][1:3, ]
## top 3 by either metric: SLIT2_ROBO2, RTN4_TNFRSF19, SLIT2_ROBO1

## Cytokine
lr.pairs[LIGAND.CLASS == 'Cytokine', ][order(lr.non0.s, decreasing = T), ][1:3, ]
lr.pairs[LIGAND.CLASS == 'Cytokine', ][order(lr.non0.p, decreasing = T), ][1:3, ]
## top 3 by either metric: CXCL12_ITGB1

# from this, make a list of Ls and Rs to run
m.to.mod <- l.to.mod <- r.to.mod <- character(0)
top.n <- 5
for(cc in c('Growth Factor', 'Matrix', 'Spatial Guidance', 'Cytokine')){
  l.to.mod <- c(l.to.mod, lr.pairs[LIGAND.CLASS == cc, ][order(lr.non0.s, decreasing = T), ][1:top.n, LIGAND])
  r.to.mod <- c(r.to.mod, lr.pairs[LIGAND.CLASS == cc, ][order(lr.non0.s, decreasing = T), ][1:top.n, RECEPTOR])
  m.to.mod <- c(m.to.mod, lr.pairs[LIGAND.CLASS == cc, ][order(lr.non0.s, decreasing = T), ][1:top.n, MECHANISM])
}
l.to.mod <- (l.to.mod |> na.omit() |> unique())
r.to.mod <- (r.to.mod |> na.omit() |> unique())
m.to.mod <- (m.to.mod |> na.omit() |> unique())

print(glue('This will run {length(l.to.mod)+length(r.to.mod)} L or R models, resulting in {length(m.to.mod)} LR pairs'))

## also construct a data vector of the total counts/infoDepth at all places
## c.dat <- colSums( l.o@assays$RNA@layers$counts) + colSums( r.o@assays$RNA@layers$counts)


###########################
## mostly code from Sam: ##
###########################

## # Look at the spatial structure to determine sampling window
## for(gg in lr.n){
##   temp <- sub2[feat == gg,] # use this as a demo gene to see spatial patterns
##   temp.p <- ggplot(temp,
##                    aes(x = x,
##                        y= y,
##                        color = MIDCount))+
##     geom_point()+
##     theme_classic() +
##     ggtitle(gg) + ylim(c(9, 30000))
##   print(temp.p)
##   readline(prompt="Press [enter] to continue")
## }

## # Take a chunk out of the middle
## x.window <- c(23000,24000)
## y.window <- c(13000,14000)
## ggplot(sub,
##        aes(x = x,
##            y= y,
##            color = MIDCount))+
##   geom_point()+
##   theme_classic()+
##   xlim(x.window)+
##   ylim(y.window)

## # Define sys.small
## sys.small <- temp[temp$x >= x.window[1] &
##                     temp$x <= x.window[2] &
##                       temp$y >= y.window[1] &
##                       temp$y <= y.window[2],]
## # Define XY coordinates
## sys.small$coordinate <- paste(sys.small$x,sys.small$y,sep=".")
## length(unique(sys.small$coordinate)) #81264 points (ok size for AOZ method?)
## sub2 <- sys.small[sys.small$feat == 'Vegfa',] # use this as a demo gene to see spatial patterns
## ggplot(sub2,
##        aes(x = x,
##            y= y,
##            color = MIDCount))+
##   geom_point()+
##   theme_classic()

## # Make wide (compute intensive)
## gc()
## # Convert to data table
## temp.data.table <- data.table::data.table(sys.small)
## # Long to Wide
## temp.wide <- dcast(temp.data.table,feat ~ coordinate,value.var = 'MIDCount',fill = 0) # fills in missing values with 0 rather than NA
## gc()

## # Give rownames
## rownames(temp.wide) <- temp.wide$feat

## # Remove geneid from data frame
## temp.wide$feat <- NULL

## # Convert NA to 0
## # temp.wide[is.na(temp.wide)] <- 0 # compute intensive / inefficient
## # temp.wide.zeros <- temp.wide %>% gtools::na.replace(0) # compute intensive / inefficient
## # temp.wide.zeros <- temp.wide %>% replace(is.na(.), 0) # compute intensive / inefficient
## #temp.wide.zeros <- temp.wide %>% mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) # from https://stackoverflow.com/questions/18142117/how-to-replace-nan-value-with-zero-in-a-huge-data-frame


## # Make into a seurat object?
## test.seurat.assay <- CreateAssay5Object(temp.wide)
## test.seurat <- CreateSeuratObject(test.seurat.assay)
## test.seurat
## #Inspect
## range(test.seurat$nFeature_RNA)
## VlnPlot(test.seurat,'nFeature_RNA',pt.size=0)
## VlnPlot(test.seurat,'nCount_RNA',pt.size=0)
## # Add spatial metadata
## test.seurat$x <- sys.small$x
## test.seurat$y <- sys.small$y

## #Inspect
## rownames(test.seurat@assays$RNA@layers$counts) <- rownames(test.seurat) # a hack!
## FOI <- 'Vegfa'
## to.plot <- test.seurat@assays$RNA@layers$counts[FOI,]
## to.plot <- cbind(to.plot,test.seurat@meta.data)

## thresh.factor <- 1

## ggplot(to.plot,
##        aes(x=x,
##            y=y,
##            color = to.plot))+
##   geom_point(size = 0.2,aes(alpha = I(ifelse(to.plot < 1, 0, 1))))+
##   theme_minimal()+
##   ggtitle(FOI)+
##   scale_color_gradient2(low="blue",mid = 'white', high="red",
##                         limits = c(0.1, max(to.plot$to.plot)*thresh.factor),
##                         midpoint = max(to.plot$to.plot)*thresh.factor*0.5,
##                         oob = scales::squish)+
##   DarkTheme()

## # Ligands only
## ligands.only <- test.seurat[rownames(test.seurat) %in% lr.mouse$Ligand.ApprovedSymbol,]
## receptors.only <- test.seurat[rownames(test.seurat) %in% lr.mouse$Receptor.ApprovedSymbol,]

## # Save for transfer
## save(ligands.only,file = 'ligands.only.for.aoz.2024-10-27.Robj')
## save(receptors.only,file = 'receptors.only.for.aoz.2024-10-27.Robj')
## save(test.seurat,file = 'test.seurat.for.aoz.2024-10-27.Robj')


## ## interactive shiny used to find boundary
## library(plotly)
## library(shiny)
## library(htmlwidgets)

## initDF <- full.raw[feat == 'Cxcl12', .(as.numeric(x), as.numeric(y))]
## setnames(initDF, c("V1", "V2"), c("x", "y")) # 45047 rows, initially

## ui <- fluidPage(
##   plotlyOutput("myPlot", height = '1000px'),
##   #verbatimTextOutput("click")
##   )

## server <- function(input, output, session) {

##   js <- "
##     function(el, x, inputName){
##       var id = el.getAttribute('id');
##       var gd = document.getElementById(id);
##       var d3 = Plotly.d3;
##       Plotly.update(id).then(attach);
##         function attach() {
##           gd.addEventListener('click', function(evt) {
##             var xaxis = gd._fullLayout.xaxis;
##             var yaxis = gd._fullLayout.yaxis;
##             var bb = evt.target.getBoundingClientRect();
##             var x = xaxis.p2d(evt.clientX - bb.left);
##             var y = yaxis.p2d(evt.clientY - bb.top);
##             var coordinates = [x, y];
##             Shiny.setInputValue(inputName, coordinates);
##           });
##         };
##   }
##   "

##   clickposition_history <- reactiveVal(initDF)

##   observeEvent(input$clickposition, {
##     clickposition_history(rbind(clickposition_history(), input$clickposition))
##   })

##   output$myPlot <- renderPlotly({
##     plot_ly(initDF[y > 9500, ], x = ~x, y = ~y, type = "scatter", mode = "markers") %>%
##       onRender(js, data = "clickposition")
##   })

##   myPlotProxy <- plotlyProxy("myPlot", session)

##   observe({
##     plotlyProxyInvoke(myPlotProxy, "restyle", list(x = list(clickposition_history()$x), y = list(clickposition_history()$y)))
##   })

##   output$click <- renderPrint({
##     clickposition_history()
##   })
## }

## shinyApp(ui, server)
