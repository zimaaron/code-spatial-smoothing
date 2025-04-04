# This function takes spatial single-cell object as input
# The user can choose the assay to pull these data from
# We will use the FANTOM5 database as a ground-truth, for now, but can be updated or customized later
# We only care about queried mechanisms that are actually in the rownames of our input dataset

# This function yields a 3-part named list as output, which includes
# a. ligand information (a matrix of normalized data, with annotated rows (genes) and columns (spots)
# b. receptor information (a matrix of normalized data, with annotation rows (genes) and columns (spots)
# c. location information (x and y coordinates for each spot)

# The ligand and receptor outputs are segrated so that they can be handled differently if desired
# We are not preserving index-matching of ligands and receptors, for now, since we can re-create that later

LRExtract <- function(object,
                      assay = "Spatial",
                      LR.database = "fantom5",
                      species = 'mouse'){
  require(NICHES)
  
  # Load up LR.database
  if(LR.database=='fantom5'){
    if(species=='mouse'){
      ground.truth <- ncomms8866_mouse
    }
    if(species=='rat'){
      ground.truth <- ncomms8866_rat
    }
    if(species=='pig'){
      ground.truth <- ncomms8866_pig
    }
    if(species=='human'){
      ground.truth <- ncomms8866_human
    }
  }
  
  # Prune to those mechanisms actually IN the dataset
  ligands <- ground.truth[ground.truth$Ligand.ApprovedSymbol %in% rownames(object),]$Ligand.ApprovedSymbol
  receptors <- ground.truth[ground.truth$Receptor.ApprovedSymbol %in% rownames(object),]$Receptor.ApprovedSymbol
  
  # Pull ligand info from assay of interest (normalized data slot) - as a full matrix
  ligand.info <- as.matrix(object@assays[[assay]]@data[ligands,])
  
  # Pull receptor info from assay of interest (normalized data slot) - as a full matrix
  receptor.info <- as.matrix(object@assays[[assay]]@data[receptors,])
  
  # Extract XY info
  location.info <- object@meta.data[,c('x','y')]
  
  # Return output
  output <- list(ligand.info,receptor.info,location.info)
  names(output) <- c('ligand.info','receptor.info','location.info')
  return(output)
}