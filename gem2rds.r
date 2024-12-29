#! /usr/bin/Rscript

library(Matrix)
library(Seurat)
library(data.table)
library(R.utils)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(argparser)

args <- arg_parser("Converting gem.gz file to RDS.")
args <- add_argument(args, "--gemf", help = "input .gem.gz file from SAW pipeline")
args <- add_argument(args, "--bin", help = "bin size, default is 200", default=200)
args <- add_argument(args, "--image_dir", help = "The path of the tissue_lowres_image.png", default = ".")
args <- add_argument(args, "--lowres", help = "scale of image low resolution, ex, 0.01", default = 0.01)
args <- add_argument(args, "--hires", help = "scale of image high resolution, ex, 0.04", default = 0.04)
argv <- parse_args(args)

gemf=argv$gemf
assay='Spatial'
binsize = argv$bin
lowres=argv$lowres
hires=argv$hires
imageDir=argv$image_dir

prefix = strsplit(basename(gemf),"\\.")[[1]][1]
outfile = paste0(prefix, ".addimg.rds") # "B03209A212.addimg.rds"

## functions: convert Gem to Seurat obj with images
gem2seurat <- function(
    gemf= gemf, # take gem.gz
    binsize = binsize,
    prefix = prefix
){
  data <- fread(gemf)
  data$x <- as.numeric(data$x)
  data$y <- as.numeric(data$y)
  
  #### Step 2 convert to different binsize --- Group counts into bins
  data$x <- trunc(data$x / binsize) * binsize
  data$y <- trunc(data$y / binsize) * binsize
  
  # reorganize the count matrix to be gene by binsize(similar to cell)
  data$MIDCount <- as.numeric(data$MIDCount)
  data <- data[, .(counts=sum(MIDCount)), by = .(geneID, x, y)]
  
  #' create sparse matrix from stereo
  data$cell <- paste0(prefix, ':', data$x, '-', data$y)
  data$geneIdx <- match(data$geneID, unique(data$geneID))
  data$cellIdx <- match(data$cell, unique(data$cell))
  
  mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts,
                      dimnames = list(unique(data$geneID), unique(data$cell))) # 43795 x 15611
    
  #### Step 3 create meta data for seurat object.  
  cell_coords <- unique(data[,c("cell", "x", "y")]) # 15611 x 3
  rownames(cell_coords) <- cell_coords$cell
  
  #### Step 4 Create seurat obj
  seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'stereo', assay = 'Spatial',names.delim = ':', meta.data = cell_coords)
  
  #### Step 5, add meta.features to counts slot
  spatial_data <- GetAssayData(seurat_spatialObj, assay = "Spatial", slot = "counts")
  n_counts <- rowSums(spatial_data)
  n_cells <- rowSums(spatial_data > 0)
  mean_umi <- round(n_counts/n_cells)
  
  meta_features <- data.frame(
      n_cells = n_cells,
      n_counts = n_counts,
      mean_umi = mean_umi
  )
  seurat_spatialObj[["Spatial"]]@meta.features <- meta_features
  
  #### add mitochondria percentage
  mtPattern <- '^mt-|^Mt-|^MT-'
  seurat_spatialObj <- PercentageFeatureSet(seurat_spatialObj, mtPattern, col.name = "percent.mito")
  
  #saveRDS(seurat_spatialObj, file = paste0(imageDir,"/",prefix,".temp.rds"))
  return(seurat_spatialObj)

}

# Add images to seurat obj created by gem2seurat.
addimg2rds <- function(obj=object,
                       binsize=binsize,
                       lowres=lowres,
                       hires=hires,
                       imageDir=imageDir,
                       outfile=outfile
){
  ## Generate tissue_positions_list.csv
  tissue_positions_list <- data.frame(row.names = rownames(obj@meta.data),
                                      tissue = 1,
                                      row = obj$y,
                                      col = obj$x,
                                      imagerow = obj$y,
                                      imagecol = obj$x)
  write.table(tissue_positions_list, paste0(imageDir,"/tissue_positions_list.csv")  , sep = ',', quote = FALSE, col.names = FALSE)

  ## Generate json
  spot_diameter_fullres=sqrt(((0.22*as.numeric(binsize))^2)*2)
  scalefactors_json <- toJSON(list(spot_diameter_fullres = spot_diameter_fullres,
                                   fiducial_diameter_fullres = 600, # take 10x Value bc stereo has no fiducial
                                   tissue_hires_scalef =hires,
                                   tissue_lowres_scalef = lowres))

  write(scalefactors_json, file =paste0(imageDir,"/scalefactors_json.json"))

  image <- Read10X_Image(imageDir,
                         image.name = "tissue_lowres_image.png",
                         assay = "Spatial",
                         slice = "slice1",
                         filter.matrix = TRUE)
  # the image id sometimes got messed up in a format of scientific. The following code to make the image-ids are in the correct format. The image id is from tissue_positions_list rownames
  imgIDs <- Cells(image)
  idx <- c() # extract which ids has format issue
  for (e in imgIDs){
    if(grepl("e\\+", e)){
      x <- which(imgIDs %in% e)
      idx <- append(idx,x)
    }
  }
  for(i in idx){ # correct image-ids
    imgIDs[i] <- as.character(format(as.numeric(imgIDs[i]), scientific = FALSE))
  }

  rownames(image@coordinates) <- imgIDs
  image <- image[Cells(x = obj)]
  DefaultAssay(image) <- 'Spatial'
  obj[['slice1']] <- image
  saveRDS(obj, file =paste0(imageDir,"/", outfile))
  return(obj)
}

obj1 <- gem2seurat(gemf = gemf,
                  binsize = binsize,
                  prefix = prefix)

obj2 <- addimg2rds(obj = obj1,
           binsize = binsize,
           lowres=lowres,
           hires=hires,
           imageDir = imageDir,
           outfile=outfile
           )

print("Generate QC plot")
png(paste0(imageDir,"/spatialFeature_QC.png"))
SpatialFeaturePlot(obj2, features = "nFeature_Spatial", stroke = 0.25,pt.size.factor = 1, image.alpha = 1, alpha = 0.2)
dev.off()








