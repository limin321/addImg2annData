library(argparser)
library(Seurat)
library(SeuratDisk)

args <- arg_parser("Convert rds to scapy h5ad")
args <- add_argument(args, "--rds_path", help = "Path of rds with image")
args <- parse_args(args)

print("Convert rds-to-scanpy-h5ad") 
#rds_path=paste0(imageDir,"/S135TL_D1.tissue_seurat.rds")
rds_path=args$rds_path
rds2h5adscanpy <- function(rdsf){
  imageDir=dirname(rdsf)
  obj = readRDS(rdsf)
  SaveH5Seurat(obj, filename =paste0(imageDir,"/tissue_sc.h5Seurat"))
  SeuratDisk::Convert(paste0(imageDir,"/tissue_sc.h5Seurat"), dest = "h5ad") # this h5ad is ready for adding images in python
  file.remove(paste0(imageDir, "/tissue_sc.h5Seurat"))
}


rds2h5adscanpy(rdsf = rds_path)
print("Image added to seurat RDS")
#print("Image is added to your seurat and scanpy h5ad is generated! Well Done!!")
quit('yes', 0)
