#! /bin/bash

## Author: Limin Chen; Jun 10th, 2024
## Add images to both seurat rds and scanpy-h5ad ann data for SAW(stereoSeq) spatial data.

set -e
# Function to display the usage of the script
## set default values
hires=4
lowres=1
registif="" # registered image from saw pipeline.
inh5ad=""
library_id="S123"
outpath=$(pwd)

display_usage(){
  echo "Usage: $0 [options]"
  echo "Options:"
  echo " -t, --registif <value>   Set regist.tif from SAW (default: ${registif})"  
  echo " -i, --inh5ad <value>     Set input h5ad from stereopy conversion (default: ${inh5ad})" 
  echo " -H, --hires <value>      (Optional) Set high resolution value (default: 4, interpret as 4/100=0.04))" 
  echo " -l, --lowres <value>     (Optional) Set low resolution value (default: 1, interpret as 1/100=0.01))" 
  echo " -d, --library_id <value>     Set library_id (default: $library_id)" 
  echo " -o, --outpath <paths>      (Optional) Set output dir path (default: ${outpath})" 
  echo " -h, --help               Display this help message"
  echo "Example:"
  echo " $0 -t ssDNA_SS200000135TL_D1_regist.tif -i S135TL_D1.tissue_seurat.h5ad -h 4 -l 1 -d "s12345" -o ."
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case "$1" in
    -t|--registif) registif="$2"; shift 2;;
    -i|--inh5ad) inh5ad="$2"; shift 2;;
    -H|--hires) hires="$2"; shift 2;;
    -l|--lowres) lowres="$2"; shift 2;;
    -d|--library_id) library_id="$2"; shift 2;;
    -o|--outpath) outpath="$2"; shift 2;;
    -h|--help) display_usage; exit 0;;
    *) echo "Error: Unknown option $key"; display_usage; exit 1;;
  esac
done


# Check if all arguments are provided
if [[ -z ${registif} || -z ${inh5ad} ]]; then
  echo "Error: Invalid argument: $key"
  display_usage
  exit 1
fi



# Your script logic here
echo "registif: ${registif}"
echo "inh5ad: ${inh5ad}"
echo "hires: ${hires}"
echo "lowres: ${lowres}"
echo "library_id: ${library_id}"
echo "outpath: ${outpath}"


#########################

# apt-get update
# apt install imagemagick -y
# rm /etc/ImageMagick-6/policy.xml
# cp policy.xml /etc/ImageMagick-6/policy.xml

# apt-get install bc
# apt-get install python3.9 -y
# pip install pandas
# pip install matplotlib
# pip install scanpy==1.9.6
# pip install --upgrade numba numpy



mkdir -p ${outpath}/addimage
chmod 777 ${outpath}/addimage

# step1: generate low and hi resolution images
convert -sample ${hires}%x${hires}% ${registif} ${outpath}/addimage/tissue_hires_image.png
convert -sample ${lowres}%x${lowres}% ${registif} ${outpath}/addimage/tissue_lowres_image.png

# step2: add lowres image to rds
hi=$(echo "scale=2; ${hires}/100" | bc)
lo=$(echo "scale=2; ${lowres}/100" | bc)
Rscript h5ad2rds.R -i ${inh5ad} -b 200 --image_dir ${outpath}/addimage/ -l ${lo} --hires ${hi}

# step2-1 generate h5ad without image from rds, which is the input for scanpy adding images
out=$(basename "${inh5ad}")
outRDS="${out%.h5ad}"
Rscript rds2h5adsc.R -r ${outpath}/addimage/${outRDS}.RDS

# step3: add low- and hi- resolution images to h5ad from step2-1, the result h5ad overwrites the input h5ad.
mkdir -p ${outpath}/addimage/spatial
chmod 777 ${outpath}/addimage/*
mv ${outpath}/addimage/tissue_*.png ${outpath}/addimage/spatial
mv ${outpath}/addimage/tissue_*.csv ${outpath}/addimage/spatial
mv ${outpath}/addimage/*.json ${outpath}/addimage/spatial
# add both hires and lowres image to h5ad 
python3 rds2annImg.py --inputDir ${outpath}/addimage --library_id ${library_id}

# bash addImage.sh -t ../ssDNA_SS200000135TL_D1_regist.tif -i ../seurat/S135TL_D1.tissue_seurat.h5ad
# registif: ../ssDNA_SS200000135TL_D1_regist.tif
# inh5ad: ../seurat/S135TL_D1.tissue_seurat.h5ad
# hires: 4
# lowres: 1
# library_id: S123
# outpath: /home/test/stereopy_demo/formatConvertion/scripts

# bash addImage.sh -t ../ssDNA_SS200000135TL_D1_regist.tif -i ../seurat/S135TL_D1.tissue_seurat.h5ad -H 6 -l 2 -d ss321 -o /home/test/stereopy_demo/formatConvertion
# registif: ../ssDNA_SS200000135TL_D1_regist.tif
# inh5ad: ../seurat/S135TL_D1.tissue_seurat.h5ad
# hires: 6
# lowres: 2
# library_id: ss321
# outpath: /home/test/stereopy_demo/formatConvertion