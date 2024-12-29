#! /bin/bash

## Author: Limin Chen; Jun 10th, 2024
## Add images to both seurat rds and scanpy-h5ad ann data for SAW(stereoSeq) spatial data.

set -e
# Function to display the usage of the script
## set default values
hires=4
lowres=1
binsize=200
registif="" # registered image from saw pipeline.
ingem=""
library_id="S123"
outpath="/home/test"

display_usage(){
  echo "Usage: $0 [options]"
  echo "Options:"
  echo " -t, --registif <value>   Set regist.tif from SAW (default: ${registif})"  
  echo " -i, --ingem <value>     Set input h5ad from stereopy conversion (default: ${ingem})" 
  echo " -b, --binsize <value>    Set binsize for stereopy conversion (default: ${binsize})" 
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
    -i|--ingem) ingem="$2"; shift 2;;
    -b|--binsize) binsize="$2"; shift 2;;
    -H|--hires) hires="$2"; shift 2;;
    -l|--lowres) lowres="$2"; shift 2;;
    -d|--library_id) library_id="$2"; shift 2;;
    -o|--outpath) outpath="$2"; shift 2;;
    -h|--help) display_usage; exit 0;;
    *) echo "Error: Unknown option $key"; display_usage; exit 1;;
  esac
done


# Check if all arguments are provided
if [[ -z ${registif} || -z ${ingem} ]]; then
  echo "Error: Invalid argument: $key"
  display_usage
  exit 1
fi



# Your script logic here
echo "registif: ${registif}"
echo "ingem: ${ingem}"
echo "binsize: ${binsize}"
echo "hires: ${hires}"
echo "lowres: ${lowres}"
echo "library_id: ${library_id}"
echo "outpath: ${outpath}"


#########################
mkdir -p ${outpath}/addimage
chmod 777 ${outpath}/addimage


# step0: format saw8 gem to standard gem for data conversion
chipID=$(echo "$(basename $ingem)" | cut -d"." -f1)
gemname=$(echo "$(basename $ingem)")


nline=$(gunzip -c ${ingem} | head -n20 | awk -F'\t' 'NR==10{print NF}')
if [ "${nline}" -eq 6 ]; then
  echo "Saw8 gem file!"
  gunzip -c ${ingem} | head -n8  > ${outpath}/addimage/header.gem
  gunzip -c ${ingem} | awk 'BEGIN {OFS="\t"} NR >= 9 {for (i=2; i<=NF; i++) printf "%s%s", $i, (i<NF?OFS:ORS)}' | awk 'BEGIN {OFS="\t"} NR==1 {$1="geneID"} 1' > ${outpath}/addimage/counts2.gem
  cat ${outpath}/addimage/header.gem ${outpath}/addimage/counts2.gem | gzip > ${outpath}/addimage/${gemname}
fi
# step1: generate low and hi resolution images
convert -sample ${hires}%x${hires}% ${registif} ${outpath}/addimage/tissue_hires_image.png
convert -sample ${lowres}%x${lowres}% ${registif} ${outpath}/addimage/tissue_lowres_image.png

# step2: convert gem to rds and add lowres image to rds
hi=$(echo "scale=2; ${hires}/100" | bc)
lo=$(echo "scale=2; ${lowres}/100" | bc)
#Rscript /usr/local/bin/h5ad2rds.R -i ${inh5ad} -b ${binsize} --image_dir ${outpath}/addimage/ -l ${lo} --hires ${hi}
if [ "${nline}" -eq 5 ]; then
  file1=${ingem}
else 
  file1=${outpath}/addimage/${gemname}
fi

Rscript /usr/local/bin/gem2rds.r --gemf ${file1} --bin ${binsize} --image_dir ${outpath}/addimage/ --lowres ${lo} --hires ${hi}

echo "Image is added to Seurat RDS file!"
# step3-1 generate h5ad without image from rds, which is the input for scanpy adding images
Rscript /usr/local/bin/rds2h5adsc.R -r ${outpath}/addimage/${chipID}.addimg.rds

# step3-2: add low- and hi- resolution images to h5ad from step3-1, the result h5ad overwrites the input h5ad.
mkdir -p ${outpath}/addimage/spatial
chmod 777 ${outpath}/addimage/*
mv ${outpath}/addimage/tissue_*.png ${outpath}/addimage/spatial
mv ${outpath}/addimage/tissue_*.csv ${outpath}/addimage/spatial
mv ${outpath}/addimage/*.json ${outpath}/addimage/spatial
# add both hires and lowres image to h5ad 
python3 /usr/local/bin/rds2annImg.py --inputDir ${outpath}/addimage --library_id ${library_id}

# clean output
if [ "$nline" -eq 6 ]; then  
  rm ${outpath}/addimage/*.gem*
fi
