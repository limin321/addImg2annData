# addImg2annData
Seurat and Scanpy are two popular spatail data analysis tool, which mainly developed based on the 10X Visium platform. Spatial data from Stereo-Seq platform has different format, which can be read-in by either Seurat or Scanpy. Though Stereopy provides code to convert the Stereo-Seq data to be seurat and Scanpy compatible ann data, the important image layer is lost after conversion. To solve this problem, I created a docker image called addImg2annData. This tool allows you to add the iamge to Stereo-seq data and generate both seurat RDS and annData for scanpy. This tool is compatible with both SAW7 and SAW8 gem.gz file. The final annData struture is the same as those you read-in 10X Visium data.

## Installation
### Where to download? 
The docker image is available [here](https://hub.docker.com/r/limin321/addimg2anndata)

**Docker**
```
docker pull limin321/addimg2anndata:0.0.1
```

**Singularity**
```
export SINGULARITY_CACHEDIR=/home/your/dir
export TMPDIR=/home/your/dir
singularity pull addimg2anndata.sif docker://limin321/addimg2anndata:0.0.1
```

## Tutorials
### Convert Stereo-Seq data to Seurat RDS and Scanpy h5ad file
#### Prepare inputs
Required inputs:
1) regist.tif. This file is the output of SAW pipeline, normally located in the `03.register/ssDNA_SS200000135TL_D1_regist.tif` for Saw7; and `outs/image/SS200000135TL_D1_ssDNA_regist.tif` for Saw8.
2) gem.gz file. For example, xxx.tissue.gem.gz file. The gem file can be converted from gef file. Please refer SAW pipeline on how to generate gem file if not provided in the standard SAW outputs.
3) bin_size. Default is bin=200. This is from previous step, but won't really affect this tool running.
4) hires. Default is 4, which means scaling the original tif image to 4%.
5) lowres. Default is 1, which means scaling the original tif image to 1%.
6) library_id. Default is "S123". Set your read library ID value.

##### Example of running
**Docker**
Use default value, assuming you have the input tif and h5ad in the local_input_dir
```
docker run --rm -v <local_input_dir>:/home/test limin321/addimg2anndata:0.0.1 bash addImage.sh -t /home/test/ssDNA_SS200000135TL_D1_regist.tif -i /home/test/S135TL_D1.tissue.gem.gz
```
Set your own value, assuming you have the input tif and gem.gz in the local_input_dir
```
docker run --rm -v /stomics_data/formatConvertion:/home/test limin321/addimg2anndata:0.0.1 bash addImage.sh -t /home/test/seurat/ssDNA_SS200000135TL_D1_regist.tif -i /home/test/seurat/S135TL_D1.tissue.gem.gz -H 6 -l 2 -d SS200000135TL -o /home/test/seurat
```
###### Convert cellbin gem.gz file.
```
 docker run --rm -v /stomics_data/liminData/project/cellbinconvert/cellbin:/home/test limin321/addimg2anndata:0.0.1 bash addImage.sh -t /home/test/SS200000135TL_DAPI_regist.tif -i /home/test/SS200000135TL.adjusted.cellbin.gem.gz -b 1 -c TRUE
 ```

###### Convert lasso data -- use bin1 to do lasso
To add image to lasso data, 1) get the json file by doing lasso in StereoMap4, making sure it is bin1; 2) run `saw reanalyze lasso ...` to get a folder called segmentation, which include the gem.gz and mask.tif files; for cellbin data, you need to run `saw convert gef2gem ...` to get the gem.gz file. Eventually, you need three files as inputs: xx.lasso.xx.gem.gz, xx.lasso.mask.tif, regist.tif.

Here are the example code:
```
# cellbin lasso data
docker run --rm -v /Users/liminchen/Limchen/apps/docker/addimgtest:/home/test limin321/addimg2anndata:0.0.1 addImage.sh -t /home/test/SS200000135TL_ssDNA_regist.tif -i /home/test/SS200000135TL.subcell.label.cellbin.gem.gz -m /home/test/SS200000135TL.lasso.cellbin.subcell.mask.tiff -b 1 -c TRUE

# square bin1 lasso data
docker run --rm -v /Users/liminchen/Limchen/apps/docker/addimgtest:/home/test limin321/addimg2anndata:0.0.1 bash addImage.sh -t /home/test/SS200000135TL_ssDNA_regist.tif -i /home/test/SS200000135TL.lasso.bin1.1084_bin1.gem.gz -m /home/test/SS200000135TL.lasso.1084_bin1.mask.tif
docke
```


**Singularity**

Here I have `ssDNA_SS200000135TL_D1_regist.tif` and `S135TL_D1.tissue.gem.gz` inputs files in the `/stomics_data/formatConvertion/seurat` folder.
```
# use default values
singularity run -B /stomics_data/formatConvertion:/home/test addimg2anndata.sif bash addImage.sh -t ./seurat/ssDNA_SS200000135TL_D1_regist.tif -i ./seurat/S135TL_D1.tissue.gem.gz

# specify your values
singularity run -B /stomics_data/formatConvertion:/home/test addimg2anndata.sif bash addImage.sh -t ./seurat/ssDNA_SS200000135TL_D1_regist.tif -i ./seurat/S135TL_D1.tissue.gem.gz -H 6 -l 2 -b 100 -d SS200000135TL -o ./seurat

# replace the <place_holder> with your data path
singularity run -B </stomics_data/formatConvertion>:/home/test addimg2anndata.sif bash addImage.sh -t <./seurat/ssDNA_SS200000135TL_D1_regist.tif> -i <./seurat/S135TL_D1.tissue.gem.gz> -H 6 -l 2 -b 100 -d SS200000135TL -o <./seurat>

```
###### convert cellbin gem.gz file.
make sure to set arguments `-b 1 -c TRUE`, see following example.
```
singularity run -B /stomics_data/liminData/project/cellbinconvert:/home/test addimg2anndata.sif bash addImage.sh -t ./cellbin/SS200000135TL_DAPI_regist.tif -i ./cellbin/SS200000135TL.adjusted.cellbin.gem.gz -H 6 -l 2 -b 1 -d SS200000135TL -c TRUE -o ./cellbin/
```

### Outputs
For successful run, you expect to have following outputs in the output path you specified (Default is the workding dir).
```
(base) [len@localhost formatConvertion]$ tree seurat/addimage/
seurat/addimage/
├── S135TL_D1.addimg.rds
├── spatial
│   ├── scalefactors_json.json
│   ├── tissue_hires_image.png
│   ├── tissue_lowres_image.png
│   └── tissue_positions_list.csv
├── spatialFeature_QC.png
└── tissue_sc.h5ad

2 directories, 7 files

```
You should see the folder `addimage`, containing 3 files and 1 folder, named spatial. <br />
The 3 files are: <br />
> *S135TL_D1.addimg.rds* -- the seurat RDS data with image. It can be loaded to R by `readRDS()` function. <br />
> *tissue_sc.h5ad* -- the scanpy h5ad data with images. It can be read in scanpy by `sc.read_h5ad()` function. <br />
> *spatialFeature_QC.png* -- this is QC plot to show the image is added to RDS successfully. The region withot tissue should be black, the image should be the background. <br />
> *spatial* -- the folder contains images related files, which mimic files from Visium Spaceranger output `spatial` folder. <br />



### Credits and Acknowledgements
Great thanks and give credits to the [Seurat](https://github.com/satijalab/seurat) team and [Scanpy](https://github.com/scverse/scanpy) team. This tool borrows codes from their achievements. 


