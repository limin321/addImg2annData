# addImg2annData
Seurat and Scanpy are two popular spatail data analysis tool, which mainly developed based on the 10X Visium platform. Spatial data from Stereo-Seq platform has different format, which can be read-in by either Seurat or Scanpy. Though Stereopy provides code to convert the Stereo-Seq data to be seurat and Scanpy compatible ann data, the important image layer is lost after conversion. To solve this problem, I created a docker image called addImg2annData. This tool allows you to add the iamge to both seurat RDS and annData for scanpy. The final annData struture is the same as those you read-in 10X Visium data.

## Installation
### Where to download? 
The docker image is available here <https://hub.docker.com/r/limin321/addimg2anndata>

**Docker**
```
docker pull limin321/addimg2anndata
```

**Singularity**
```
singularity pull addimg2anndata.sif docker://limin321/addimg2anndata:0.0.0
```

## Tutorials
*Convert Stereo-Seq data to Seurat RDS and Scanpy h5ad file*


