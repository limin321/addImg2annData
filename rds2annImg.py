#! /usr/bin/env python
from __future__ import annotations
import json
import argparse
from pathlib import Path
from typing import TYPE_CHECKING
import pandas as pd
from matplotlib.image import imread
import scanpy as sc


def main():
    parser = argparse.ArgumentParser(description="Add images to annData")
    parser.add_argument('--inputDir', type=str, required=True, help="Dir path where the h5ad file and /spatial folder locate")
    parser.add_argument('--library_id', type=str, required=True, help="The library ID of the experiment.")
    args = parser.parse_args()

    # input dir path, where it contains: S135TLD1.h5ad  spatial(scalefactors_json.json  tissue_hires_image.png  tissue_lowres_image.png  tissue_positions_list.csv)
    inputDir = args.inputDir
    library_id = args.library_id

    h5ad_path = list(Path(inputDir).glob('*.h5ad'))[0]
    adata=sc.read_h5ad(h5ad_path) # converted from rds by seurat

    # remove all meta data from R to make the exact same structure as reading from visium data
    adata.obs = adata.obs[[]]
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][library_id] = dict()

    # prepare image-related files
    path1 = Path(inputDir + "/spatial/")
    tissue_positions_file = (path1 / "tissue_positions_list.csv")
    files = dict(
        tissue_positions_file=tissue_positions_file,
        scalefactors_json_file=path1 / "scalefactors_json.json",
        hires_image=path1 / "tissue_hires_image.png",
        lowres_image=path1 /"tissue_lowres_image.png"
    )

    # add hires and lowres image
    adata.uns["spatial"][library_id]["images"] = dict()
    for res in ["hires", "lowres"]:
        try:
            adata.uns["spatial"][library_id]["images"][res] = imread(
                str(files[f"{res}_image"]) # this search the file name pattern.
            )
        except Exception:
            raise OSError(f"Could not find '{res}_image'")

    # read json scalefactors
    adata.uns["spatial"][library_id]["scalefactors"] = json.loads(
        files["scalefactors_json_file"].read_bytes()
    )

    attrs = dict()
    adata.uns["spatial"][library_id]["metadata"] = {
        k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
        for k in ("chemistry_description", "software_version")
        if k in attrs
    }

    # read coordinates
    positions = pd.read_csv(
        files["tissue_positions_file"],
        header=0 if tissue_positions_file.name == "tissue_positions.csv" else None,
        index_col=0
    )

    positions.columns = [
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]

    # merge positions to adata
    adata.obs.index = adata.obs.index.astype(str)
    positions.index = positions.index.astype(str)
    adata.obs = adata.obs.join(positions, how="left")
    adata.obsm["spatial"] = adata.obs[
    ["pxl_row_in_fullres","pxl_col_in_fullres"]
    ].to_numpy()

    adata.obs.drop(
        columns=["pxl_row_in_fullres","pxl_col_in_fullres"],
        inplace=True
    )

    # put image path in uns
    source_image_path = inputDir + "/spatial/tissue_hires_image.png"
    if source_image_path is not None:
        # get an absolute path
        source_image_path = str(Path(source_image_path).resolve())
        adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
            source_image_path
        )

    adata.write(inputDir + '/tissue_sc.h5ad')
    print("Finished!!")

if __name__ == "__main__":
    main()

## Example:
#  python rds2annImg.py --inputDir /stomics_data/liminData/stereopy_demo/formatConvertion/scanpy/input --library_id SS20000135TLD1
