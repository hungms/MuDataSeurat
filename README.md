# MuDataSeurat

**This is a fork from zqfang's [forked version](https://github.com/zqfang/MuDataSeurat) of [MuDataSeurat](https://github.com/PMBio/MuDataSeurat)**

Please refer to the original repo for more details

## Why this fork ?

This fork enables streamline conversion of a multimodal Seurat object into AnnData/MuData objects which can be uploaded onto the [CELLXGENE_VIP]() server for the Calado lab. This fork corrects for the class assignments of Seurat V5's metadata, enabling integer class to be viewed as a numeric vector on the server.

## Installation

Please install the main branch

```R
remotes::install_github("hungms/MuDataSeurat", dependencies = F)
```

## Example Usage

### Convert Seurat to H5AD
`MuDataSeurat::SeuratToH5AD()` allows you to convert Seurat object with multiple reductions and/or have ADT expression data.  

An input is required for the `umap` argument for successful conversion.  

`umap`, `reductions`, `snn`, `nn`, arguments allow you specify which reductions and nearest-neighbour graphs to keep.  

```R
MuDataSeurat::SeuratToH5AD(seuratobj, assay = "RNA", umap = "umap", h5ad = "tests/seuratobj.h5ad", adt_name = "ADT", columns = NULL, reductions = c("pca", "phate"), snn = "RNA_snn", nn = "RNA_nn")
```

### Read H5AD to Seurat
```R
ReadH5AD()
ReadH5MU()
```
You may also use the native support of anndata in R: `anndataR::read_h5ad`
