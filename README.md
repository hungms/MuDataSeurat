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

## Usage
### 


### Export to H5AD, H5MU

MuDataSeurat only export 3 layers: `count`, `data`, `scale.data`

Therefore, You need `JoinLayers` for each modality first with seurat v5.

```R
library(tidyverse)
library(MuDataSeurat)
SeuratToH5AD(seuratobj, h5ad = "tests/seuratobj.h5ad", umap = "umap", pca = "pca", snn = "RNA_snn", nn = "RNA_nn", assay = "RNA", adt = NULL)
```

### Read H5AD to Seurat

```R
ReadH5AD()
ReadH5MU()
```
You may also use the native support of anndata in R: `anndataR::read_h5ad`
