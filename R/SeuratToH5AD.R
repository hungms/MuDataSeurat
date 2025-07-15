#' convert_seurat_to_anndata
#'
#' Convert Seurat to AnnData with SeuratDisk::Convert()
#' @param obj Seurat object
#' @param h5ad path for output .h5ad file
#' @param columns a vector of metadata columns to keep. Defaults to NULL to keep all columns
#' @param pca reduction name for X_pca. Defaults to "pca"
#' @param umap reduction name for X_umap. Defaults to "umap"
#' @param assay assay name for gene expression. Defaults to "RNA"
#' @param adt assay name for ADT counts if any
#' @param return_genes if TRUE, return genes from assays named c("CC", "BCR", "TCR", "MHC") to current assay. Defaults to TRUE
#' @import stringr
#' @export
SeuratToH5AD <- function(obj, h5ad, assay = "RNA", adt = NULL, columns = NULL, pca = NULL, umap = NULL, snn = NULL, nn = NULL){

    # downgrade to v4
    options(Seurat.object.assay.version = "v3")
    if(!is.object(obj)){
        if(file.exists(obj)){
            obj <- qread(paste0(obj))}}
    
    # preselect columns if necessary
    if(!length(columns) > 0){
        message(paste0("selecting all columns from metadata"))
        columns <- colnames(obj@meta.data)}
    else{
	    message(paste0("selecting columns from metadata: ", paste0(columns, collapse = ", ")))
        obj@meta.data <- obj@meta.data[,columns]}
    
    # convert non-numeric columns to character class
    classlist <- c()
    for(i in seq_along(colnames(obj@meta.data))){
        if(is.integer(obj@meta.data[[i]]) | is.numeric(obj@meta.data[[i]])){
            obj@meta.data[[i]] <- ifelse(is.na(obj@meta.data[[i]]), 0, obj@meta.data[[i]])
            obj@meta.data[[i]] <- as.numeric(obj@meta.data[[i]])}
        if(!is.numeric(obj@meta.data[[i]])){
            obj@meta.data[[i]] <- ifelse(is.na(obj@meta.data[[i]]), "NA", obj@meta.data[[i]])
            obj@meta.data[[i]] <- as.character(obj@meta.data[[i]])}}
    
    # add adt features to metadata
    if(!is.null(adt)){
        stopifnot(adt %in% names(obj@assays))
        adt_features <- paste0(tolower(adt), "_", rownames(obj[[adt]]))
        adt <- FetchData(obj, vars = adt_features, slot = "counts")
        obj@meta.data[colnames(adt)] <- NULL
        obj@meta.data <- cbind(obj@meta.data, adt)}

    # convert assay to V3
    for(i in names(obj@assays)){
        if(!str_detect(i, "SCT|integrated")){
            obj[[i]] <- as(object = obj[[i]], Class = "Assay")}}

    # store PCA/UMAP/neighbours
    DefaultAssay(obj) <- assay

    if(!is.null(pca)){
        stopifnot(pca %in% names(obj@reductions))
        obj@reductions[["pca"]] <- obj@reductions[[pca]]
        obj@reductions[["pca"]]@assay.used <- assay
        }
    if(!is.null(umap)){
        stopifnot(umap %in% names(obj@reductions))
        obj@reductions[["umap"]] <- obj@reductions[[umap]]
        obj@reductions[["umap"]]@assay.used <- assay
        }
    if(!is.null(snn)){
        stopifnot(snn %in% names(obj@graphs))
        obj@graphs[[paste0(assay, "_snn")]] <- obj@graphs[[snn]]
        }
    if(!is.null(nn)){
        stopifnot(nn %in% names(obj@graphs))
        obj@graphs[[paste0(assay, "_nn")]] <- obj@graphs[[nn]]
        }

    # write to h5ad
    MuDataSeurat::WriteH5AD(obj, h5ad, assay=assay, scale.data = F) #https://github.com/zqfang/MuDataSeurat #remotes::install_github("zqfang/MuDataSeurat", dependencies = F)
    options(Seurat.object.assay.version = "v5")
    }