#' convert_seurat_to_anndata
#'
#' Convert Seurat to AnnData with SeuratDisk::Convert()
#' @param obj Seurat object
#' @param assay assay name for gene expression. Defaults to "RNA"
#' @param adt_name assay name for adt_name counts if any
#' @param columns a vector of metadata columns to keep. Defaults to NULL to keep all columns
#' @param pca reduction name for X_pca. Defaults to "pca"
#' @param umap reduction name for X_umap. Defaults to "umap"
#' @param snn graph name for SNN. Defaults to "RNA_snn"
#' @param nn graph name for NN. Defaults to "RNA_nn"
#' @param h5ad path for output .h5ad file
#' @import stringr
#' @export
SeuratToH5AD <- function(obj, assay = "RNA", h5ad, umap, adt_name = NULL, columns = NULL, reductions = NULL, snn = NULL, nn = NULL){

    # umap must be present
    stopifnot(umap %in% names(obj@reductions))

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
    
    # add adt_name features to metadata
    if(!is.null(adt_name)){
        stopifnot(adt_name %in% names(obj@assays))
        adt_features <- paste0(tolower(adt_name), "_", rownames(obj[[adt_name]]))
        adt_counts <- FetchData(obj, vars = adt_features, slot = "counts")
        obj@meta.data[colnames(adt_counts)] <- NULL
        obj@meta.data <- cbind(obj@meta.data, adt_counts)}

    # convert assay to V3
    for(i in names(obj@assays)){
        if(!str_detect(i, "SCT|integrated")){
            obj[[i]] <- as(object = obj[[i]], Class = "Assay")}}

    # store PCA/UMAP/neighbours
    DefaultAssay(obj) <- assay

    obj@reductions[["umap"]] <- obj@reductions[[umap]]
    obj@reductions[["umap"]]@assay.used <- assay

    if(!is.null(reductions)){
        stopifnot(all(reductions %in% names(obj@reductions)))
        for(i in reductions){
            obj@reductions[[i]] <- obj@reductions[[i]]
            obj@reductions[[i]]@assay.used <- assay}
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