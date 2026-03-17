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
SeuratToH5AD <- function(obj, assay = "RNA", h5ad, umap, adt_name = NULL,
                          columns = NULL, reductions = NULL, snn = NULL, nn = NULL) {

    stopifnot(umap %in% names(obj@reductions))
    stopifnot(assay %in% names(obj@assays))

    options(Seurat.object.assay.version = "v3")
    on.exit(options(Seurat.object.assay.version = "v5"), add = TRUE)

    # Subset and coerce metadata columns
    if (!is.null(columns)) {
        message("selecting columns from metadata: ", paste(columns, collapse = ", "))
        obj@meta.data <- obj@meta.data[, columns, drop = FALSE]
    } else {
        message("selecting all columns from metadata")
    }
    for (i in seq_along(colnames(obj@meta.data))) {
        col <- obj@meta.data[[i]]
        if (is.integer(col) || is.numeric(col)) {
            obj@meta.data[[i]] <- as.numeric(ifelse(is.na(col), 0, col))
        } else {
            obj@meta.data[[i]] <- as.character(ifelse(is.na(col), "NA", col))
        }
    }

    # Append ADT counts to metadata
    if (!is.null(adt_name)) {
        stopifnot(adt_name %in% names(obj@assays))
        adt_features <- paste0(tolower(adt_name), "_", rownames(obj[[adt_name]]))
        adt_counts   <- FetchData(obj, vars = adt_features, slot = "counts")
        obj@meta.data[colnames(adt_counts)] <- NULL
        obj@meta.data <- cbind(obj@meta.data, adt_counts)
    }

    # Downgrade all v5 assays to v3 Assay
    for (i in names(obj@assays)) {
        if (!str_detect(i, "SCT|integrated")) {
            message("converting ", i, " to Assay")
            assay_obj <- obj[[i]]
            if (inherits(assay_obj, "Assay5")) {
                obj[[i]] <- CreateAssayObject(counts = LayerData(assay_obj, layer = "counts"))
            } else {
                obj[[i]] <- as(assay_obj, "Assay")
            }
        }
    }

    # Set the active umap reduction and optionally update others
    obj@reductions[["umap"]] <- obj@reductions[[umap]]
    obj@reductions[["umap"]]@assay.used <- assay
    if (!is.null(reductions)) {
        stopifnot(all(reductions %in% names(obj@reductions)))
        for (i in reductions) {
            obj@reductions[[i]]@assay.used <- assay
        }
    }

    # Register graph slots
    if (!is.null(snn)) {
        stopifnot(snn %in% names(obj@graphs))
        obj@graphs[[paste0(assay, "_snn")]] <- obj@graphs[[snn]]
    }
    if (!is.null(nn)) {
        stopifnot(nn %in% names(obj@graphs))
        obj@graphs[[paste0(assay, "_nn")]] <- obj@graphs[[nn]]
    }

    MuDataSeurat::WriteH5AD(obj, h5ad, assay = assay, scale.data = FALSE)
}