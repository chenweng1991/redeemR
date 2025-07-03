#' convert_mitotracing_redeemR
#'
#' This function is to convert the older version mitotracing into the redeemR
#' @param mitotracing mitotracing object
#' @return redeemR class
#' @export
convert_mitotracing_redeemR<-function(mitotracing,addAssignedVariant=T, V.fitered=F){
ob<-new("redeemR")
ob@GTsummary.filtered<-mitotracing@GTsummary.filtered
ob@CellMeta<-mitotracing@CellMeta
if (V.fitered){
    ob@V.fitered<-mitotracing@V.fitered.list[[1]]
}
ob@HomoVariants<-""
ob@UniqueV<-mitotracing@UniqueV
ob@Cts.Mtx<-mitotracing@Cts.Mtx
ob@Cts.Mtx.bi<-mitotracing@Cts.Mtx.bi
ob@TREE<-mitotracing@TREE
ob@DistObjects<-mitotracing@DistObjects
ob@para<-mitotracing@para
ob@Seurat<-mitotracing@Seurat
ob@Ctx.Mtx.depth<-mitotracing@Ctx.Mtx.depth
if (addAssignedVariant){
ob@AssignedVariant <- mitotracing@AssignedVariant
}    
return(ob)    
}

## Internal function to convert the variant names, implemented by convert_variant
convert_variant_1 <- function(input_string) {
  if (grepl("^Variants[0-9]+[A-Za-z]{2}$", input_string)) {
    # Convert from "Variants10000GA" to "10000_G_A"
    number <- gsub("Variants([0-9]+)[A-Za-z]{2}$", "\\1", input_string)
    letters <- gsub("^Variants[0-9]+([A-Za-z]{2})$", "\\1", input_string)
    output_string <- paste0(number, "_", substr(letters, 1, 1), "_", substr(letters, 2, 2))
  } else if (grepl("^[0-9]+_[A-Za-z]_[A-Za-z]$", input_string)) {
    # Convert from "10000_G_A" to "Variants10000GA"
    parts <- strsplit(input_string, "_")[[1]]
    output_string <- paste0("Variants", parts[1], parts[2], parts[3])
  } else {
    stop("Input string format is not recognized.")
  }
  return(output_string)
}

#' convert_variant
#'
#' @param x this is a vector of variant names, either way for example from '10000_G_A' to/from 'Variants10000GA'
#' @return a vector of strings 
#' @export
convert_variant <- function(x){
    res<-sapply(x,convert_variant_1)
    return(as.character(res))
}


#' Convert a redeemR variant & depth matrix into long format for downstream analysis
#'
#' This function filters a variant count matrix (`mat_var`) to include only cells in the provided
#' whitelist and variants observed in at least two cells, applies the same filtering to the corresponding
#' depth matrix (`mat_depth`), and then melts both into a single long-format `data.table`.  The output
#' contains one row per cell–variant pair with its count (`a`) and depth (`d`).
#' 
#' Included in redeemR-2.0  2025-07-01
#' @param sample Character; an identifier for this dataset (printed to the console).
#' @param mat_var A numeric matrix of variant UMI counts (cells × variants).
#' @param mat_depth A numeric matrix of total depth (same dimensions and row/column names as `mat_var`).
#' @param cell_whitelist A character vector of cell barcodes to retain (rows of `mat_var`/`mat_depth`).
#'
#' @return A `data.table` with columns:
#'   - `cell`: cell barcode  
#'   - `variant`: variant identifier  
#'   - `a`: UMI count (filtered)  
#'   - `d`: corresponding depth  
#'
#' @details  
#' 1. Retains only rows in `mat_var` whose rownames appear in `cell_whitelist`.  
#' 2. Drops variants observed in fewer than two cells.  
#' 3. Ensures both matrices share the same dimensions before melting.  
#' 4. Returns a combined long table for joint count/depth analysis.  
#'
#' @import data.table
#' @importFrom glue glue
#' @export
convert_redeem_matrix_long <- function(sample="",mat_var, mat_depth, cell_whitelist){
    library(data.table)
    print(glue("{sample}=========="))
    # Start filtering
    print(glue("{length(row.names(mat_var))} total cells in this redeem dataset"))
    print(glue("{sum(row.names(mat_var) %in% cell_whitelist)} are HSCs in the HSC whitelist"))
    mat_var_filtered<-mat_var[row.names(mat_var) %in% cell_whitelist,]
    mat_var_filtered <- mat_var_filtered[,colSums(mat_var_filtered)>=2]
    mat_var_filtered <- mat_var_filtered[rowSums(mat_var_filtered)>0,]
    # Keep the same rows and columns for depth mat
    mat_depth_filtered<-mat_depth[row.names(mat_var_filtered),colnames(mat_var_filtered)]
    print(glue("The minimum variant has at least {min(colSums(mat_var_filtered))} cells sharing it;
    the minimum ecll has at least {min(rowSums(mat_var_filtered))} variants"))
    print(glue("Finally, after filtering, {nrow(mat_var_filtered)} cells, {ncol(mat_var_filtered)} variants"))    

    # Convert sparse matrix to dense matrix - var and depth
    dense_mat_var <- as.matrix(mat_var_filtered)
    dense_mat_depth <- as.matrix(mat_depth_filtered)

    # Convert to data.table (long format) -- var
    dt_var <- as.data.table(as.table(dense_mat_var))
    setnames(dt_var, c("cell", "variant", "a"))
    # Convert to data.table (long format) -- depth
    dt_depth <- as.data.table(as.table(dense_mat_depth))
    setnames(dt_depth, c("cell", "variant", "d"))
    # Combine variant counts and depth
    print(glue("The variant has {nrow(dt_var)}; The depth has {nrow(dt_depth)}. They should be the same.\n"))

    dt_var_depth<-cbind(dt_var,d=dt_depth$d)
    return(dt_var_depth)
}