
#' @details
#' The main functions you will need to use are CreateCLIPPRObject() and runCLIPPR(clippr_object).
#' For additional details on running the analysis step by step, please refer to the example vignette.
#' @aliases clippr-package
"_PACKAGE"


#' The clippr Class
#'
#'
#' The clippr Class
#' The clippr object is .....
#' @name clippr
#' @rdname clippr
#' @aliases clippr-class
#' @exportClass clippr

clippr <- methods::setClass(
  "clippr", 
  slots = c(
    bulkdata = "ANY", 
    bulk_normalized_data="ANY", 
    bulkdata_ss = "ANY", 
    seurat_obj="ANY",
    casperObj="ANY",
    BULKFTS="ANY", 
    SCELLFTS="ANY", 
    BULKFTS_RF="ANY", 
    BULKFTS_ADJ="ANY",
    PRED="ANY",
    CNV_PRED="ANY", 
    META_PRED="ANY",
    chr_mean="ANY",
    chr_median="ANY")
)
