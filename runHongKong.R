library(readxl)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
setwd("/Users/akdes/Library/CloudStorage/Box-Box/baylor//arya_akash_rnaseq_classifier/")
load("/Users/akdes/Library/CloudStorage/Box-Box/baylor/akash/expression/data/09192022_GE_RAW_Normalized_N330PrimarySamples.rda")

setwd("/Users/akdes/Library/CloudStorage/Box-Box/baylor//arya_akash_rnaseq_classifier/")
load("/Users/akdes/Library/CloudStorage/Box-Box/baylor/meningioma_raleigh_hongkong/Inhouse_Hongkong_RNA_Seq_Combined_v2.rda")
data_GS <- data
colnames(data_GS) <- colnames(data_GS)
ensembl <- nth(tstrsplit(rownames(data_GS), split = "\\."), n = 1)
ids <- data.frame(SYMBOL = mapIds(org.Hs.eg.db, keys = ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
ids$ensembl <- rownames(ids)
data_GS$SYMBOL <- ids$SYMBOL
data_GS <- data_GS[match(na.omit(unique(data_GS$SYMBOL)), data_GS$SYMBOL), ]
rownames(data_GS) <- data_GS$SYMBOL
data_GS <- data_GS[, !colnames(data_GS) %in% "SYMBOL"]
colnames(data_GS) <- colnames(data)


normalized_data_GS <- normalized_data
rownames(samplesheet)[305:608] <- make.unique(samples[match(rownames(samplesheet)[305:608], samples$ExpressionID),"Sample.Name"])
colnames(normalized_data_GS)[305:608] <-  make.unique(samples[match(colnames(normalized_data_GS)[305:608], samples$ExpressionID),"Sample.Name"])
colnames(data_GS)[305:608] <-  make.unique(samples[match(colnames(data_GS)[305:608], samples$ExpressionID),"Sample.Name"])


samples <- samplesheet
Consistent_Classifications <- read_excel("/Users/akdes/Library/CloudStorage/Box-Box/baylor//arya_akash_rnaseq_classifier/data/Consistent_Classifications.xlsx")
assign("Consistent_Classifications", Consistent_Classifications, envir = globalenv())
samples$Train <- rep("no", length(samples$Batch))
samples$Train[rownames(samplesheet) %in% Consistent_Classifications$`Expression ID`] <- "yes"

samples$CNV <- rep("no", length(samples$Batch))
samples$CNV[samples$Train == "yes" & samples$class == "A"] <- "yes"
rownames(samples) <- samples$ExpressionID

load("/Users/akdes/Library/CloudStorage/Box-Box/baylor/scell_meningioma_ucsf_bcm_sean/Rdata/091622_SCELL_UCSF_BCM_seuratObj_ann.rda")
new.cluster.ids <- c(
    "other", # 0
    "tumor", # 1
    "microglia", # 2
    "tumor", # 3
    "T-cell", # 4
    "other", # 5
    "EC", # 6
    "mesenchymal", # 7
    "other", # 8
    "microglia", # 9
    "other", # 10
    "other", # 11
    "other", # 12
    "T-cell", # 13
    "EC", # 14
    "other", # 15
    "other", # 16
    "other", # 17
    "other", # 18
    "other", # 19
    "other", # 20
    "other", # 21
    "other", # 22
    "other", # 23
    "other"
) # 24

names(new.cluster.ids) <- levels(seuratObj)
seuratObj <- RenameIdents(seuratObj, new.cluster.ids)
sub <- subset(seuratObj, cells = colnames(seuratObj)[seuratObj$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3", "MSC5", "MSC1", "MSC4")])
sub$class <- rep("", length(sub$seurat_clusters))
sub$class[sub$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3")] <- "C"
sub$class[sub$orig.ident %in% c("MSC5")] <- "B"
sub$class[sub$orig.ident %in% c("MSC1", "MSC4")] <- "A"
sub <- subset(sub, cells = colnames(sub)[Idents(sub) %in% c("tumor", "microglia", "T-cell", "EC", "mesenchymal")])
sub$cellType <- Idents(sub)
rownames(samples) <- samples$ExpressionID

bulk_count_data <- data_GS
bulk_normalized_data <- normalized_data_GS
samples <- samples
seurat_obj <- sub

save(list=c("bulk_count_data", "bulk_normalized_data", "samples", "seurat_obj"), file="HongKong_Meningioma_Data.Rdata")

library(readxl)
library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(Seurat)
setwd("/Users/akdes/Library/CloudStorage/Box-Box/baylor//arya_akash_rnaseq_classifier/")
load("HongKong_Meningioma_Data.Rdata")
library(CLIPPR)
##### bulkdata is raw data with gene symbol row ids
##### bulkdata is normalizedd data

# seurat_obj seurat_obj$cellType should be cell type, seurat_obj$class should be class
# bulk_count_data matrix with gene symbols
# bulk_normalized_data matrix with gene symbols,  
# samples$Train  no or yes
# samples$CNV no or yes for CNV control 


samples$class <- factor(samples$class)
rownames(samples) <- colnames(bulk_count_data)
object <- CreateCLIPPRObject(
    bulkdata = bulk_count_data,
    bulk_normalized_data = bulk_normalized_data,
    bulkdata_ss = samples,
    seurat_obj = seurat_obj
)
samples$class <- factor(samples$class)
rownames(samples) <- colnames(bulk_count_data)
project_id <- "HongKong_Meningioma"
object <- runFeatureSelection(object)
object <- extractCNVSignal(object, project_id=project_id)
load(paste0(project_id, "_BULK_finalChrMat_thr_1.rda"))
finalChrMat <- t(finalChrMat)
del <- data.frame(apply(finalChrMat[, ], 1, function(x) paste(names(which(x<(0))), collapse=",")))
amp <- data.frame(apply(finalChrMat[,  ], 1, function(x) paste(names(which(x>(0))), collapse=",")))

all(rownames(samples)==rownames(amp))

all(rownames(samples)==rownames(del))

object <- runCNVClassifier(object, fts = c("chr1p","chr14q", "chr22q"))
object <- runClassifier(object)
object <- runMetaModel(object, modelNames = c("RF.Pred_ECMarkers", "RF.Pred_bulk", "RF.Pred_mesenchymalMarkers"), CNV = T)


oob_error_bulk_scell <- do.call(cbind, lapply(object@PRED, function(x) x$Model$confusion[, "class.error"]))
oob_error_cnv <- object@CNV_PRED$RF.Pred_CNV_chrmean$Model$confusion[, "class.error"]
oob_error_meta <- object@META_PRED$Model$confusion[, 4]


predictions_bulk_scell <- do.call(cbind, lapply(object@PRED, function(x) x$Predictions[, -1]))
predictions_cnv <- object@CNV_PRED$RF.Pred_CNV_chrmean$Predictions[, -1]
colnames(predictions_cnv) <- paste0("RF.Pred_CNV_", colnames(predictions_cnv))
predictions_meta <- object@META_PRED$Predictions[, -1]
colnames(predictions_meta) <- paste0("RF.Pred_METAMODEL_", colnames(predictions_meta))

predictions_bulk_scell_class <- do.call(cbind, lapply(object@PRED, function(x) x$Class))
predictions_cnv_class <- object@CNV_PRED$RF.Pred_CNV_chrmean$Class
predictions_meta_class <- object@META_PRED$Class

results <- data.frame(
    ID = rownames(samples), isTrain = samples$Train, chr_mean = t(object@chr_mean[c("chr1p", "chr14q","chr22q"), ]), amp=amp, del=del, class = samples$class, predictions_meta_class,
    predictions_bulk_scell_class, predictions_cnv_class, predictions_meta, predictions_bulk_scell, predictions_cnv
)
rownames(results) <- rownames(samples)

results_1st <- results

openxlsx::write.xlsx(results, file = "/Users/akdes/HongKong_CLIPPR_results.xlsx")

casperObj <- object@casperObj
cnv.scale=3
assignInNamespace(x = "draw_matrix", value = draw_matrix2, 
    ns = asNamespace("pheatmap"))
assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", 
    ns = asNamespace("pheatmap"))
breaks <- seq(-1, 1, by = 0.2)
color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))

idx <- cumsum(table(casperObj@annotation.filt$Chr)[as.character(1:22)])
xlabel <- rep("", length(rownames(casperObj@control.normalized[[3]])))
half <- round(table(casperObj@annotation.filt$Chr)[as.character(1:22)]/2)[-1]
xpos <- c(half[1], (idx[-22] + half))
xlabel[xpos] <- 1:22

data2 <- (casperObj@control.normalized[[cnv.scale]])


ann <- data.frame(class=results$predictions_meta_class)
rownames(ann) <- rownames(samples)
ann <- ann[order(ann[,1]),, drop=F]

p<- pheatmap(t(data2[ ,rownames(ann)]), cluster_cols = F, cluster_rows = F, gaps_col = idx, 
    color = color, breaks = breaks, labels_col = xlabel,annotation_row =ann,
    show_rownames = T, filename = "heatmap.pdf",fontsize_row = 6)


pdf("heatmap.pdf", width=10, height=50)
p
dev.off()

results[match(soi,rownames(results)), ]



samples$Train_1st <- samples$Train

samples$Train_2nd <- (apply((results)[,6:13], 1, function(x) length(unique(x))==1 )==TRUE)
samples$Train[which(samples$Train_2nd)] <- "no"
samples$Train[which(samples$Train_2nd)] <- "yes"
samples$CNV <- "no"
samples$CNV[samples$Train == "yes" & samples$class == "A"] <- "yes"
rownames(samples) <- samples$ExpressionID


object <- CreateCLIPPRObject(
    bulkdata = bulk_count_data,
    bulk_normalized_data = bulk_normalized_data,
    bulkdata_ss = samples,
    seurat_obj = seurat_obj
)

object <- runFeatureSelection(object)
object <- extractCNVSignal(object)
object <- runCNVClassifier(object, fts = c("chr1p", "chr22q"))
object <- runClassifier(object)
object <- runMetaModel(object, modelNames = c("RF.Pred_ECMarkers", "RF.Pred_bulk", "RF.Pred_mesenchymalMarkers"), CNV = T)


oob_error_bulk_scell <- do.call(cbind, lapply(object@PRED, function(x) x$Model$confusion[, "class.error"]))
oob_error_cnv <- object@CNV_PRED$RF.Pred_CNV_chrmean$Model$confusion[, "class.error"]
oob_error_meta <- object@META_PRED$Model$confusion[, 4]


predictions_bulk_scell <- do.call(cbind, lapply(object@PRED, function(x) x$Predictions[, -1]))
predictions_cnv <- object@CNV_PRED$RF.Pred_CNV_chrmean$Predictions[, -1]
colnames(predictions_cnv) <- paste0("RF.Pred_CNV_", colnames(predictions_cnv))
predictions_meta <- object@META_PRED$Predictions[, -1]
colnames(predictions_meta) <- paste0("RF.Pred_METAMODEL_", colnames(predictions_meta))

predictions_bulk_scell_class <- do.call(cbind, lapply(object@PRED, function(x) x$Class))
predictions_cnv_class <- object@CNV_PRED$RF.Pred_CNV_chrmean$Class
predictions_meta_class <- object@META_PRED$Class

results <- data.frame(
    ID = rownames(samples), isTrain = samples$Train, chr_mean = t(object@chr_mean[c("chr1p", "chr22q"), ]), class = samples$class, predictions_meta_class,
    predictions_bulk_scell_class, predictions_cnv_class, predictions_meta, predictions_bulk_scell, predictions_cnv
)
rownames(results) <- rownames(samples)

openxlsx::write.xlsx(results, file = "/Users/akdes/CLIPPR_results_2nditeration.xlsx")


soi <- c("15-01-017-R3",
 "16-01-049-R1" ,"17-01-013-R1" ,  "17-01-035R",   "18-01-013R"   , "Meningiomas_TL-20-468276_20-01-022",   "DCC33R"    ,   "DCC34R"     ,  "DCC37R")


results[match(soi,rownames(results)), ]
