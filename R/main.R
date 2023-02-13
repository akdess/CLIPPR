
CreateCLIPPRObject <- function(bulkdata, bulk_normalized_data,
                               bulkdata_ss,
                               seurat_obj) {
  object <- new(Class = "clippr", seurat_obj = seurat_obj, bulkdata_ss = bulkdata_ss, bulkdata = bulkdata, bulk_normalized_data = bulk_normalized_data)

  return(object)
}

runFeatureSelection <- function(object) {
  seurat_obj <- object@seurat_obj
  bulkdata_ss <- object@bulkdata_ss
  bulkdata <- object@bulkdata

  object <- runBULKFS(object)
  ct <- unique(seurat_obj$cellType)
  scMRKS <- list()
  for (i in 1:length(ct)) {
    sub <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj$cellType %in% ct[i]])
    message(paste0("Finding single cell markers for ", ct[i], " single cell data."))
    scMRKS[[paste0(ct[i], "Markers")]] <- runSCELLFS(sub, pervar = 75, numtree = 500, mrkr = TRUE, numfts = 500, mpct = .25, LFCT = 1)
  }
  object@SCELLFTS <- scMRKS
  object
}

extractCNVSignal <- function(objec, project_id) {
  data <- object@bulkdata
  bulkdata_ss <- object@bulkdata_ss

  annotation <- generateAnnotation(id_type = "hgnc_symbol", genes = rownames(data), ishg19 = F, centromere)
  data <- data[match(annotation$Gene, rownames(data)), ]
  control.sample.ids <- rownames(bulkdata_ss)[which(bulkdata_ss$CNV == "yes")]

  casperObj <- CreateCasperObject(
    raw.data = data, loh.name.mapping = NULL, sequencing.type = "bulk",
    cnv.scale = 3, loh.scale = 3, expr.cutoff = 1, filter = "mean", matrix.type = "raw",
    annotation = annotation, method = "iterative", loh = NULL,
    control.sample.ids = as.character(control.sample.ids), cytoband = cytoband
  )

  smoothed <- na.omit((casperObj@control.normalized[[3]]))
  chr_mean <- apply(smoothed[, ], 2, function(x) unlist(lapply(split(x, casperObj@annotation.filt$cytoband), mean)))
  chr_median <- apply(smoothed[, ], 2, function(x) unlist(lapply(split(x, casperObj@annotation.filt$cytoband), mean)))
  rownames(chr_mean) <- paste0("chr", rownames(chr_mean))
  rownames(chr_median) <- paste0("chr", rownames(chr_median))

  object@chr_mean <- chr_mean
  object@chr_median <- chr_mean

  casperObj <- runCaSpERWithoutLOH(casperObj, project=paste0(project_id, "_BULK")) 
 # del <- data.frame(apply(finalChrMat[, ], 1, function(x) paste(names(which(x<(0))), collapse=",")))
 # amp <- data.frame(apply(finalChrMat[,  ], 1, function(x) paste(names(which(x>(0))), collapse=",")))
  object@casperObj <- casperObj

  object
}

runCNVClassifier <- function(object, fts) {
  PRED <- list()
  PRED[[paste0("RF.Pred_CNV_chrmean")]] <- RFClass(object, data = object@chr_mean, fts = fts, ntree = 5000)
  PRED[[paste0("RF.Pred_CNV_chrmedian")]] <- RFClass(object, data = object@chr_mean, fts = fts, ntree = 5000)
  object@CNV_PRED <- PRED
  object
}

runClassifier <- function(object) {
  allFTS <- c(list(
    "bulk" = object@BULKFTS$symbol
  ), lapply(object@SCELLFTS, function(x) x$gene))
  # Classification Random Forest
  PRED <- list()

  for (i in 1:length(allFTS)) {
    PRED[[paste0("RF.Pred_", names(allFTS)[i])]] <- RFClass(object, data = object@bulk_normalized_data, fts = allFTS[[i]], ntree = 5000)
  }

  object@PRED <- PRED
  object
}

runMetaModel <- function(object, modelNames = c("RF.Pred_ECMarkers", "RF.Pred_bulk", "RF.Pred_mesenchymalMarkers"), CNV = T) {
  pred <- do.call(cbind, lapply(object@PRED[modelNames], function(x) x$Predictions[, -1]))

  if (CNV) {
    predCNV <- do.call(cbind, lapply(object@CNV_PRED["RF.Pred_CNV_chrmean"], function(x) x$Predictions[, -1]))
    pred <- cbind(pred, predCNV)
  }
  pred <- data.frame(pred)
  rownames(pred) <- colnames(object@bulk_normalized_data)


  META_PRED <- RFClass(object, data = t(pred), fts = colnames(pred), ntree = 5000)
  object@META_PRED <- META_PRED
  object
}

# DESeq feature selection
runBULKFS <- function(object, L2FC = 2, SignExcl = FALSE) {
  bulkdata_ss <- object@bulkdata_ss
  bulkdata <- object@bulkdata

  bulkdata <- bulkdata[, bulkdata_ss$Train == "yes"]
  bulkdata_ss <- bulkdata_ss[bulkdata_ss$Train == "yes", ]

  if (!is.null(bulkdata_ss$Batch)) {
    design <- "~Batch+class"
  } else {
    (
      design <- "~class"
    )
  }

  rds <- DESeq(DESeqDataSetFromMatrix(countData = (bulkdata), colData = bulkdata_ss, design = as.formula(design)))

  preDEGs <- list()
  x <- levels(bulkdata_ss$class)
  # computes all contrasts between values in factors of interest
  for (i in 1:length(x)) {
    for (ii in i:length(x)) {
      if (i == ii) {
        next
      }
      preDEGs[[paste0("res_", x[i], "vs", x[ii])]] <- results(rds, contrast = c("class", x[i], x[ii]), pAdjustMethod = "fdr")
    }
  }

  # Clean up results. Apply filters for p-adj., L2FC, and map ENSID to gene symbol.
  for (i in 1:length(preDEGs)) {
    # 3a Convert to data table for easier proccesing
    preDEGs[[i]] <- setDT(as.data.frame(preDEGs[[i]]), keep.rownames = TRUE)[]
    # 3b Trim ENSID for mapping to gene symbol
    #  preDEGs[[i]]$rn <- substr(preDEGs[[i]]$rn,0,15)
    # 3c Map to gene symbol
    preDEGs[[i]]$symbol <- preDEGs[[i]]$rn
    # Filters
    # Filters for significant DEGs
    preDEGs[[i]] <- preDEGs[[i]][padj < .001]
    # Filters for L2FC if a threshold is specified.
    if (is.numeric(L2FC)) {
      preDEGs[[i]] <- preDEGs[[i]][abs(preDEGs[[i]]$log2FoldChange) > L2FC]
    }
  }

  # Collects intersections of contrasts for a given FOI value.
  DEGs <- list()
  for (i in 1:length(x)) {
    # Store all DEGs that came up for a factor across all contrasts
    for (ii in i:length(x)) {
      if (i == ii) {
        next
      }
      # Stores DEGs for factor value of loop i
      DEGs[[x[i]]] <- rbind(DEGs[[x[i]]], preDEGs[[paste0("res_", x[i], "vs", x[ii])]])
      DEGs[[x[ii]]] <- rbind(DEGs[[x[ii]]], preDEGs[[paste0("res_", x[i], "vs", x[ii])]])
    }
  }

  # Filter for DEGs that show up in all of an FOI value's contrasts
  for (i in 1:length(DEGs)) {
    # Frequency counts of each DEG
    FREQ <- as.data.frame(table(DEGs[[x[i]]]$symbol))
    # Create list of DEGs that came up in all contrasts. Frequency =
    # number of contrasts gene came up in. Levels-1 = number of total contrasts.
    FREQ <- FREQ[FREQ$Freq == (length(x) - 1), 1]
    # Filter for DEGs that only showed up in some contrasts.
    DEGs[[x[i]]] <- DEGs[[x[i]]][DEGs[[x[i]]]$symbol %in% FREQ, ]
    # Add FOI value information to genes that pass this filter.
    DEGs[[x[i]]]$class <- x[i]
  }

  # Removes genes that came up in contrasts that did not include FOI value.
  nonx <- list()
  # Collects all DEGs from contrasts that exclude FOI value for current loop
  for (i in 1:length(DEGs)) {
    # Collect all DEGs in contrasts excluding the FOI value of loop i.
    # Iterates through all contrasts
    for (ii in 1:length(DEGs)) {
      # Skips contrasts that include FOI value of loop i
      if (i == ii) {
        next
      }
      for (zz in 1:length(DEGs)) {
        # Skips contrasts that include FOI value of loop i
        if (i == zz) {
          next
        }
        # Collects DEGs for contrast which excludes FOI value of loop i.
        nonx[[x[i]]] <- rbind(nonx[[x[i]]], preDEGs[[paste0("res_", x[ii], "vs", x[zz])]])
      }
    }

    # Filter out DEGs that showed up in all contrasts of FOI value of loop i, that also showed up in contrasts that excluded FOI value of loop i.
    DEGs[[x[i]]] <- DEGs[[x[i]]][which(!(DEGs[[x[i]]]$symbol %in% nonx[[x[i]]]$symbol))]
  }
  # Aggregate DEGs
  FTS <- data.frame()
  for (i in 1:length(DEGs)) {
    FTS <- rbind(FTS, cbind(DEGs[[i]]))
  }
  # Check for any genes that show up as distinct in more than one group.
  # Find all unique combinations of FOI and gene.
  # Genes that passed the filter when they shouldn't will show up as duplicates in symbol column.
  UNI <- unique(FTS[, c("class", "symbol")])
  # Remove DEGs from aggregated features that are duplicated in the symbol column of unique.
  FTS <- FTS[which(!(FTS$symbol %in% UNI$symbol[which(!isUnique(UNI$symbol))]))]
  # Now FTS contains DESeq data for each gene from each contrast for the FOI value for which it is unique.I.E. genes are entered more than once.

  # Remove excess entries.
  FTS <- distinct(FTS, symbol, .keep_all = TRUE)

  object@BULKFTS <- FTS
  object
}


# #Functions to accompany modularized scripts
# #DESeq feature selection
# blkDEG_ADJ <- function(object)
# {
#   bulkdata_ss <- object@bulkdata_ss
#   bulkdata <- object@bulkdata

#   bulkdata_ss$class <- as.factor(bulkdata_ss$class)
#   #Checks if batch effect has been specified
#     if(!is.null(bulkdata_ss$Batch)){
#       bulkdata_ss$Batch <- as.factor(bulkdata_ss$Batch)
#     }
#   #1) Create variables to separate out cohorts by factor of interest (FOI) in DESeq analysis
#     #Get different values for factor of interest
#       x <- levels(bulkdata_ss$class)
#     #For each value in factor of interest, create new column indicating binary status of sample for given value.
#     for (i in 1:length(x)) {
#       #Create xADJ column which indicates Y/N status for value in factor of interst.
#       samp[,paste0(x[i],"ADJ")] <- as.factor(ifelse(samp$FOI == x[i], 'Y','N'))
#     }

#   #2) Run DESeq Analysis
#     coldata <- samp
#     #initialize list to store results
#     DESResList <- list()
#     for (i in 1:length(x)) {
#       #Establish design of DESeq analysis
#         des <- paste0(x[i],"ADJ")
#       #Include batch effect if one is specified.
#         if(!is.null(Batch)){
#           des <-paste0("~Batch + ",des)
#         }
#       #Perform DESeq analysis and store result in list.
#       DESResList[[paste0("res_",x[i])]] <- results(DESeq(DESeqDataSetFromMatrix(countData = bulkdata, colData = bulkdata_ss, design= as.formula(des))))
#     }

#   #3) Clean up results
#     for (i in 1:length(x)){
#       #3a Convert to data table for easier proccesing
#         DESResList[[paste0("res_",x[i])]] <-setDT(as.data.frame(DESResList[[paste0("res_",x[i])]]), keep.rownames = TRUE)[]
#       #3b Trim ENSID for mapping to gene symbol
#         DESResList[[paste0("res_",x[i])]]$rn <- substr(DESResList[[paste0("res_",x[i])]]$rn,0,15)
#       #3c Map to gene symbol
#         DESResList[[paste0("res_",x[i])]]$symbol <- mapIds(org.Hs.eg.db,keys=DESResList[[paste0("res_",x[i])]]$rn,column="SYMBOL",keytype="ENSEMBL",multiVals="first")
#     }

#     #4) Select features from DEG results
#     for (i in 1:length(x)){
#       DESResList[[paste0(x[i],"_fts.pre")]] <- dplyr::filter(DESResList[[paste0("res_",x[i])]],
#                                                              #Filter for padj <.05
#                                                               DESResList[[paste0("res_",x[i])]]$padj<.05 &
#                                                              #Filter for abs(L2FC) > 2.5
#                                                               abs(DESResList[[paste0("res_",x[i])]]$log2FoldChange)>L2FC)
#     }

#     #3e Filters out features that are found in more than one group
#       if(SignExcl){
#         #Creates two datasets for each factor value to store positive and negative DEGs
#         for (i in 1:length(x)){
#           DESResList[[paste0("pre",x[i],"_ftspos")]] <- setDT(dplyr::filter(DESResList[[paste0(x[i],"_fts.pre")]],DESResList[[paste0(x[i],"_fts.pre")]]$log2FoldChange >0))
#           DESResList[[paste0("pre",x[i],"_ftsneg")]] <- setDT(dplyr::filter(DESResList[[paste0(x[i],"_fts.pre")]],DESResList[[paste0(x[i],"_fts.pre")]]$log2FoldChange <0))
#         }
#         #Filters out DEGs that have the same sign in multiple multiple factor values
#         for (i in 1:length(x)){
#           #create vectors to store DEGs of pos/neg sign in other factor values
#           otherPos <- c()
#           otherNeg <- c()
#           for (ii in 1:length(x)){
#             if(i == ii){next}
#             #stores genes of pos/neg sign of other factor values
#             otherPos <- append(otherPos,DESResList[[paste0("pre",x[ii],"ftspos")]]$symbol)
#             otherNeg <- append(otherneg,DESResList[[paste0("pre",x[ii],"ftsneg")]]$symbol)
#           }
#           #Filter out DEGs for factor value of interest that match the sign of DEGs in other factor values.
#           DESResList[[paste0(x[i],"_ftspos")]] <- setDT(dplyr::filter(DESResList[[paste0("pre",x[i],"ftspos")]], DESResList[[paste0("pre",x[i],"ftspos")]]$symbol %in% otherPos==FALSE))
#           DESResList[[paste0(x[i],"_ftsneg")]] <- setDT(dplyr::filter(DESResList[[paste0("pre",x[i],"ftspos")]], DESResList[[paste0("pre",x[i],"ftspos")]]$symbol %in% otherNeg==FALSE))
#         }
#         message("Returing DEGs that are type exclusive, considering L2FC direction.")
#         #Aggregate genes
#         FTS <- c()
#         for (i in 1:length(x)) {
#           FTs <- append(FTS, DESResList[[paste0(x[i],"_ftspos")]])
#           FTs <- append(FTS, DESResList[[paste0(x[i],"_ftsneg")]])
#         }
#         FTS<- unique(FTS)
#       }else{
#         for (i in 1:length(x)){
#           #Store significant DEGs from values outside of value of interst in vector
#           #Create empty vector for storage
#           otherDEGs <- c()
#           #Iterate through DEGs for all values in factor
#           for (ii in 1:length(x)){
#             #skips storing factor value DEGs when on factor value of interest for given loop.
#             if(x[i] == x[ii]){next}
#             #Stores DEGs when on factor value that is not the factor value of interst.
#             otherDEGs <- append(otherDEGs,DESResList[[paste0(x[ii],"_fts.pre")]]$symbol)
#           }
#           #Filter DEGs for factor value of interst
#           DESResList[[paste0(x[i],"_fts")]] <- setDT(dplyr::filter(DESResList[[paste0(x[i],"_fts.pre")]],
#                                                                    #For DEGs that are not DEGs in other factor values.
#                                                                    DESResList[[paste0(x[i],"_fts.pre")]]$symbol %in% c(otherDEGs)==FALSE))
#         }
#         message("Returning DEGs that are type exclusive,ignoring L2FC direction.")
#         #Aggregate genes
#         FTS <- c()
#         for (i in 1:length(x)){
#           FTS <- append(FTS,  DESResList[[paste0(x[i],"_fts")]]$symbol)
#         }
#         FTS <-unique(FTS)
#       }
#       object@BULKFTS_ADJ <- FTS
#       object
# }


# #Bulk Random Forest Feature Selection
# runBULKFS_RF  <- function(object, numvar = 10000, numtree = 5000, numfts = 5000){

#   bulkdata_ss <- object@bulkdata_ss
#   bulkdata <- object@bulkdata

#   normdatsd <- apply(bulkdata,1,sd)
#   #returns index positions of sds ordered from greatest to least
#     ord <- order(normdatsd,decreasing = TRUE)
#   #Sorts normalized sequencing data by sd
#     normdatsd <- bulkdata[ord,]
#   #takes sequencing data for top 10,0000 most variable genes.
#     normdatsd <- normdatsd[1:numvar,]
#     row.names(normdatsd) <- paste0("V", 1:nrow(normdatsd)) #appends index to each gene

#   #5b.2)Produce feature selection RF####
#   #Create random forest forfeature selection
#     message("Building Random Forest for feature selection.")
#     rf.varsel <- randomForest(as.factor(as.character(bulkdata_ss$class))~.,data = t(bulkdata),ntree=numtree,importance=TRUE)
#   #Store feature importance as a seperate variable
#     imp.meandecrease <- rf.varsel$importance[,dim(rf.varsel$importance)[2]-1] #stores permutation variable importance

#   #5b.3)Extracts top "n" most important features#####
#   if(numfts > length(imp.meandecrease)){
#     numfts <- length(imp.meandecrease)
#   }
#   or <- order(imp.meandecrease) #creates vector of indices of features ordered from most to least important.
#   datnamed <- bulkdata[ord,]
#   FTS <- datnamed[or[1:numfts],] #Creates important fts sequencing dataset from sequencing for entire cohort.
#   FTS <- rownames(FTS)
#   object@BULKFTS_RF <- FTS
# }

# Single Cell Feature Selection
runSCELLFS <- function(seurat_obj, pervar = 75, numtree = 500, numfts = 50, mrkr = FALSE, mpct = .9, LFCT = 1.5) {
  if (mrkr) {
    message("Selecting features from type specific markers in single cell data.")
    Idents(seurat_obj) <- "class"
    # finds markers for each group in factor of interest.
    MRKRS <- FindAllMarkers(seurat_obj, only.pos = F, min.pct = mpct, logfc.threshold = LFCT)
    MRKRS <- MRKRS[unique(MRKRS$gene), ]
    if (nrow(MRKRS) < numfts) {
      message(paste0("Requested ", numfts, " features.There were only ", nrow(MRKRS), " found. Returning all of them."))
      FTS <- MRKRS # get all info about markers, in addition to symbol.
      # FTS <- MRKRS$gene
    } else {
      # Filter training data for biomarkers
      message("Initiating Random Forest Feature Selection on identified markers.")
      TrainDat <- as.matrix(GetAssayData(seurat_obj, slot = "scale.data"))
      TrainDat <- TrainDat[which(rownames(TrainDat) %in% rownames(MRKRS)), ]
      # Collect metadata for training data
      TrainDatmet <- seurat_obj@meta.data
      # Produce feature selection tree
      # Calculate standard deviations of each gene
      Traindatsd <- apply(TrainDat, 1, sd)
      # returns index positions of sds ordered from greatest to least
      ord <- order(Traindatsd, decreasing = TRUE)
      # Sorts normalized sequencing data by sd
      Traindatsd <- TrainDat[ord, ]
      # takes sequencing data for top 75% most variable genes.
      Traindatsd <- Traindatsd[1:(nrow(Traindatsd) * (pervar / 100)), ]
      # Create random forest for feature selection
      TrainDatmet$class <- as.factor(TrainDatmet$class)
      rownames(Traindatsd) <- janitor::make_clean_names(rownames(Traindatsd))
      # Create feature selection tree
      rf.varsel_tum <- randomForest(TrainDatmet$class ~ ., data = t(Traindatsd), ntree = numtree, importance = TRUE)
      # Extract top numfts most important features
      # stores permutation variable importance
      imp.meandecrease <- rf.varsel_tum$importance[, dim(rf.varsel_tum$importance)[2] - 1]
      # Gets order of features from most to least important
      or <- order(imp.meandecrease)
      # Pulls number of desired features from list.
      RFSFscMRKRS <- as.data.frame(TrainDat[or[1:numfts], ])
      FTS <- rownames(RFSFscMRKRS)
    }
  } else {
    message("Selecting features from all single cell data")
    # Produce feature selection tree
    # Establish training data
    TrainDat <- as.matrix(GetAssayData(seurat_obj, slot = "scale.data"))
    # Collect metadata for training data
    TrainDatmet <- seurat_obj@meta.data
    # Calculate standard deviations of each gene
    Traindatsd <- apply(TrainDat, 1, sd)
    # returns index positions of sds ordered from greatest to least
    ord <- order(Traindatsd, decreasing = TRUE)
    # Sorts normalized sequencing data by sd
    Traindatsd <- TrainDat[ord, ]
    # takes sequencing data for top pervar% most variable genes.
    Traindatsd <- Traindatsd[1:(nrow(Traindatsd) * (pervar / 100)), ]
    # Create random forest for feature selection
    TrainDatmet$class <- as.factor(TrainDatmet$class)
    rownames(Traindatsd) <- janitor::make_clean_names(rownames(Traindatsd))
    # Create feature selection tree
    rf.varsel <- randomForest(TrainDatmet$class ~ ., data = t(Traindatsd), ntree = numtree, importance = TRUE)
    # Extract top n most important features
    # stores permutation variable importance
    imp.meandecrease <- rf.varsel$importance[, dim(rf.varsel$importance)[2] - 1]
    # creates vector of indices of features ordered from most to least important.
    or <- order(imp.meandecrease)
    scRFSF <- as.data.frame(TrainDat[or[1:numfts], ])
    # Creates important fts sequencing dataset from sequencing for entire cohort.
    scRFSF$rn <- rownames(scRFSF)
    # Aggregate Random Forest selected Features
    FTS <- scRFSF$rn
  }
  return(FTS)
}


# Trains random forest classifier
RFClass <- function(object, data, fts, ntree = 5000) {
  # Create Training Data

  bulkdata_ss <- object@bulkdata_ss
  bulk_normalized_data <- data

  Train <- makeTrain(dat = bulk_normalized_data, ss = bulkdata_ss, fts)
  # Build RF on data for consistently classified samples
  RFmodel <- randomForest(as.factor(as.character(Train$samp$class)) ~ ., data = t(Train$datTrain), ntree = 5000, importance = TRUE)
  # Produce Predictions from model trained on consistently classified samples
  # Produce predictions
  RES <- data.table::as.data.table(predict(RFmodel, newdata = data.frame(t(Train$datTest)), type = "prob"), keep.rownames = TRUE)

  Class <- apply(as.matrix(RES[, -1]), 1, function(x) names(x)[which(x == max(x))])
  # Associate predictions with sample ID
  # Aggregate model and results for export
  RFResult <- list("Model" = RFmodel, "Predictions" = RES, "Class" = Class)
  RFResult
}

# #Trains Neural Network
# NNClass <- function(object, fts, n=20, maxI=10000, dec=.001, MNW = 2000){
# #Create training data

#   bulkdata_ss <- object@bulkdata_ss
#   bulkdata <- object@bulkdata
#   Train <- makeTrain(dat=bulkdata,ss=bulkdata_ss,fts)
# #Train neural network
#   NNmodel <- nnet(as.factor(as.character(Train$samp$class))~.,data = t(Train$datTrain), size=n, maxit=maxI, decay=dec, MaxNWts = MNW)
# #Generate predictions for neural network
#   RES <- data.table::as.data.table(predict(NNmodel, newdata=data.frame(t(Train$datTest)), type="class"))
# #Aggregate model and results for export
#   NMResult <- list("Model" = NNmodel,"Predictions" = RES)
# }

# #Trains SVM
# SVMClass <- function(object, fts, krn = "radial"){
# #Create Training Data
#  bulkdata_ss <- object@bulkdata_ss
#   bulkdata <- object@bulkdata
#   Train <- makeTrain(dat=bulkdata,ss=bulkdata_ss,fts)
# #
# #Build SVM on data for consistently classified samples.
#   SVMmodel <- svm(as.factor(as.character(Train$samp$class))~.,data = t(Train$datTrain), kernel = krn)
# #Produce Predictions from model trained on consistently classified samples
#   #Produce predictions
#     RES <- data.table::as.data.table(predict(SVMmodel, newdata=data.frame(t(Train$datTest)), type="prob"), keep.rownames = TRUE)

# #Aggregate model and results for export
#   SVMResult <- list("Model" = SVMmodel,"Predictions" = RES)
# }

# Creates Training Data that is normalized data filtered for data corresponding to selected features from consistently classified samples.
makeTrain <- function(dat, fts, ss) {
  # dat -> data from which training data is produced
  # fts -> features that are used to create training data

  # Filters normalized data for features
  fndata <- dat[which(rownames(dat) %in% fts), ]
  Datatest <- fndata
  rownames(Datatest) <- paste0("V", 1:nrow(Datatest))
  # Filters feature data for consistently classified samples
  # Consistently classified sample data
  samplesCC <- ss[ss$Train == "yes", ]
  # Data for consistently classified samples
  Datatrain <- Datatest[, rownames(samplesCC)]
  row.names(Datatrain) <- paste0("V", 1:nrow(Datatrain))
  # Returns sample info and seq. data for training data
  Train <- list(samp = samplesCC, datTrain = Datatrain, datTest = Datatest)
}


# Function for producing signature expression matrix from sc data.
SigExpMatMkr <- function(scobj, genes, cellTypes) {
  # Pull count data.
  expMat <- as.matrix(GetAssayData(scobj, slot = "counts"))
  # Pull single cell data corresponding to genes of interest.
  fexpMat <- expMat[genes, ]
  # Append cell type ID and gene name to data as data frame.
  x <- data.frame(rbind(c("gene", as.vector(Idents(scobj))), cbind(rownames(fexpMat), fexpMat)))
  x[, c(1, which(x[1, ] %in% cellTypes))]
}


# LOOCV <- function(object)
# {
# LOOCVres <- list()
#   for (i in names(Fts)){
#     LOOCVres[[i]] <- LOOCV(LOOCVsamps = samplesCC, FTS = Fts[[i]],xdat = as.data.frame(normalized_data,row.names = rownames(normalized_data)))
#   }
# }

# CombROC(LOOCVres)
#      MetModdat <- cbind(Mods[["Bulkmodel"]]$Predictions[,c("A","B","C")],
#                        Mods[["Mesenchymemodel"]]$Predictions[,c("A","B","C")],
#                        Mods[["Endothelialmodel"]]$Predictions[,c("A","B","C")]
#                       )
#     rownames(MetModdat) <- Mods[["Bulkmodel"]]$Predictions$rn
#     colnames(MetModdat) <- paste0("V",c(1:ncol(MetModdat)))
#   #Filter metamodel training data for consistently classified samples.
#     MetModdatTrain <- MetModdat[rownames(MetModdat) %in% samplesCC$ExpressionID,]
#   #Train metamodel
#     MetMod <- randomForest(as.factor(as.character(samplesCC$class))~.,data = MetModdatTrain,ntree=5000,importance=TRUE)
#   #Predict remaining samples with metamodel
#     MetModPred <- as.data.frame(predict(MetMod, newdata= data.frame(MetModdat), type = "prob"), row.names = rownames(MetModdat))
#     #Export predictions
#     # write_xlsx(cbind("rn" = samples$Sample.Name[which(samples$ExpressionID == rownames(MetModdat))],MetModPred),"C:/Users/arya       shetty/Desktop/RNA.Classifier_FinalMetaModelPredictions.xlsx")
#   #Leave one out cross validation
#       rownames(MetModdat) <- samples$ExpressionID
#       MetModdat <- MetModdat[which(rownames(MetModdat) %in% samplesCC$ExpressionID),]
#            n <- nrow(MetModdat)
#                  LOOCVResults <- data.frame(Sample = rep(NA,n),
#                                  Test = rep(NA,n),
#                                  Pred = rep(NA,n),
#                                  Check = rep(NA,n))
#       for (i in 1:n){
#         #Prep model inputs
#          modDat <- MetModdat[-i,]
#          modClass <- as.factor(as.character(samplesCC$class[-i]))
#         #Train Random Forest
#             model <- randomForest(modClass~.,data = modDat, ntree=5000,
#                                   importance=TRUE)
#         #Predict left out sample
#           LOOCVResults$Sample[i]<- samplesCC$ExpressionID[i]
#           LOOCVResults$Test[i] <- levels(modClass)[samplesCC$class[i]]
#           LOOCVResults$Pred[i] <- levels(modClass)[predict(model,
#                                                           newdata=MetModdat[i,],
#                                                           type="class")]
#           LOOCVResults[i,levels(as.factor(samplesCC$class))] <- predict(model,
#                                                                      newdata=MetModdat[i,],
#                                                                      type="prob")
#         #Tally/record result
#           LOOCVResults$Check[i] <- ifelse(LOOCVResults$Test[i] == LOOCVResults$Pred[i], "Y","N")
#       }


#     LOOCVres[["MetaModel"]] <- LOOCVResults
#     CombROC(LOOCVres)

# library(tidytext)
# #Iterate through LOOCVres
# AUCs <- CalcAUC(LOOCVres)
# AUCs$Model[which(!(AUCs$Model %in% c("Bulk","MetaModel")))] <- paste0("Single-Cell ",AUCs$Model[which(!(AUCs$Model %in% c("Bulk","MetaModel")))])
# mycol <- ggsci::pal_d3("category10")(8)
# mycolmap <- c("Single-Cell Tumor" = mycol[1],
#               "Single-Cell Microglia" = mycol[2],
#               "Single-Cell Tcell" = mycol[3],
#               "Single-Cell Endothelial" = mycol[4],
#               "Single-Cell Mesenchyme" = mycol[5],
#               "MetaModel" = mycol[6],
#               "Bulk" = mycol[8],
#               "Single-Cell Immune" = mycol[7])
# for (i in levels(factor(AUCs$Class))) {
#   print(
#       ggplot(dat = AUCs[which(AUCs$Class == i),],
#              mapping = aes(fill = Model,
#              y = AUC, x = reorder(Model,AUC)))+
#         geom_bar(stat = 'identity', position = 'dodge') +
#         ggtitle(paste0("Class ", i, " AUCs"), ) +
#         scale_fill_manual(values = mycolmap) +
#         theme_bw() +
#         theme(legend.position = "none",
#               plot.title = element_text(size = 25),
#               axis.text.y = element_text(size =20),
#               axis.text.x = element_text(size =15),
#               axis.title.x = element_text(size =25),
#               axis.title.y = element_blank())+
#         coord_flip()
# )
# }
