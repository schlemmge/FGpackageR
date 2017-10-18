
#' Read a count matrix from disk.
#' 
#' @param file Name of file containing the count matrix.
#' @param path Path to count matrix file.
#' @param row.names Column number for rownames of count matrix; defaults to 1.
#' @param header Logical indicating whether header line is present in count matrix file; defaults to TRUE.
#' @param sep Column separator for count matrix file; defaults to "\t".
#' @param ... Further arguments passed to read.table function that is used for reading the count matrix.
#' @return Count data as sparse Matrix object.
readCountMatrix <- function(file, path = NULL, row.names = 1, header = TRUE, sep = "\t", ...) {
  if (!is.null(path)) {
    currentPath <- getwd()
    setwd(path)
  }
  denseData <- read.table(file, 
                          row.names = row.names,
                          header = header,
                          ...)
  denseData <- as.matrix(denseData)
  sparseData <- Matrix(denseData, 
                       sparse = TRUE)
  if (!is.null(path)) setwd(currentPath)
  return(sparseData) 
}


#' Assign cellIds to count matrix and provide mapping between cellIds and matrix column names.
#' 
#' @param countMatrix Sparse Matrix object with read/UMI counts per gene and cell.
#' @return List with elements 'countMatrix' and 'metadata'.
assignCellIDs <- function(countMatrix) {
  returnList <- list(countMatrix = countMatrix,
                     metadata = data.frame(cellId = seq(from = 0,
                                                        to = ncol(countMatrix)-1,
                                                        by = 1),
                                           cellName = colnames(countMatrix),
                                           stringsAsFactors = FALSE))
  colnames(returnList$countMatrix) <- returnList$metadata$cellId
  return(returnList)
}


#' Prepare a metadata table from count matrix.
#' 
#' @param countMatrix Sparse Matrix object that holds read count and meta data.
#' @param lines Vector of row numbers or row names that hold metadata in the count matrix. Defaults to 1:5.
#' @return Data frame holding the transposed matrix subset.
getMetaFromMatrix <- function(countMatrix, lines = 1:5) {
  metaMatrix <- as.matrix(countMatrix[lines,])
  metaDf <- data.frame(cellId = seq(from = 0, 
                                    to = ncol(countMatrix)-1, 
                                    by = 1),
                       t(metaMatrix) )
  return(metaDf)
}


#' Prepare a metadata table from count matrix column names.
#' 
#' @param countMatrix Sparse Matrix object that holds read count and meta data.
#' @param sep Separator character for metadata in column names; defaults to "_".
#' @return Data frame with metadata extrcted from count matrix column names.
getMetaFromColnames <- function(countMatrix, sep = "_") {
  splitNames <- strsplit(colnames(countMatrix), split = sep)
  columnCount <- min( sapply(splitNames, length) )
  metaDf <- data.frame(cellId = seq(from = 0, 
                                    to = ncol(countMatrix)-1, 
                                    by = 1),
                       row.names = colnames(countMatrix))
  for (i in 1:columnCount) {
    metaDf[, paste("V", i, sep = "")] <- sapply(splitNames, "[[", i)
  }
  return(metaDf)
}


#' Map gene IDs to Entrez IDs
#' 
#' @param countMatrix Sparse Matrix object with read/UMI counts per gene and cell.
#' @param IDtype ID type of count matrix rownames; defaults to "SYMBOL".
#' @param annoDB OrgDb object for mapping of gene IDs.
#' @return List with elements 'entrez', 'non_entrez', each of which are subsets of the count matrix with uniquely mapping and multiply mapping or unmappable Entrez IDs, as well as 'log' with mapping information.
entrezMapping <- function(countMatrix, IDtype = "SYMBOL", annoDB) {
  returnList <- list(entrez = NULL,
                     non_entrez = NULL,
                     log = data.frame(originalID = row.names(countMatrix),
                                      stringsAsFactors = FALSE))
  mappedIDs <- mapIds(annoDB,
                      keys = returnList$log$originalID,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "list")
  multiMappingIndex <- sapply(mappedIDs, length)>1
  unassignedIndex <- sapply(sapply(mappedIDs, na.omit), length)==0
  returnList$log$mappedID <- sapply(mappedIDs , function(x) paste(x, collapse = " // "))
  duplicatedTargetIDs <- unique(returnList$log$mappedID[duplicated(returnList$log$mappedID[!unassignedIndex])])
  duplicatedTargetIndex <- if (length(duplicatedTargetIDs)>0) {
    returnList$log$mappedID %in% duplicatedTargetIDs
  } else {
    rep(FALSE, length(returnList$log$mappedID))
  }
  returnList$log$mappingLog <- "successfully mapped to unique Entrez ID"
  returnList$log$mappingLog[multiMappingIndex] <- "mapped to multiple Entrez IDs"
  returnList$log$mappingLog[unassignedIndex] <- "no Entrez ID defined"
  returnList$log$mappingLog[duplicatedTargetIndex] <- "multiple original IDs mapped to same Entrez ID"
  entrezIndex <- !multiMappingIndex & !unassignedIndex & !duplicatedTargetIndex
  returnList$entrez <- countMatrix[entrezIndex,]
  returnList$non_entrez <- countMatrix[!entrezIndex,]
  return(returnList)
}


#' Generate a dataset description from a minimum set of information.
#' 
#' @param manifest
#' @return Character vector containing the dataset description displayed in the FASTGenomics data store.
makeDatasetDescription <- function(manifest) {
  returnText <- paste("This dataset comprises transcriptomes from ", 
                      sampleSize, 
                      " single cells generated with ", 
                      technology, 
                      " from ", 
                      tissue, 
                      ". Results have been published in ", 
                      insertLink(paperID, paperURL), 
                      "; the data has been retrieved from ", 
                      insertLink(datasetID, datasetURL), 
                      ".",
                      sep="")
  return(returnText)
}


#' Create manifest structure for FASTGenomics data package.
#' 
#' @return List with elements required in FASTGenomics manifest file.
makeRawManifest <- function() {
  manifest <- list(schema_version = 2.1,
                   data = list(cell_metadata = list(file = "cell_metadata.tsv",
                                                    organism = 10090,
                                                    batch_column = NULL),
                               gene_metadata = list(file = "gene_metadata.tsv"),
                               expression_data = list(file = "expression_data.tsv"),
                               supplemental = list(unconsidered_genes = list(expression_data = list(file = NULL),
                                                                             gene_metadata = list(file = NULL)))),
                   metadata = list(title = "",
                                   technology = "",
                                   version = 1,
                                   contact = "",
                                   description = "",
                                   short_description = "",
                                   processing = list(notes = "",
                                                     tools = ""),
                                   image = "image.png"))
  return(manifest)
}


#' Make Markdown-formatted link to website.
#' 
#' @param url Character vector containing the link destination.
#' @param title Character vector containing the link name.
#' @return The named link in Markdown format.
#' @examples
#' makeMdlink(url = "https://fastgenomics.org",
#'            title = "FASTGenomics website")
makeMdLink <- function(url, title) {
  paste("[", title, "](", url, " \"", title, "\")", sep = "")
}


#' 
initFGpackage <- function() {
  
}