#'  Get Flybase annotation ID file given the version
#'
#' The function takes Flybase version as an input, and download the the Flybase annotation ID file. 
#' Then if filter out non-melanogaster genes and copy information from annotation.ID to secondary.annotation.ID if latter is empty
#' 
#' The function accesses the FlyBase FTP site, so requires internet-connection.

#' @param x a vector of gene symbols.
#' @param version FlyBase ID version.
#' @keywords flybase
#' @export
#' @examples
#' GeneConverter(x, version="FB2016_04")
GetFlybaseAnnotation <- function(version) {
  temp <- tempfile()
  download.file(paste("ftp://ftp.flybase.net/releases/", version, "/precomputed_files/genes/fbgn_annotation_ID.tsv.gz", sep=""), temp, method="curl", quiet = T)
  fbgn.annotation.ID.table <- read.delim( gzfile(temp), blank.lines.skip = T, comment.char = "#", header=F, stringsAsFactors = F)
  fbgn.annotation.ID.table <- fbgn.annotation.ID.table[ !grepl("\\", fbgn.annotation.ID.table$V1, fixed=T), ] # removing non-melanogastor genes
  colnames(x = fbgn.annotation.ID.table) <- c("gene.symbol", "primary.FBgn", "secondary.FBgn", "annotation.ID", "secondary.annotation.ID")
  # By merging on the annotation.ID were are missing some genes (non-coding genes e.g.: noe) because they have different annotation ID in newer Flybase versions (CG -> CR)
  # Solution to copy annotation.ID to secondary.annotation.ID if secondary.annotation.ID is empyt and then merge on secondary.annotation.ID
  fbgn.annotation.ID.table$secondary.annotation.ID[fbgn.annotation.ID.table$secondary.annotation.ID==""] <- fbgn.annotation.ID.table$annotation.ID[fbgn.annotation.ID.table$secondary.annotation.ID==""]
  
  unlink(temp, force=T)
  return (fbgn.annotation.ID.table)
}

#'  Converting Gene Symbols or FlyBase IDs to a certain version Flybase version.
#'
#' The function takes Gene Symbols or FlyBase IDs (e.g. FBgn#######) as an input, and converts it into certain versions of IDs.
#' The function accesses the FlyBase FTP site, so requires internet-connection.
#' FlyBase IDs for genes that are split into multiple genes will be concatenated with two colons (::).  Genes that does not have matching IDs will be shown as "unknown".

#' @param x a vector of gene symbols.
#' @param version FlyBase ID version for the updated result. It should be either version numbers (e.g. 6.12) or FlyBase Release Numbers (e.g. FB2017_04).  The default is "current".  
#' @keywords flybase
#' @export
#' @examples
#' GeneConverter(x, version=6.12)
#' GeneConverter(x, version="FB2016_04")

GeneConverter <- function(x, version, verbose = F) {
  if ( !requireNamespace("readr", quietly = T) ) {
    stop("readr should be installed.", call. = F)
  }
  
  if ( !requireNamespace("dplyr", quietly = T) ) {
    stop("dplyr should be installed.", call. = F)
  }
  
  if ( !requireNamespace("tidyr", quietly = T) ) {
    stop("tidyr should be installed.", call. = F)
  }
  
  #########################################################
  # Detect the Flybase release version from gene signature
  #########################################################
  
  # Download the conversion table
  library(readr)
  all.conversion.tables <- read_tsv(file = "https://raw.github.com/mase5/flybaseR/master/data/20190604_flybase_all_releases_conversion_table.tsv.gz"
                                    , col_names = T
                                    , quote = ''
                                    , progress = verbose)
  
  # Detect the most likely Flybase release from the given gene list x
  library(dplyr)
  library(tidyr)
  coverage.table <- all.conversion.tables %>% 
    group_by(Flybase.version) %>% 
    summarize(coverage = sum(x %in% gene.symbol)/length(x = x)) %>%
    arrange(desc(coverage))
  
  tmp <- strsplit(x = coverage.table$Flybase.version[1], split = "_", fixed = 2)[[1]]
  source.fb.release.version <- paste0(tmp[2],"_",tmp[3])
  
  message(paste0("#################################################"))
  message(paste0("# Converting genes from ", source.fb.release.version, " --> ", version))
  message(paste0("#################################################"))
  
  #########################################
  # Get the gene ID list file from FlyBase
  #########################################
  
  
  if(verbose) {
    message(paste0("[SOURCE] - Downloading fbgn_annotation_ID file from Flybase (version ", source.fb.release.version,").."))
  }
  source.fbgn.annotation.ID.table <- GetFlybaseAnnotation(version = source.fb.release.version)
  if(verbose) {
    message(paste0("[TARGET] - Downloading fbgn_annotation_ID file from Flybase (version ", version,").."))
  }
  target.fbgn.annotation.ID.table <- GetFlybaseAnnotation(version = version)
  
  # If we don't unlist, some genes are still escaping while they shouldn't (e.g.: Piezo)
  # Therefore we will unlist the secondary annotation ID values for each row in the Flybase annotation ID table (both for source and target) 
  # and then merge with the secondary annotation ID
  source.fbgn.annotation.ID.table.ul <- do.call(what = "rbind", args = apply(X = source.fbgn.annotation.ID.table, MARGIN = 1, FUN = function(row) {
    row <- as.data.frame(x = t(x = row), stringsAsFactors = F)
    return (data.frame(row[, -c(ncol(x = row))], secondary.annotation.ID=unlist(strsplit(as.character(row[["secondary.annotation.ID"]]),",")), row.names = NULL, stringsAsFactors = F))
  }))
  target.fbgn.annotation.ID.table.ul <- do.call(what = "rbind", args = apply(X = target.fbgn.annotation.ID.table, MARGIN = 1, FUN = function(row) {
    row <- as.data.frame(x = t(x = row), stringsAsFactors = F)
    return (data.frame(row[, -c(ncol(x = row))], secondary.annotation.ID=unlist(strsplit(as.character(row[["secondary.annotation.ID"]]),",")), row.names = NULL, stringsAsFactors = F))
  }))
  
  #########################################
  # Create the conversion table
  #########################################
  conversion.table <- merge(x = source.fbgn.annotation.ID.table.ul
                            , y = target.fbgn.annotation.ID.table.ul
                            , by = "secondary.annotation.ID"
                            , suffixes = c(".source", ".target"))
  gene.signature.mask <- conversion.table$gene.symbol.source %in% x
  
  conversion.table.filtered <- conversion.table %>% 
    filter(gene.signature.mask) %>%
    filter(!duplicated(gene.symbol.source))
  
  genes.recovered <- as.data.frame(conversion.table.filtered)$gene.symbol.target
  genes.not.recovered <- x[!(x %in% conversion.table.filtered$gene.symbol.source)]
  n.genes.recovered <- length(x = genes.recovered)
  
  genes.not.recovered.final <- c()
  
  message("Genes that were not recovered: \n")
  # Further investigate the genes not recovered
  for(i in genes.not.recovered) {
    source.annotation.ID <- source.fbgn.annotation.ID.table$annotation.ID[source.fbgn.annotation.ID.table$gene.symbol==i]
    target.annotation.ID <- target.fbgn.annotation.ID.table$annotation.ID[target.fbgn.annotation.ID.table$gene.symbol==i]
    if(length(x = target.annotation.ID) == 0) {
      genes.not.recovered.final <- c(genes.not.recovered.final, i)
    } else {
      genes.recovered <- c(genes.recovered, i)
    }
  }
  
  n.genes.recovered <- length(x = genes.recovered)
  n.genes.not.recovered <- length(x = x) - n.genes.recovered
  percent.genes.recovered <- n.genes.recovered/length(x = x)*100
  
  message(paste0("Total number of genes recovered: ", n.genes.recovered, " (", round(x = percent.genes.recovered, digits = 2), "). \nTotal genes that were not recovered (", n.genes.not.recovered,"): \n"))
  for(i in genes.not.recovered.final) {
    message(paste0(i, " (possibly withdrawn)"))
  }
  
  return (genes.recovered)
}