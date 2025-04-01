#' Retrieve genes corresponding to the organism of interest
#'
#' @param count_matrix Matrix or dgCMatrix containing the counts from
#'  single cell RNA sequencing data.
#' @param organism String specifying the organism of origin of the input
#'  data where there are two standard options,
#'
#'  * "Hsap"
#'  * "Mmus"
#'
#'  If a user wishes to use a non-standard organism they must input a list
#'  containing strings for the patterns to match mitochondrial and ribosomal
#'  genes of the organism. If available, nuclear-encoded genes that are likely
#'  retained in the nucleus, such as in nuclear speckles, must also
#'  be specified. An example for humans is below,
#'
#'  * organism = c(mito_pattern = "^MT-",
#'                 ribo_pattern = "^(RPS|RPL)",
#'                 nuclear <- c("NEAT1","XIST", "MALAT1")
#'
#' * Default is "Hsap"
#' @return List containing the indices of the count matrix corresponding to
#'  mitochondrial, non-mitochondrial, and ribosomal gene sets.
#' @keywords internal
get_organism_indices <- function(
    count_matrix,
    organism
){
  if (organism == "Hsap") {
    mito_pattern <- "^MT-"
    ribo_pattern <- "^(RPS|RPL)"
    nuclear <- c("FIRRE", "NEAT1","XIST", "MALAT1", "MIAT",
                 "MEG3", "KCNQ1OT1", "HOXA11-AS", "FTX")
    MALAT1 <- "MALAT1"
  }

  if (organism == "Mmus") {
    mito_pattern <- "^mt-"
    ribo_pattern <- "^(rps|rpl)"
    nuclear <- c("Firre", "Neat1","Xist", "Malat1", "Miat",
                 "Meg3", "Kcnq1ot1", "Hoxa11-as", "Ftx")
    MALAT1 <- "Malat1"
  }

  # Allow for user specification for non-standard organism
  if (!organism %in% c("Hsap", "Mmus")) {
    mito_pattern <- organism$mito_pattern
    ribo_pattern <- organism$ribo_pattern
    nuclear <- organism$nuclear
    MALAT1 <- organism$MALAT1
  }

  # Isolate gene set indices (consistent across cells, not subsetting the matrix)
  mito_idx <- grep(mito_pattern, rownames(count_matrix), ignore.case = FALSE)
  nucl_idx <- which(rownames(count_matrix) %in% nuclear)
  mito_idx <- c(mito_idx, nucl_idx)
  non_mito_idx <- setdiff(seq_len(nrow(count_matrix)), mito_idx)
  ribo_idx <- grep(ribo_pattern, rownames(count_matrix), ignore.case = FALSE)

  return(list(
    mito_idx = mito_idx,
    ribo_idx = ribo_idx,
    non_mito_idx = non_mito_idx,
    mito_pattern =  mito_pattern,
    ribo_pattern = ribo_pattern,
    MALAT1 = MALAT1
  ))

}

utils::globalVariables(c(
  "Features", "New_Features", "New_MitoProp", "New_RiboProp", "RiboProp",
  "Original_Features", "Original_MitoProp", "Original_RiboProp", "Ribo. prop"
))

