#' Read the alignment output hit table .txt file from a TCR sequence BLAST query
#'
#' Read in the alignment hit table file output by BLAST following a query of TCR sequences.
#' This function reads and reformats the alignment results for downstream use. It expects the format
#' output by the BLAST Download > Alignment > Hit Table(text).
#' 
#' @param filename character string, the path and file name to read.
#' @param identity_threshold numeric, the threshold for \% identity match between the query and the subject sequences. Defaults to 100, can be reduced to allow less strict matching.
#' @param length_threshold numeric, the threshold for \% of the total length of the query sequence that must be found in the alignment. Defaults to 100, can be reduced to allow less strict matching.
#' @param query_seq_in_name logical, whether the alignment query sequence is included in the fasta query name. Alignment length thresholding requires the query sequence, which is assumed to be at the end of the query name.
#' @export
#' @return a data frame containing the alignment hits, optionally filtered by identity and length thresholds
#' @usage \code{
#'   read_BLAST_TCR_align_hits(
#'     filename, identity_threshold=100, length_threshold=100,
#'     query_seq_in_name=TRUE)}
read_BLAST_TCR_align_hits <-
  function(
    filename, identity_threshold=100, length_threshold=100,
    query_seq_in_name=TRUE) {
    col_names <- 
      grep("Fields", readLines(filename),
           value=TRUE)[1] %>%
      str_replace("# Fields: ", "") %>%
      str_split(pattern=", ", simplify=TRUE) %>%
      str_replace_all("%", "perc") %>%
      standardize_names()
    hits <-
      read.table(filename, sep="\t", col.names=col_names, quote="")
    hits <- # make identity threshold cut
      hits[hits[,"perc_identity"] >= identity_threshold,]
    
    # make alignment length threshold cut
    if (query_seq_in_name) {
      hits[,"query_sequence"] <-
        str_extract(hits[,"query_acc_ver"], "[A-Z]+$")
      hits <- 
        hits[
          with(hits,
               (alignment_length / nchar(query_sequence) * 100) >= length_threshold),]
    } else {
      warning("Cannot apply alignment length thresholding without query sequence.")
    }
    
    hits
  }
