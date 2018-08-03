#' Write out TCR junction sequences to a FASTA file for BLAST query
#'
#' Write out TCR junctions to a FASTA file, to enable submission for BLAST query against known
#' protein sequences. This allows identification of known TCRs from published sequences.
#' @param tcrs data frame, the TCrs to be included. Must include all the columns specified in \code{cols_for_name} and \code{junction_col}.
#' @param filename character string, the path and file name to output. Defaults to "tcr_junctions.fasta"
#' @param pos_control character string containing a positive control TCR sequence that will be prepended to the fasta file. Used to ensure that BLAST search parameters are working properly. Can be set to NULL to not include positive control.
#' @param cols_for_name character vector, the columns in \code{tcrs} to be included in the fasta query name. Values will be concatenated with "_" as separator to generate the query name. To use alignment results with \code{read_BLAST_TCR_align_hits} with length thresholding, include the junction at the end of the query name.
#' @param junction_col character string, the name of the column in \code{tcrs} to use for the CDR3 junction sequence. Must be included in \code{tcrs}.
#' @param unique_only logical, whether to include only unique combinations of \code{cols_for_name} in the output.
#' @export
#' @usage \code{
#' write_TCR_junction_fasta(
#'   tcrs, filename="tcr_junctions.fasta",
#'   pos_control=">flu_1_TRBV_CAGAGSQGNLIF CAGAGSQGNLIF",
#'   cols_for_name=c("sample", "cln_count", "v_gene", "j_gene", "junction"),
#'   junction_col="junction",
#'   unique_only=TRUE)}
write_TCR_junction_fasta <-
  function(tcrs, filename="tcr_junctions.fasta",
           pos_control=">flu_1_TRBV_CAGAGSQGNLIF CAGAGSQGNLIF",
           cols_for_name=c("sample", "cln_count", "v_gene", "j_gene", "junction"),
           junction_col="junction",
           unique_only=TRUE) {
    
    # check for junction_col in tcrs object
    if (!(junction_col %in% colnames(tcrs))) {
      stop("Junction column '", junction_col, "' not found in tcrs object.")
    }
    # check if specified columns are present in TCRs object
    if (!all(cols_for_name %in% colnames(tcrs))) {
      warning(
        "Some columns specified for sample name were not found in tcrs object. ",
        "The following columns have been dropped:",
        paste(cols_for_name[!(cols_for_name %in% colnames(tcrs))]))
      cols_for_name <- cols_for_name[cols_for_name %in% colnames(tcrs)]
    }
    
    tcrs <- tcrs[,union(cols_for_name, junction_col)]
    if (unique_only) tcrs <- tcrs[!duplicated(tcrs[,cols_for_name]),]
      
  sink(filename); on.exit(sink())
  cat(stringr::str_replace(pos_control, " ", "\r\n"), "\r\n", sep="")
  for (i in 1:nrow(tcrs))
    cat(">", paste(tcrs[i,cols_for_name], collapse="_"), "\r\n",
        tcrs[i,junction_col], "\r\n", sep="")
  
}
