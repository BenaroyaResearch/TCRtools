#' Write out TCR junction sequences to a FASTA file for BLAST query
#'
#' Write out TCR junctions to a FASTA file, to enable submission for BLAST query against known
#' protein sequences. This allows identification of known TCRs from published sequences.
#' @param tcrs data frame, the TCrs to be included. Must include all the columns specified in \code{cols_for_name} and \code{junction_col}.
#' @param filename character string, the path and file name to output. Defaults to "tcr_junctions.fasta"
#' @param pos_control character string containing a positive control TCR sequence that will be prepended to the fasta file. Used to ensure that BLAST search parameters are working properly. Can be set to NULL to not include positive control. The default value should give a perfect match to protein accession 1OGA_D.
#' @param sample_col character, the name of the column in \code{tcrs} to be included in the fasta query name. Separate from \code{cols_for_name} to enable easier editing.
#' @param cols_for_name character vector, the columns in \code{tcrs} to be included in the fasta query name. Values will be concatenated with "_" as separator to generate the query name. To use alignment results with \code{read_BLAST_TCR_align_hits} with length thresholding, include the junction at the end of the query name.
#' @param junction_col character string, the name of the column in \code{tcrs} to use for the CDR3 junction sequence. Must be included in \code{tcrs}.
#' @param unique_only logical, whether to include only unique combinations of \code{cols_for_name} in the output.
#' @export
#' @usage
#' write_TCR_junction_fasta(
#'   tcrs, filename="tcr_junctions.fasta",
#'   pos_control=">flu_1_TRBV_CAGAGSQGNLIF CAGAGSQGNLIF",
#'   sample_col="libid",
#'   cols_for_name=c(sample_col, "cln_count", "v_gene", "j_gene", "junction"),
#'   junction_col="junction",
#'   unique_only=TRUE)
#' @details
#' This function outputs a file designed to be used for a BLAST query on
#' https://blast.ncbi.nlmh.nih.gov. On that site, use protein blast (blastp), under "Choose File",
#' upload the file that is output by this function. Make sure that the non-redundant protein sequences
#' database (nr) is selected, set organism to "Homo sapiens (taxid:9606)", and under Algorith parameters
#' set Expect threshold to 1000 (to account for the length of the sequences being queried). If you
#' included the default pos_control, you should see it match to protein accession 1OGA_D.
write_TCR_junction_fasta <-
  function(tcrs, filename="tcr_junctions.fasta",
           pos_control=">flu_1_TRBV_CAGAGSQGNLIF CAGAGSQGNLIF",
           sample_col="libid",
           cols_for_name=c(sample_col, "cln_count", "v_gene", "j_gene", "junction"),
           junction_col="junction",
           unique_only=TRUE) {

    if (!is.data.frame(tcrs)) stop("Input tcrs object must be a data frame")
    tcrs <- as.data.frame(tcrs)

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
