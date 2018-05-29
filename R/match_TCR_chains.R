#' Match the TCR junctions and (optionally) V and J genes between single-cell libraries
#'
#' Determine matching TCR chain sequences among a set of single-cell TCR sequences.
#' purpose is to show clonotype sharing among cells within and between samples. Samples can optionally
#' be colored by grouping variables, and there are multiple options for determining TCR clonality.
#' 
#' @param tcrs a data frame containing the TCR sequences. Should include an identifier column, plus columns for junction and optionally the V and J genes.
#' @param id_col number or name of the column containing the identifiers
#' @param junction_col nmber or name of the column containing the TCR junctions
#' @param match_V_gene,match_J_gene logical, whether to include the V/J genes when matching TCR chains. Default to TRUE. If set to FALSE, the V/J gene identities are ignored when matching chains.
#' @param V_gene_col,J_gene_col nmber or name of the columns containing the V and J genes. Ignored if match_V_gene and/or match_J_gene are set to FALSE.
#' @export
#' @return A data frame containing the TCR chain information (junction and optionally V and J genes) and the identifiers for the two cells with matching chains. Identifier column names in the returned data frame are "tcr1" and "tcr2".
#' 
#' @usage \code{
#' match_TCR_chains(
#'      tcrs,
#'      id_col="libid",
#'      junction_col="junction"
#'      match_V_gene=TRUE, match_J_gene=TRUE,
#'      V_gene_col="v_gene", J_gene_col="j_gene")}
match_TCR_chains <-
  function(tcrs,
           id_col="libid",
           junction_col="junction",
           match_V_gene=TRUE, match_J_gene=TRUE,
           V_gene_col="v_gene", J_gene_col="j_gene") {
    
    # convert numeric column identifiers to column names
    for (i in c("id_col", "junction_col", "V_gene_col", "J_gene_col"))
      if (is.numeric(get(i))) assign(i, colnames(tcrs)[get(i)])
    
    # define columns to match chains by
    match_cols <- junction_col
    if (match_V_gene) match_cols <- c(match_cols, V_gene_col)
    if (match_J_gene) match_cols <- c(match_cols, J_gene_col)
    
    # cut to only ID, junction, and (optionally) V and/or J gene columns
    tcrs <- tcrs[,c(id_col, match_cols)]
    
    # remove duplicates (generally due to issues in mixcr)
    tcrs <- unique(tcrs)
    
    # match cells by shared junctions
    tcrs.match <- merge(tcrs, tcrs, by=match_cols)
    
    # remove self-matches
    tcrs.match <-
      tcrs.match[tcrs.match[,paste0(id_col, ".x")] != tcrs.match[,paste0(id_col, ".y")],]
    
    # rename identifer columns
    colnames(tcrs.match)[colnames(tcrs.match)==paste0(id_col, ".x")] <- "tcr1"
    colnames(tcrs.match)[colnames(tcrs.match)==paste0(id_col, ".y")] <- "tcr2"
    
    tcrs.match
  }
