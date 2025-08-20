#' Tabulate the shared TCRs chains in pairs of single-cell libraries
#'
#' Aggregate a set of TCR chain matches to generate a table with pairs of libraries and the number
#' of chains shared between the libraries.
#'
#' @param tcr_chain_matches a data frame containing the TCR chain matches. Should contain, at minimum, the identifiers for the two TCRs. Often the output of \code{match_TCR_chains}.
#' @param tcr1_col,tcr2_col numbers or names of the columns containing the identifiers
#' @importFrom stats aggregate na.omit
#' @importFrom utils read.table
#' @importFrom rlang sym
#' @import dplyr
#' @export
#' @return A data frame containing the two TCR identifiers and the number of chains shared between them.
#'
#' @usage
#' tabulate_shared_TCR_chains(
#'      tcr_chain_matches,
#'      tcr1_col = "tcr1", tcr2_col = "tcr2")
tabulate_shared_TCR_chains <-
  function(tcr_chain_matches,
           tcr1_col = "tcr1", tcr2_col = "tcr2") {
    checkmate::assert(
      checkmate::check_data_frame(tcr_chain_matches)
    )

    if (nrow(tcr_chain_matches) == 0) {
      tcr_chains_shared <-
        cbind(
          tcr_chain_matches[, c(tcr1_col, tcr2_col)],
          "num_shared_chains" = integer())
    } else {
      # convert numeric column identifiers to column names
      for (currTcrColName in c("tcr1_col", "tcr2_col"))
        if (is.numeric(get(currTcrColName))) assign(currTcrColName, colnames(tcr_chain_matches)[get(currTcrColName)])

      # extract all the matches
      tcr_chains_shared <-
        dplyr::summarise(
          dplyr::group_by(tcr_chain_matches, !!rlang::sym(tcr1_col), !!rlang::sym(tcr2_col)),
          num_shared_chains = n())
      tcr_chains_shared <- as.data.frame(tcr_chains_shared)
          
      # # make the order of the cells consistent, and drop the duplicates
      # for (currRow in seq_len(nrow(tcr_chains_shared))) {
      #   tcr_chains_shared[currRow, c(tcr1_col, tcr2_col)] <-
      #     tcr_chains_shared[currRow, c(tcr1_col, tcr2_col)] %>%
      #     unlist() %>%
      #     sort()
      # }
      # tcr_chains_shared <- distinct(tcr_chains_shared)
      
      # drop duplicates where tcr1_col == tcr2_col, alternate approach
      tcr_chains_shared <-
        mutate(rowwise(tcr_chains_shared),
               tcrs_combined = paste(sort(c(!!rlang::sym(tcr1_col), !!rlang::sym(tcr2_col))), collapse = ""))
      tcr_chains_shared <- as.data.frame(ungroup(tcr_chains_shared))
      tcr_chains_shared <-
        tcr_chains_shared[!duplicated(tcr_chains_shared$tcrs_combined),]
      tcr_chains_shared$tcrs_combined <- NULL
    }

    tcr_chains_shared
  }
