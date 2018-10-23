#' Tabulate the shared TCRs chains in pairs of single-cell libraries
#'
#' Aggregate a set of TCR chain matches to generate a table with pairs of libraries and the number
#' of chains shared between the libraries.
#' 
#' @param tcr_chain_matches a data frame containing the TCR chain matches. Should contain, at minimum, the identifiers for the two TCRs. Often the output of \code{match_TCR_chains}.
#' @param tcr1_col,tcr2_col numbers or names of the columns containing the identifiers
#' @export
#' @return A data frame containing the two TCR identifiers and the number of chains shared between them.
#' 
#' @usage \code{tabulate_shared_TCR_chains(
#'      tcr_chain_matches,
#'      tcr1_col="tcr1", tcr2_col="tcr2")}
tabulate_shared_TCR_chains <-
  function(tcr_chain_matches,
           tcr1_col="tcr1", tcr2_col="tcr2") {
    assert(
      check_data_frame(tcr_chain_matches)
    )
    
    # convert numeric column identifiers to column names
    for (i in c("tcr1_col", "tcr2_col"))
      if (is.numeric(get(i))) assign(i, colnames(tcrs)[get(i)])
    
    # extract all the matches
    tcr_chains_shared <-
      aggregate(
        list(num_shared_chains=rep(1, nrow(tcr_chain_matches))),
        tcr_chain_matches[,c(tcr1_col, tcr2_col)], length)
    
    # make the order consistent, and drop the duplicates
    for (i in 1:nrow(tcr_chains_shared)) {
      tcr_chains_shared[i, c(tcr1_col, tcr2_col)] <-
        tcr_chains_shared[i, c(tcr1_col, tcr2_col)] %>%
        unlist() %>%
        sort()
    }
    tcr_chains_shared <- unique(tcr_chains_shared)
    
    tcr_chains_shared
  }
