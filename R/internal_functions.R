##########################
### INTERNAL FUNCTIONS ###
##########################

# ----------------------------------------------------------------------------------------

#' @title .corum.reshape
#'
#' @description
#' Reshapes the CORUM database downloads files for gene entrenchment analyses
#'
#' @param corum_table Either the path to the CORUM database file or a data.frame with the same structure.
#'
#' @return Data.frame with 4 columns: complex.id, complex.name, organism, protein.members.
#'
#' @importFrom data.table fread
#' @importFrom stringr str_split
#'
#' @keywords internal


.corum.reshape =
  function(corum_table) {

    ## import annotations
    if ("data.frame" %in% class(corum_table)) {
      corum = as.data.frame(corum_table) }
    else {
      corum = data.table::fread(corum_table, data.table = FALSE)
    }

    ## reshape database
    members = stringr::str_split(corum$subunits_gene_name, pattern = ";")
    n = lengths(members)

    corum_list = data.frame(complex.id = rep(corum$complex_id, n),
                            complex.name = rep(corum$complex_name, n),
                            organism = rep(corum$organism, n),
                            protein.members = unlist(members))

    return(corum_list)

  } # END .corum.reshape


# ----------------------------------------------------------------------------------------
