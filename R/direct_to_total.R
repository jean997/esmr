#' Solves for total effects from direct effects using (I - B_dir)^{-1} - I
#'
#' @param B_dir
#' @param restrict_dag
#'
#' @return
#' @export
#'
#' @examples
direct_to_total <- function(B_dir, restrict_dag = TRUE) {
  n <- nrow(B_dir)
  B_total <- solve(diag(n) - B_dir) - diag(n)
  if(!all(diag(B_total) == 0) && restrict_dag){
    stop("Failed to compute total effects from direct. Check that supplied B_dir corresponds to a valid DAG.\n")
  }
  return(B_total)
}

#' Solves for total effects from direct effects using I - (I + B_tot)^{-1}
#' Assumes that B_tot is a valid adjacency matrix for a DAG and that spectral radius is less than one
#'
#' @param B_tot
#' @param restrict_dag Should the function fail if B_tot is not a DAG? Otherwise, will check that spectral radius max(abs(eigenvalues)) < 1
#'
#' @return
#' @export
total_to_direct <- function(B_tot, restrict_dag = TRUE){
  n <- nrow(B_tot)
  B_dir <- diag(n) - solve(diag(n) + B_tot)
  if(!all(diag(B_tot) == 0)){
    stop("Failed to compute total effects from direct. Check that supplied B_tot corresponds to a valid DAG.\n")
  }
  return(B_dir)
}

#' Compute direct to total effect adjacency matrix
#' Does not assume that B is a valid DAG
#'
#' @param B Direct effect adjacency matrix
#' @param remove_diag Logical indicating whether to remove the diagonal
#' @examples
#' X <- structure(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0), dim = c(4L,4L))
#' print(direct_to_total_adj(X))
direct_to_total_adj <- function(B, remove_diag = TRUE) {
  B_tot <- B
  B_step <- B
  p <- ncol(B)
  for (i in 2:p) {
    B_step <- B_step %*% B
    # Logical OR to keep all previous edges plus n-step transition edges
    B_tot_next <- (B_tot | B_step) + 0
    if (all(B_tot_next == B_tot)) {
      # No new edges added, so we are done
      break
    }
    B_tot <- B_tot_next
  }
  # Remove the diagonal
  if (remove_diag) diag(B_tot) <- 0
  return(B_tot)
}
