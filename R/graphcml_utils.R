#' Pulled from graph-cML Lin2023
total_to_dir_projection_gcml <- function(G_tot, maxit = 100) {
  iter <- 0
  G_dir <- G_tot %*% solve(diag(nrow(G_tot))+G_tot)
  G_dir0 <- G_dir
  converge <- 0
  while(any(abs(diag(G_dir)) > 1e-4) & iter < maxit){
    G_dir <- G_tot %*% solve(diag(nrow(G_tot)) + G_tot)
    GG <- G_dir * t(G_tot)
    diag(GG) <- 0
    diag(G_tot) <- rowSums(GG)
    iter <- iter + 1;
  }
  if(iter>=maxit){
    G_dir <- G_dir0
    converge <- 1
    warning("Did not converge after", maxit)
  }
  G_dir
}
