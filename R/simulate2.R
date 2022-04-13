#'@title Simulate multivariate GWAS data
#'@param N Either a single number representing GWAS sample size for all studies or an K+2
#'vector giving sample size for X, Y, and then each of Z_1 to Z_K
#'@param J Number of variants to simulate
#'@param taux_xz,tau_yz Length K vectors. Effects between each of Z_1, ..., Z_K and X and Y respectively.
#'@param dir_xz,dir_yz Direction of effects between each of Z_1, ..., Z_K and X and Y respectively.
#' A value of 1 indicates an effect from Z_k to X. A value of -1 indicates an effect from X to Z_k.
#'@param h2_x,h2_y Scalars giving hertiabiltiy of each of X and Y.
#'@param h2_z Length K vector giving heritability of each of Z_1 to Z_K.
#'@param pi_x,pi_y Scalrs between 0 and 1 giving the proportion of non-zero effects for each of X and Y
#'@param pi_z Length K vector giving proporiton of non-zero effects for each of Z_1 to Z_K.
#'@param gamma Causal effect of X on Y.
#'@details This function simulates data for traits X, Y, and Z_1, ..., Z_k.
#'tau_xz and tau_yz give effects of Z_1, ..., Z_K on (or from) X and Y respectively in units of proportion of variance
#'explained. For example, if
#'tau_xz = c(0.2, 0.15, 0.15), tau_yz = c(0.3, 0.1, 0.4), dir_xz = c(1, -1, -1), dir_yz = c(1, 1, -1)
#'then there are three variables Z_1, Z_2, Z_3. Z_1 is a common cause of X and Y
#'(indicated by a 1 in both dir_xz and dir_yz).
#'Z_1 explains 20% of the variance of X and 30% of the variance of Y.
#'Z_2 is a mediator between X and Y. There is an effect of X on Z_2 (the second entry of dir_xz is -1) explainng 15% of
#'the variance of Z_1. There is also an effect of Z_1 on Y explaining 10% of the variance of Y.
#'Z_3 is a common child of X and Y. X explains 15% of the variation of Z_3 while Y explains 40% of the variance of Z_3.
#'@export
sim_mv <- function( tau_xz, tau_yz, dir_xz, dir_yz, gamma,
                    h2_x, h2_y, h2_z, N, J,
                    pi_x, pi_y, pi_z){
  p <- length(tau_xz)
  stopifnot(length(tau_yz) == p)
  stopifnot(length(dir_xz) == p)
  stopifnot(length(dir_yz) == p)
  stopifnot(length(h2_z) == p)
  stopifnot(length(pi_z) == p)
  #stopifnot(all(tau_yz >= 0) & all(tau_xz >= 0))

  if(any(dir_xz == 1 & dir_yz == -1)){
    stop("No cycles allowed")
  }

  tau_xz <- sqrt(abs(tau_xz))*sign(tau_xz)
  tau_yz <- sqrt(abs(tau_yz))*sign(tau_yz)
  gamma <- sqrt(abs(gamma))*sign(gamma)

  # Direct Effects
  G <- matrix(0, nrow = p+2, ncol = p+2)
  G[1,2] <- gamma
  G[which(dir_xz == 1) + 2, 1] <- tau_xz[dir_xz ==1]
  G[1, which(dir_xz == -1) + 2] <- tau_xz[dir_xz == -1]

  G[which(dir_yz == 1) + 2, 2] <- tau_yz[dir_yz ==1]
  G[2, which(dir_yz == -1) + 2] <- tau_yz[dir_yz == -1]

  G_t <- direct_to_total(G)


  # Compute Input and Direct Heritability
  h2 <- c(h2_x, h2_y, h2_z)
  C <- colSums(h2*(G_t)^2)
  H <- t(G_t^2)
  input_h2 <- solve(diag(p+2) + H) %*% matrix(C, ncol = 1) %>% as.vector()

  direct_h2 <- h2-input_h2
  if(any(direct_h2 < 0)) stop("Input variance greater than 1 for at least one variable.")
  G_t <- G_t*sqrt(direct_h2)
  diag(G_t) <- sqrt(direct_h2)

  F_mat <- t(G_t)

  dat <- sim_sumstats_lf(F_mat = F_mat,
                         N = N, J = J, h_2_trait = c(h2_x, h2_y, h2_z),
                         omega = rep(1, p + 2),
                         pi_L = c(pi_x, pi_y, pi_z),
                         overlap_prop = 0,
                         h_2_factor = rep(1, p+2),
                         pi_theta = 1)
  dat$F_mat_init <- F_mat
  dat$total_effects <- t(dat$F_mat)/diag(dat$F_mat)
  dat$direct_effects <- G
  diag(dat$direct_effects) <- 1
  dat$B <- dat$Z * dat$se_beta_hat
  return(dat)
}

path_one_step <- function(graph, paths){
  paths_new <- map(paths, function(p){
    start <- p[1]
    into <- graph[,start]
    if(all(into == 0)) return(NULL)
    new_paths <- lapply(which(into!=0), function(i){ c(i, p)})
    return(new_paths)
  })
  paths_new <- unlist(paths_new, recursive = FALSE)
  return(paths_new)
}

all_paths <- function(graph, start, end){
  p <- list(c(end))
  my_paths <- list()
  while(length(p)> 0){

    p <- path_one_step(graph, p)

    starts <- map(p, 1) %>% unlist()
    if(any(starts == start)){
      ix <- which(starts == start)
      my_paths <- append(my_paths, p[ix])
      p <- p[-ix]
    }
  }
  return(my_paths)
}

get_total_effect <- function(graph, start, end){
  paths <- all_paths(graph, start, end)
  if(length(paths) == 0) return(0)
  map(paths, function(x){
    e <- sapply(seq_along(x[-length(x)]), function(i){graph[x[i], x[i+1]]})
    return(prod(e))
  }) %>% unlist() %>%
    sum() %>%
    return()
}

direct_to_total <- function(graph){
  n = nrow(graph)
  G_total <- expand.grid(start = seq(n), end = seq(n))
  G_total$value <- map2(G_total$start, G_total$end, function(x, y){
    get_total_effect(graph, x, y)
  }) %>% unlist()
  G_total <- G_total %>% pivot_wider(names_from = end, values_from = value) %>%
    select(-start) %>% as.matrix()
  return(G_total)
}
