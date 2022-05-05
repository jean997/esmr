
#'@export
eb_mr <- function(beta_hat_Y, se_Y, beta_hat_X, se_X,
                  R = NULL,
                  ebnm_fn = flashier:::as.ebnm.fn(prior_family = "point_normal", optmethod = "nlm"),
                  max_iter = 100,
                  seed = 123, sigma_beta = Inf,
                  tol = 1e-7){

  set.seed(seed)

  n <- length(beta_hat_Y)
  stopifnot(length(se_Y) == n)

  if("matrix" %in% class(beta_hat_X)){
    stopifnot(nrow(beta_hat_X) == n)
    stopifnot(nrow(se_X) == n)
    p <- ncol(beta_hat_X) + 1
  }else{
    stopifnot(length(beta_hat_X) == n)
    stopifnot(length(se_X) == n)
    p <- 2
  }
  if(missing(R) | is.null(R)){
    R <- diag(p)
  }
  Y <- cbind(beta_hat_Y, beta_hat_X)
  S <- cbind(se_Y, se_X)

  s_equal <- apply(S, 2, function(x){all(x == x[1])}) %>% all()
  if(s_equal){
    s <- S[1,]
    omega <- solve(diag(s) %*% R %*% diag(s))
  }else{
    omega <- apply(S, 1, function(s){
      solve(diag(s) %*% R %*% diag(s))
    }, simplify = FALSE)
  }

  # Initialize alpha and gamma at ebnm solutions assuming beta = 0
  if(s_equal){
    l_init <- map(seq(p), function(j){
      ebnm_res <- ebnm_fn(x = Y[,j],
                          s = S[1,j],
                          g_init = NULL, fix_g= FALSE, output = output_all())

      ebnm_res$KL <-  (ebnm_res$log_likelihood
                       - flashier:::normal.means.loglik(Y[,j],
                                       S[,j],
                                       ebnm_res$posterior$mean,
                                       ebnm_res$posterior$second_moment))
      ebnm_res})
  }else{
    l_init <- map(seq(p), function(j){
      ebnm_res <- ebnm_fn(x = Y[,j],
                          s = S[,j],
                          g_init = NULL, fix_g= FALSE, output = output_all())

      ebnm_res$KL <-  (ebnm_res$log_likelihood
                       - flashier:::normal.means.loglik(Y[,j],
                                                        S[,j],
                                                        ebnm_res$posterior$mean,
                                                        ebnm_res$posterior$second_moment))
      ebnm_res})
  }
  obj <- map(l_init, "KL") %>% unlist() %>% sum()
  cat("0: ", obj, "\n")

  # Start with beta update

  lbar <- l_init %>% map(function(l){l$posterior$mean}) %>% do.call(cbind, .)
  l2bar <- l_init %>% map(function(l){l$posterior$second_moment}) %>% do.call(cbind, .)

  beta_m <- rep(0, p-1)
  i <- 1

  make_fbar <- function(beta_m, beta_s){
    p <- length(beta_m) + 1
    fbar <- f2bar <- diag(p)
    fbar[1,-1] <- beta_m
    f2bar[1,-1] <- beta_s^2 + beta_m^2
    return(list(fbar = fbar, f2bar = f2bar))
  }
  check <- 1
  while(i < max_iter & check > tol){
    # beta update
    beta_upd <- map(seq(p-1), function(k){
      update_beta_k(Y, k, lbar, l2bar, omega, beta_m, sigma_beta = sigma_beta)
    })
    beta_m <- map(beta_upd, "m") %>% unlist()
    beta_s <- map(beta_upd, "s") %>% unlist()

    ff <- make_fbar(beta_m, beta_s)
    fbar <- ff$fbar
    f2bar <-  ff$f2bar

    Ybark <- map(seq(p), function(kk){lbar[,kk,drop = FALSE] %*% t(fbar[,kk,drop = FALSE])})
    R_k <- map(seq(p), function(kk){
      Y - (Ybark[-kk] %>% do.call(`+`, .))
    })

    # alpha and gamma updates
    l_update <- map(seq(p), function(j){
      update_l_k(R_k[[j]], fbar[,j], f2bar[,j], omega, ebnm_fn)})

    lbar <- l_update %>% map(function(l){l$posterior$mean}) %>% do.call(cbind, .)
    l2bar <- l_update %>% map(function(l){l$posterior$second_moment}) %>% do.call(cbind, .)

    obj_old <- obj
    obj <- map(l_update, "KL") %>% unlist() %>% sum()
    check <- obj - obj_old
    #if(check < 0) warning("Objective increased on iteration ", i, ".\n")
    check <- abs(check)
    cat(i, ": ", obj, " ", beta_m, " ", beta_s, "\n")
    i <- i + 1
  }
  result <- list(beta_m = beta_m, beta_s = beta_s)
  return(result)
}


calc_ll <- function(Y, lbar, l2bar, fbar, f2bar, omega){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  stopifnot(ncol(fbar) ==p)


  #First col of lbar is alpha the rest are gammas
  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) != "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }

  Ybar <- lbar %*% t(fbar)
  Y2bar <- l2bar %*% t(f2bar)
  R <- Y - Ybar
  R2 <- Y^2 -2*Y*Ybar + Y2bar

  # Calculate E(log p(Y | l, f, omega))

  if(s_equal){
    ll_t1 <- (-1/2)*t(t(R2)*diag(omega)) %>% sum()
    o <- omega
    diag(o) <- 0
    ll_t2 <- (-1/2)*(R %*% o %*% t(R) %>% sum())
  }else{
    OD <- map(omega, function(o){diag(o)}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    ll_t1 <- (-1/2)*R2*OD %>% sum()
    ll_t2 <- (-1/2)*(map(seq(n), function(i){
      o <- omega[[i]]
      diag(o) <- 0
      sum(R[i,,drop = FALSE] %*% o %*% t(R[1,,drop = FALSE]))
    }) %>% unist() %>% sum())
  }

  #log likelihood minus constant depending on omega
  ll <- ll_t1 + ll_t2
}

### Testing

## Simulate Some Data


if(FALSE){
G <- matrix(c(0, sqrt(0.25), 0, 0), nrow = 2, byrow = TRUE)
colnames(G) <- row.names(G) <- c("X", "Y")
G

set.seed(100)
sim_dat1 <- sim_mv(G = G,
                   N = 80000, J = 50000,
                   h2 = c(0.3, 0.4),
                   pi = 800/50000)


res <- eb_mr(beta_hat_Y = sim_dat1$beta_hat[,2],
             se_Y = sim_dat1$se_beta_hat[,2],
             beta_hat_X = sim_dat1$beta_hat[,1],
             se_X = sim_dat1$se_beta_hat[,1],
             R = NULL)

R_E <- diag(2)
R_E[1,2] <- R_E[2,1] <- -0.8
G[1,2] <- 0.5
sim_dat2 <- sim_mv(G = G,
                   N = 20000, J = 50000,
                   h2 = c(0.3, 0.4),
                   pi = 800/50000,
                   R_E = R_E, overlap_prop = 0)

res <- eb_mr(beta_hat_Y = sim_dat2$beta_hat[,2],
             se_Y = sim_dat2$se_beta_hat[,2],
             beta_hat_X = sim_dat2$beta_hat[,1],
             se_X = sim_dat2$se_beta_hat[,1],
             max_iter = 30,
             R = sim_dat2$R)




# Initialize mB and sB at 0 and v times the variance of the IVW estimator
library(TwoSampleMR)
library(dplyr)
library(mr.raps)
est_eff <- data.frame(bhat_x = sim_dat2$beta_hat[,1],
                      bhat_y = sim_dat2$beta_hat[,2],
                      se_bhat_x = sim_dat2$se_beta_hat[,1],
                      se_bhat_y = sim_dat2$se_beta_hat[,2]) %>%
  mutate(p_val_x = 2*pnorm(-abs(bhat_x/se_bhat_x)))

# Variants we select for single-variable MR
sv_mr_inst <- filter(est_eff, p_val_x < 5e-8)
mr_res1 <- mr_ivw(b_exp = sv_mr_inst$bhat_x, b_out = sv_mr_inst$bhat_y,
                  se_exp = sv_mr_inst$se_bhat_x, se_out = sv_mr_inst$se_bhat_y)

# > mr_res1
# $b
# [1] 0.4957664
#
# $se
# [1] 0.008556995

# Does not run
mr_res_rapsl2 <- mr.raps::mr.raps(b_exp = sv_mr_inst$bhat_x, b_out = sv_mr_inst$bhat_y,
                  se_exp = sv_mr_inst$se_bhat_x, se_out = sv_mr_inst$se_bhat_y,
                  over.dispersion = FALSE, loss.function = "l2")

mr_res_rapshuber <- mr.raps::mr.raps(b_exp = sv_mr_inst$bhat_x, b_out = sv_mr_inst$bhat_y,
                                  se_exp = sv_mr_inst$se_bhat_x, se_out = sv_mr_inst$se_bhat_y,
                                  over.dispersion = FALSE, loss.function = "huber")



}
