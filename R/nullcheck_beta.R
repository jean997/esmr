nullcheck_beta <- function(dat, k, max_iter, tol = default_precision(dim(dat$Y))){

  dat_k <- zero_beta(dat, k)
  dat_k <- ebmr_solve(dat_k, max_iter, tol)
  obj1 <- dat$obj[length(dat$obj)]
  obj2 <- dat_k$obj[length(dat_k$obj)]
  obj_diff <- obj1 - obj2
  if(obj_diff < tol) return(dat_k)
  return(dat)
}

zero_beta <- function(dat, k){
  dat$beta$beta_j <- dat$beta$beta_j[-k]
  dat$beta$beta_k <- dat$beta$beta_k[-k]
  dat$beta$beta_m <- dat$beta$beta_m[-k]
  dat$beta$beta_s <- dat$beta$beta_s[-k]
  dat$f_fun <- make_f_fun(dat$p, dat$beta$beta_j, dat$beta$beta_k)
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  if(length(dat$beta_j) == 0) dat$fix_beta = TRUE
  return(dat)
}
