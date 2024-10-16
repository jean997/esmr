update_tau <- function(dat, tau_min = 0, tau_max = 1){
  s <- 1/(tau_max - tau_min)
  new_tau <- optimize(ell_tau, interval = c(tau_min*s, tau_max*s), scale = s, dat = dat, maximum = TRUE)
  dat$tau <- new_tau$maximum/s
  dat$omega <- get_omega_tau(dat$sigma, dat$tau, dat$ld_scores, dat$RE)
  dat$f <- make_f(dat)
  return(dat)
}

ell_tau <- function(tau, scale, dat){
  tau <- tau/scale
  dat$tau <- tau
  dat$omega <- get_omega_tau(dat$sigma, tau, dat$ld_scores, dat$RE)
  dat$f <- make_f(dat)
  #log_det <- sapply(dat$omega, function(o){log(det(o))}) %>% sum()
  with(dat, calc_ell2(Y, l$abar, l$a2bar, f$fgbar,  omega,  s_equal)) #+ 0.5*log_det

}

