update_beta_sequential_future <- function(dat, jj, fix_beta = NULL){
  p <- dat$p

  beta_j <- dat$beta$beta_j
  beta_k <- dat$beta$beta_k

  fbar_o <- dat$f$fbar_o
  f2bar_o <- dat$f$f2bar_o
  beta_m <- dat$beta$beta_m
  beta_s <- dat$beta$beta_s

  if(missing(jj)){
    coords <- seq(length(beta_j))
  }else{
    coords <- jj
  }
  if(!is.null(fix_beta)){
    if(length(fix_beta) == 1){
      fix_beta <- rep(fix_beta, length(coords))
    }else if(length(fix_beta) != length(coords)){
      stop("fix_beta and expected to have length 1 or ", length(coords), ". Found ", length(fix_beta), ".\n")
    }
    stopifnot(class(fix_beta) == "logical")
    coords <- coords[!fix_beta]
  }

  for(i in coords){
    k <- beta_k[i]
    j <- beta_j[i]

    R_k <- dat$Y - dat$l$lbar[,-k,drop=FALSE]%*%t(fbar_o[,-k,drop=FALSE])
    b <- update_beta_k(R_k = R_k, j=j, k=k,
                       lbar=dat$l$lbar, l2bar=dat$l$l2bar,
                       omega = dat$omega, fbar = fbar_o,
                       sigma_beta = dat$sigma_beta)
    beta_m[i] <- b$m
    beta_s[i] <- b$s
    f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
    fbar_o <- f$fbar_o
    f2bar_o <- f$f2bar_o
  }
  dat$beta$beta_m <- beta_m
  dat$beta$beta_s <- beta_s
  dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  return(dat)
}
