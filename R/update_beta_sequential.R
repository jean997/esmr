
update_beta_sequential <- function(dat){

  coords <- seq(length(dat$beta$beta_j))
  coords <- coords[!dat$beta$fix_beta]

  for(i in coords){
    k <- dat$beta$beta_k[i]
    j <- dat$beta$beta_j[i]

    #b <- update_beta_k(j,k,dat)
    b <- update_beta_joint(dat, j = j, ix = k)
    dat$beta$beta_m[i] <- b$m
    #dat$beta$beta_s[i] <- b$s
    dat$beta$beta_s[i] <- sqrt(b$S)
    dat$f <- make_f(dat)

  }
  return(dat)
}
