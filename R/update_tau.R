update_tau <- function(dat){
  R <- with(dat, Y - l$lbar %*% t(f$fbar))
  tau <- dat$n/colSums(R^2)

  if("matrix" %in% class(dat$omega)){
    s_equal <- TRUE
    d <- pmax(diag(dat$omega_given), tau)
    diag(dat$omega) <- d
  }else{
    s_equal <- FALSE
    dat$omega <- map(seq_along(dat$omega), function(i){
      d <- pmax(diag(dat$omega_given[[i]]), tau)
      diag(dat$omega[[i]]) <- d
    })
  }
  return(dat)
}
