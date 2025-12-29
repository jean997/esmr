remove_worst_factor <- function(dat){
  if(!dat$is_factors){
    stop("This function should only be called for a factors object.\n")
  }
  info_abar <- colSums(dat$l$abar^2)
  worst_abar <- which.min(info_abar[-1]) + 1
  remove_suggest <- which.max(abs(dat$G[,worst_abar]))
  if(remove_suggest == 2){ ## Don't remove primary exposure
    return(dat)
  }

  dat$factors_matrix <- dat$factors_matrix[,-(remove_suggest - 2), drop = FALSE]
  dat$k <- dat$k - 1
  dat$l$lbar <- dat$l$lbar[, -remove_suggest, drop = FALSE]
  dat$l$l2bar <- dat$l$l2bar[, -remove_suggest, drop = FALSE]
  dat$l$abar <- dat$l$abar[, -remove_suggest, drop = FALSE]
  dat$l$a2bar <- dat$l$a2bar[, -remove_suggest, drop = FALSE]
  dat$l$g_hat <- dat$l$g_hat[-remove_suggest, drop = FALSE]

  i <- which(dat$beta$beta_k == remove_suggest)
  dat$beta$beta_j <- dat$beta$beta_j[-i]
  dat$beta$beta_k <- dat$beta$beta_k[-i]
  dat$beta$beta_m <- dat$beta$beta_m[-i]
  dat$beta$beta_s <- dat$beta$beta_s[-i]
  dat$beta$V <- dat$beta$V[-i, -i, drop = FALSE]
  dat$beta$fix_beta <- dat$beta$fix_beta[-i]

  dat$beta$beta_k[dat$beta$beta_k > remove_suggest] <- dat$beta$beta_k[dat$beta$beta_k > remove_suggest] -1
  dat$f <- make_f_factors(dat)
  dat$G <- dat$G[-remove_suggest, -remove_suggest, drop = FALSE]

  return(dat)

}
