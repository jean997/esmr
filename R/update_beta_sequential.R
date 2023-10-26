#'@export
update_beta_k <- function(j, k, dat){
  s_equal <- check_equal_omega(dat$omega)

  Va <- dat$l$a2bar - (dat$l$abar^2)



  if(s_equal){

    A <- t(dat$l$abar) %*% dat$l$abar + diag(colSums(Va))
    Astar <- dat$G %*% A %*% t(dat$G) # E[t(L) %*% L]

    s2inv <- (1/dat$sigma_beta^2) + Astar[k,k]*dat$omega[j,j]
    s2 <- 1/s2inv

    m_pt1 <- sum(dat$l$lbar[,k]*t(dat$omega[j,,drop =FALSE]%*%t(dat$Y)))
    m_pt2 <- sum(Astar[k, -k] * (dat$omega[j,] %*% dat$f$fbar[,-k]))
    m_pt3 <- Astar[k,k]*sum(dat$f$fbar[-j,k]*dat$omega[j,-j])

    m <- (m_pt1 - m_pt2 - m_pt3) * s2
  }else{
    stop("Not implemented yet.")

    ojj <- map(dat$omega, function(o){o[j,j]}) %>% unlist()

    s2inv <- (1/dat$sigma_beta^2) + sum(dat$l$l2bar[,k]*ojj)
    s2 <- 1/s2inv

    Oj <- map(omega, function(o){o[j,]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    Ojmj <- map(omega, function(o){o[j,-j]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    m_over_s2_pt1 <- sum(lbar[,k]*diag(Oj%*%t(R_k)))

    m_over_s2_pt2 <- sum(l2bar[,k]*t(t(fgbar[-j,k,drop=FALSE])%*%t(Ojmj)))

    m <- (m_over_s2_pt1 - m_over_s2_pt2) * s2
  }

  return(list(s = sqrt(s2), m = m))

}

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
    dat$beta$beta_s[i] <- b$S
    dat$f <- make_f(dat)

  }
  return(dat)
}
