#'@export
update_beta_k <- function(j, k, dat){
  s_equal <- check_equal_omega(dat$omega)

  Vl <- dat$l$l2bar_o - (dat$l$lbar_o^2)

  if(s_equal){
    A <- t(dat$l$lbar_o) %*% dat$l$lbar_o + diag(colSums(Vl))
    Astar <- dat$G %*% A %*% t(dat$G)

    s2inv <- (1/dat$sigma_beta^2) + Astar[k,k]*dat$omega[j,j]
    s2 <- 1/s2inv

    m_pt1 <- sum(dat$l$lbar[,k]*t(dat$omega[j,,drop =FALSE]%*%t(dat$Y)))
    m_pt2 <- sum(Astar[k, -k] * (dat$omega[j,] %*% dat$f$fbar_o[,-k]))
    m_pt3 <- Astar[k,k]*sum(dat$f$fbar_o[-j,k]*dat$omega[j,-j])

    m <- (m_pt1 - m_pt2 - m_pt3) * s2
  }else{
    stop("Not implemented yet.")
    ojj <- map(dat$omega, function(o){o[j,j]}) %>% unlist()

    s2inv <- (1/dat$sigma_beta^2) + sum(dat$l$l2bar[,k]*ojj)
    s2 <- 1/s2inv

    Oj <- map(omega, function(o){o[j,]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    Ojmj <- map(omega, function(o){o[j,-j]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    m_over_s2_pt1 <- sum(lbar[,k]*diag(Oj%*%t(R_k)))

    m_over_s2_pt2 <- sum(l2bar[,k]*t(t(fbar[-j,k,drop=FALSE])%*%t(Ojmj)))

    m <- (m_over_s2_pt1 - m_over_s2_pt2) * s2
  }

  return(list(s = sqrt(s2), m = m))

}

update_beta_sequential <- function(dat){
  #p <- dat$p

  #beta_j <- dat$beta$beta_j
  #beta_k <- dat$beta$beta_k

  # fbar_o <- dat$f$fbar_o
  # f2bar_o <- dat$f$f2bar_o
  # beta_m <- dat$beta$beta_m
  # beta_s <- dat$beta$beta_s

  coords <- seq(length(dat$beta$beta_j))
  coords <- coords[!dat$beta$fix_beta]

  for(i in coords){
    k <- dat$beta$beta_k[i]
    j <- dat$beta$beta_j[i]

    b <- update_beta_k_future(j,k,dat)
    dat$beta$beta_m[i] <- b$m
    dat$beta$beta_s[i] <- b$s
    dat$f <- make_f(dat)
    # fbar_o <- f$fbar_o
    # f2bar_o <- f$f2bar_o
  }
  # dat$beta$beta_m <- beta_m
  # dat$beta$beta_s <- beta_s
  #dat$f <- dat$f_fun(dat$beta$beta_m, dat$beta$beta_s)
  return(dat)
}
