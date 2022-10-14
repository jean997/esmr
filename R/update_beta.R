
#'@export
update_beta_k <- function(R_k, j, k, lbar, l2bar, omega, fbar, sigma_beta = Inf){
  n <- nrow(R_k)
  p <- ncol(R_k)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  stopifnot(nrow(fbar) == p & ncol(fbar) == p)
  stopifnot(j <= p & k <= p)
  #First col of lbar is alpha the rest are gammas
  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    cat(class(omega), "\n")
    stopifnot(class(omega) == "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }

  if(s_equal){
    s2inv <- (1/sigma_beta^2) + sum(l2bar[,k])*omega[j,j]
    s2 <- 1/s2inv


    m_over_s2_pt1 <- sum(lbar[,k]*t(omega[j,,drop =FALSE]%*%t(R_k)))
    m_over_s2_pt2 <- sum(l2bar[,k]*(sum(fbar[-j,k]*omega[j,-j])))

    m <- (m_over_s2_pt1 - m_over_s2_pt2) * s2
  }else{
    ojj <- map(omega, function(o){o[j,j]}) %>% unlist()

    s2inv <- (1/sigma_beta^2) + sum(l2bar[,k]*ojj)
    s2 <- 1/s2inv

    Oj <- map(omega, function(o){o[j,]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    Ojmj <- map(omega, function(o){o[j,-j]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)
    m_over_s2_pt1 <- sum(lbar[,k]*diag(Oj%*%t(R_k)))

    m_over_s2_pt2 <- sum(l2bar[,k]*t(t(fbar[-j,k,drop=FALSE])%*%t(Ojmj)))

    m <- (m_over_s2_pt1 - m_over_s2_pt2) * s2
  }

  return(list(s = sqrt(s2), m = m))

}
