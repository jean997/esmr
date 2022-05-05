
#'@export
update_beta_k <- function(Y, k, lbar, l2bar, omega, beta_bar, sigma_beta = Inf){
  n <- nrow(Y)
  p <- ncol(Y)
  stopifnot(nrow(lbar) == n & ncol(lbar) == p)
  stopifnot(nrow(l2bar) == n & ncol(l2bar) == p)
  stopifnot(length(beta_bar) == p-1)
  stopifnot(k <= p-1)
  #First col of lbar is alpha the rest are gammas
  if("matrix" %in% class(omega)){
    s_equal <- TRUE
  }else{
    stopifnot(class(omega) != "list")
    stopifnot(length(omega) == n)
    s_equal <- FALSE
  }
  beta_bar <- matrix(beta_bar, nrow = p-1)
  if(s_equal){
    s2inv <- (1/sigma_beta^2) + sum(l2bar[,k+1])*omega[1,1]
    s2 <- 1/s2inv

    if(p ==2){
      m_over_s2_pt1 <- omega[1,1]*sum(lbar[,k+1]*(Y[,1] - lbar[,1]))
      m_over_s2_pt2 <- omega[1,k+1]*sum(lbar[,k+1]*Y[,k+1] - l2bar[,k+1])
    }else{
      m_over_s2_pt1 <- omega[1,1]*sum(lbar[,k+1]*(Y[,1] - lbar[,1]) -
                                        (lbar[,-c(1,k+1)] %*% beta_bar[,-k, drop = FALSE])*lbar[,k+1])
      m_over_s2_pt2 <- omega[1,k+1]*sum(lbar[,k+1]*Y[,k+1] - l2bar[,k+1]) +
                       sum(omega[1,-c(1, k+1)]*colSums(lbar[,(k+1), drop = FALSE] * ( Y[, -c(1, k +1)] - lbar[, -c(1, k + 1)] )))
    }
    m <- (m_over_s2_pt1 + m_over_s2_pt2) * s2
  }else{
    o11 <- map(omega, function(o){o[1,1]}) %>% unlist()
    o1k1 <- map(omega, function(o){o[1,k+1]}) %>% unlist()
    s2inv <- (1/sigma_beta^2) + sum(l2bar[,k+1]*o11)
    s2 <- 1/s2inv
    if(p == 2){
      m_over_s2_pt1 <- sum((lbar[,k+1]*(Y[,1] - lbar[,1]) )*o11)
      m_over_s2_pt2 <- sum((lbar[,k+1]*Y[,k+1] - l2bar[,k+1])*o1k1)
    }else{
      O1nK <- map(omega, function(o){o[1,-c(1, k+1)]}) %>% unlist() %>% matrix(nrow = n, byrow = TRUE)

      m_over_s2_pt1 <- sum(o11*(lbar[,k+1]*(Y[,1] - lbar[,1]) -
                                        (lbar[,-c(1,k+1)] %*% beta_bar[,-k, drop = FALSE])*lbar[,k+1]))
      m_over_s2_pt2 <- sum(o1k1*(lbar[,k+1]*Y[,k+1] - l2bar[,k+1])) +
        sum(  (lbar[,(k+1), drop = FALSE] * ( Y[, -c(1, k +1)] - lbar[, -c(1, k + 1)] )) * O1nK  )
    }
    m <- (m_over_s2_pt1 + m_over_s2_pt2)*s2
  }

  return(list(s = sqrt(s2), m = m))

}
