get_sigma <- function(R, S, s_equal, any_missing){
  p <- ncol(S)
  R_is_id <- is.null(R) | all(R == diag(p))

  if(s_equal & R_is_id){
    s <- S[1,]
    sigma <- diag(s^2)
  }else if(s_equal){
    s <- S[1,]
    sigma <- t(t(s * R) * s)
  }else if(R_is_id & !any_missing){
    sigma <- apply(S, 1, function(s){
      diag(s^2, nrow = p)
    }, simplify = FALSE)
  }else if(!any_missing){
    sigma <- apply(S, 1, function(s){
      t(t(s * R) * s)
    }, simplify = FALSE)
  }else{
    pat <- data.frame(is.na(S))
    pat$pat <- apply(pat, 1, function(x){paste0(x, collapse = "-")})
    pat_sum <- pat %>% group_by_all() %>% summarise(n = n())
    pat$which_pat <- match(pat$pat, pat_sum$pat)
    p <- ncol(S)
    z <- matrix(NA, nrow = p, ncol = p)
    sigma <- map(seq(nrow(pat_sum)), function(i){
      ixT <- which(pat_sum[i,] == TRUE)
      ixF <- which(pat_sum[i,] == FALSE)
      ii <- which(pat$which_pat == i)
      if(length(ixT) == 0){
        sig <- map(ii, function(j){
          s <- S[j,]
          t(t(s * R) * s)
        })
        return(sig)
      }
      myR <- R[-ixT,-ixT, drop = FALSE]
      cat(dim(myR), "\n")
      pn <- nrow(myR)
      sig <- map(ii, function(j){
        s <- S[j,-ixT]
        o <- t(t(s * myR) * s)
        om <- z
        om[ixF, ixF] <- o
        om
      })
      return(sig)
    }) %>% unlist(recursive = FALSE)
    ix1 <- sapply(seq(nrow(pat_sum)), function(i){which(pat$which_pat == i)}) %>% unlist()
    sigma <- sigma[match(seq(nrow(S)), ix1)]
  }
  return(sigma)
}




get_omega_tau <- function(sigma, tau, ld_scores, RE){
  sig_equal <- check_equal_omega(sigma)
  p <- length(ld_scores)


  if(sig_equal){
    omega <- lapply(1:p, function(i){
      solve(sigma + tau*ld_scores[i]*RE)
    })
  }else{
    omega <- lapply(1:p, function(i){
      solve(sigma[[i]] + tau*ld_scores[i]*RE)
    })
  }

  return(omega)
}
