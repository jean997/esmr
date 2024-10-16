
get_top_combinations <- function(x, max_logsumexp, nmax = Inf){
  c <- apply(x, 1, which.max)

  combination_list <- list(c)
  v <- sapply(1:nrow(x), function(i){x[i,c[i]]}) %>% sum()
  values <- c(v)
  l <- Inf
  i <- 1
  done <- FALSE
  #check_ix <- c(1)
  while(!done){
    #i <- check_ix[1]
    #check_ix <- check_ix[-1]
    #cat(i, " ")
    v <- values[i]
    c <- combination_list[[i]]
    new_comb <- best_valid_one_move(c, x, combination_list)
    losses <- new_comb[[2]]
    #cat(losses, "\n")
    if(all(is.na(losses))){
      if(i == length(values)) done <- TRUE
    }else{
      new_combs <- new_comb[[1]][!is.na(losses)]
      losses <- losses[!is.na(losses)]
      combination_list <- c(combination_list, new_combs)
      new_vals <- v + losses
      values <- c(values, new_vals)
      o <- order(values, decreasing = T)
      values <- values[o]
      combination_list <- combination_list[o]
    }
    i <- i + 1
    w <- matrixStats::logSumExp(values[1:i])
    if(i > nmax | w > max_logsumexp){
      done <- TRUE
    }
  }
  top_combinations <- combination_list[1:i] %>% unlist() %>%
                         matrix(ncol = nrow(x), byrow = T)
  return(list(combs = top_combinations, values = values[1:i]))
}

best_valid_one_move <- function(c, x, combination_list){
  n <- length(c)
  best_one_move <- lapply(1:n, function(i){
    cix <- lapply(combination_list, function(cc){
      if(all(cc[-i] == c[-i])){
        return(cc[i])
      }else{
        return(NULL)
      }
    }) %>% unlist()
    cix <- c(c[i], cix)

    xx <- x[i,]
    xx[cix] <- -Inf
    l <- xx-x[i, c[i]]

    if(!any(is.finite(l))){
      return(c(NA, NA))
    }
    newc <- c
    newc[i] <- which.max(l)
    return(list(newc, max(l)))
  })
  combinations <- map(best_one_move, 1)
  losses <- map(best_one_move, 2) %>% unlist()
  return(list(combinations, losses))
}

