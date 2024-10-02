check_numeric <- function(x, n){
  if(is.null(x)) return(x)
  if("matrix" %in% class(x) | "data.frame" %in% class(x)){
    if(ncol(x) > 1) stop(paste0(deparse(substitute(x)), " must be a numeric vector or one column array."))
    x <- as.numeric(x[,1])
  }else if(!"numeric" %in% class(x)){
    stop(paste0(deparse(substitute(x)), " must be a numeric vector or one column array."))
  }
  if(!missing(n)){
    if(length(x) != n) stop(paste0("Expected ", deparse(substitute(x)), " to have length ", n, ", found ", length(x), "\n"))
  }
  return(x)
}

check_matrix <- function(x,  n, p){
  if(is.null(x)) return(x)
  if("data.frame" %in% class(x)){
    cat("Coercing ", deparse(substitute(x)), " to matrix.\n")
    x <- as.matrix(x)
  }else if("numeric" %in% class(x)){
    x <- as.matrix(x, ncol = 1)
  }else if(!"matrix" %in% class(x)){
    stop(paste0(deparse(substitute(x)), " must be a numeric vector, matrix, or data.frame."))
  }
  if(!missing(n)){
    if(nrow(x) != n) stop(paste0("Expected ", deparse(substitute(x)), " to have ", n, " rows, found ", nrow(x), "\n"))
  }
  if(!missing(p)){
    if(ncol(x) != p) stop(paste0("Expected ", deparse(substitute(x)), " to have ", p, " columns, found ", ncol(x), "\n"))
  }
  return(x)
}

check_R <- function(R, tol = 1e-8){
  if(is.null(R)) return(R)
  if(!Matrix::isSymmetric(R)){
    stop("R is not symmetric.\n")
  }
  evR <- eigen(R, only.values = TRUE)$values
  if(!all(evR > tol)){
    stop("R is not positive definite.\n")
  }
  #if(!all(diag(R) == 1)){
  if(!all.equal(diag(R), rep(1, nrow(R), use.names = F) )){
    warning("R is not a correlation matrix.\n")
  }
  return(R)
}

check_missing <- function(Y, S){
  missing_ix <- which(is.na(Y))
  any_missing <- length(missing_ix) > 0
  p <- ncol(Y)
  n <- nrow(Y)
  if(!any_missing){
    s_equal <- apply(S, 2, function(x){all(x == x[1])}) |> all()
    return(list(Y = Y, S = S, n = n, p = p,
                s_equal = s_equal,  any_missing = FALSE))
  }

  if(any(is.na(S[-missing_ix]))) stop("Found missing SEs for non-missing effect estimates.\n")

  nmiss_r <- rowSums(is.na(Y))
  nmiss_c <- colSums(is.na(Y))
  if(any(nmiss_r == p)) stop("Data cannot have variants missing for all traits.\n")
  if(any(nmiss_c == n)) stop("At least one trait is missing for all variants.\n")
  Y[missing_ix] <- 0
  S[missing_ix] <- NA
  return(list(Y = Y, S = S, n = n, p = p,
              s_equal = FALSE, any_missing = TRUE))
}


check_equal_omega <- function(omega){
  if("matrix" %in% class(omega)){
    #check_matrix(omega, "omega", p, p)
    return(TRUE)
  }else{
    stopifnot(class(omega) == "list")
    return(FALSE)
  }
}

check_B_template <- function(B, p){
  B <- check_matrix(B, p, p)

  if(!all(B %in% c(0, 1))) stop("Direct effects template should contain only 0 and 1 entries.")
  #if(!all(diag(B) ==0)){
  if(!all.equal(diag(B), rep(0, p))){
    stop("Direct effects template must have 0s on the diagonal.")
  }
  B_tot <- tryCatch(direct_to_total(B), error = function(e){
    stop("Check that supplied template corresponds to a valid DAG.\n")
  })
  #if(!all(diag(B_tot) == 0)){
  if(!all.equal(diag(B_tot), rep(0, p))){
    stop("Check that supplied template corresponds to a valid DAG.\n")
  }

  B_tot[!B_tot == 0] <- 1
  which_tot_u <- which(B_tot != 0 & B != 0, arr.ind = TRUE)
  which_tot_c <- which(B_tot != 0 & B == 0, arr.ind = TRUE)
  return(list(B_dir = B, B_tot = B_tot, which_tot_u = which_tot_u, which_tot_c = which_tot_c ))
}

