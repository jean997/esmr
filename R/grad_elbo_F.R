calc_grad_hess_B <- function(
    Y,
    abar,
    a2bar,
    bgbar,
    omega,
    s_equal = TRUE
) {

  n <- nrow(Y)
  p <- ncol(Y)
  k <- ncol(abar)

  stopifnot(ncol(bgbar) == k)
  stopifnot(nrow(bgbar) == p)

  varabar <- a2bar - abar^2

  grad <- matrix(0, nrow = p, ncol = k)
  hdiag <- matrix(0, nrow = p, ncol = k)

  # ==========================================================
  # Equal omega case
  # ==========================================================

  if (s_equal) {

    OB <- omega %*% bgbar   # p x k
    YO <- Y %*% omega       # n x p

    for (j in seq_len(k)) {

      mj <- abar[, j]       # length n
      vj <- varabar[, j]    # length n

      # ------------------------------------------------------
      # term1
      #
      # sum_r y_r^T Omega_.i * m_rj
      #
      # -> p-vector
      # ------------------------------------------------------

      term1 <- colSums(YO * mj)

      # ------------------------------------------------------
      # term2
      #
      # variance correction
      #
      # v_j * (Omega B_.j)
      #
      # -> p-vector
      # ------------------------------------------------------

      term2 <- sum(vj) * OB[, j]

      # ------------------------------------------------------
      # term3
      #
      # mean interaction
      #
      # m_j * (m^T B^T Omega_.i)
      #
      # -> scalar × p-vector
      # ------------------------------------------------------

      proj_j <- as.vector(abar %*% OB[, j])

      scalar_part <- sum(mj * proj_j)

      term3 <- rep(scalar_part, p)

      # ------------------------------------------------------

      grad[, j] <- term1 - term2 - term3

      # diagonal Hessian

      hdiag[, j] <- -sum(a2bar[, j]) * diag(omega)
    }

  } else {

    # ========================================================
    # Unequal omega case
    # ========================================================

    for (r in seq_len(n)) {

      yi <- Y[r, ]              # length p
      li <- abar[r, ]           # length k
      li2 <- a2bar[r, ]         # length k
      vi <- li2 - li^2

      Oi <- omega[[r]]          # p x p
      OB <- Oi %*% bgbar        # p x k

      for (j in seq_len(k)) {

        # ----------------------------------------------------
        # term1
        #
        # li[j] * (yi^T Oi)
        #
        # -> p-vector
        # ----------------------------------------------------

        term1 <- li[j] * as.vector(yi %*% Oi)

        # ----------------------------------------------------

        term2 <- vi[j] * OB[, j]

        # ----------------------------------------------------
        # scalar correction
        # ----------------------------------------------------

        scalar_part <- li[j] * sum(li * crossprod(bgbar, Oi %*% diag(p))[, 1])

        term3 <- rep(scalar_part, p)

        # ----------------------------------------------------

        grad[, j] <- grad[, j] + term1 - term2 - term3

        hdiag[, j] <- hdiag[, j] -
          li2[j] * diag(Oi)
      }
    }
  }

  # ==========================================================
  # Pairwise Hessian helper
  # ==========================================================

  hessian_pair <- function(i, j, i2, j2) {

    if (s_equal) {

      val <- -sum(
        abar[, j] * abar[, j2]
      ) * omega[i, i2]

    } else {

      val <- 0

      for (r in seq_len(n)) {
        val <- val -
          abar[r, j] *
          abar[r, j2] *
          omega[[r]][i, i2]
      }
    }

    val
  }

  list(
    grad = grad,
    hdiag = hdiag,
    hessian_pair = hessian_pair
  )
}

calc_grad_upper <- function(
    bgbar,
    lbar,
    l2bar,
    Y,
    omega,
    s_equal = TRUE
) {

  n <- nrow(Y)
  p <- ncol(Y)

  grad_upper <- matrix(0, nrow = p, ncol = p)

  # For each forbidden edge (i,j) where i < j in ordering coordinates,
  # compute the gradient of ELBO w.r.t. B[i,j]
  #
  # The gradient measures signal for swapping positions i and j
  # Based on how much E[L_i] correlates with unexplained variance in Y_j
  # weighted by precision (Omega)

  # Compute residuals: Y - E[L] B^T
  residuals <- Y - lbar %*% t(bgbar)

  if (s_equal) {
    # For constant precision matrix across samples:
    # grad[i,j] ∝ E[L_i]^T Omega residuals[j]

    OmegaT <- omega  # p x p precision matrix

    # For each pair (i,j) with i < j (forbidden edge in current ordering)
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {

        # Precision-weighted covariance: how much L_i explains residuals_j
        # when properly accounting for the marginal precision structure
        weighted_residuals <- as.vector(OmegaT %*% residuals[, j])
        cov_term <- sum(lbar[, i] * weighted_residuals)

        grad_upper[i, j] <- cov_term
      }
    }

  } else {

    # For varying precision matrices (unequal omega):
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {

        cov_term <- 0
        for (r in 1:n) {
          # Weight residual by precision structure
          weighted_resid <- as.vector(omega[[r]] %*% matrix(residuals[r, j], nrow = p))
          cov_term <- cov_term + lbar[r, i] * weighted_resid[i]
        }

        grad_upper[i, j] <- cov_term / n
      }
    }
  }

  grad_upper
}

nesmr_with_grad <- function(
    dat,
    order = seq_len(ncol(dat$Y))
) {

  # ----------------------------------------------------------
  # Basic dimensions
  # ----------------------------------------------------------

  stopifnot(length(order) == ncol(dat$Y))

  Y <- dat$Y[, order, drop = FALSE]

  abar <- dat$l$abar
  a2bar <- dat$l$a2bar

  fgbar <- dat$f$fgbar[order, , drop = FALSE]

  omega <- dat$omega
  omega_logdet <- dat$omega_logdet
  s_equal <- check_equal_omega(omega)

  # ----------------------------------------------------------
  # Expected log likelihood
  # ----------------------------------------------------------

  ll <- calc_ell2(
    Y = Y,
    abar = abar,
    a2bar = a2bar,
    fgbar = fgbar,
    omega = omega,
    omega_logdet = omega_logdet,
    s_equal = s_equal
  )

  # ----------------------------------------------------------
  # Full ELBO
  #
  # matches:
  #
  # obj = ll + dat$l$kl + dat$beta$kl
  # ----------------------------------------------------------

  kl_l <- dat$l$kl

  if (is.null(dat$beta$kl)) {
    kl_beta <- 0
  } else {
    kl_beta <- dat$beta$kl
  }

  elbo <- ll + kl_l + kl_beta

  # ----------------------------------------------------------
  # Gradient + Hessian wrt B (direct effect matrix)
  # ----------------------------------------------------------

  gh <- calc_grad_hess_B(
    Y = Y,
    abar = abar,
    a2bar = a2bar,
    bgbar = fgbar,
    omega = omega,
    s_equal = s_equal
  )

  # ----------------------------------------------------------
  # Forbidden-edge gradients for topological swaps
  #
  # grad_upper[i,j] = gradient of ELBO w.r.t. B[i,j]
  # for i < j in the current ordering (forbidden edges)
  #
  # These measure how much allowing each forbidden edge would
  # improve the ELBO - used to guide topological search
  # ----------------------------------------------------------

  grad_upper <- calc_grad_upper(
    bgbar = fgbar,
    lbar = abar,
    l2bar = a2bar,
    Y = Y,
    omega = omega,
    s_equal = s_equal
  )

  # ----------------------------------------------------------
  # Optional: current lower-triangular beta matrix
  #
  # useful for debugging
  # ----------------------------------------------------------

  beta_mat <- NULL

  if (!is.null(dat$beta$beta_m)) {

    p <- ncol(dat$Y)

    beta_mat <- matrix(0, p, p)

    #
    # beta_j = column
    # beta_k = row
    #

    idx <- cbind(
      dat$beta$beta_k,
      dat$beta$beta_j
    )

    beta_mat[idx] <- dat$beta$beta_m

    beta_mat <- beta_mat[order, order, drop = FALSE]
  }

  # ----------------------------------------------------------
  # Return object
  # ----------------------------------------------------------

  out <- list(
    order = order,

    # objective pieces
    ll = ll,
    kl_l = kl_l,
    kl_beta = kl_beta,
    elbo = elbo,

    # derivatives
    grad = gh$grad,
    hdiag = gh$hdiag,
    grad_upper = grad_upper,
    hessian_pair = gh$hessian_pair,

    # optional diagnostics
    beta_mat = beta_mat
  )

  class(out) <- "nesmr_grad"

  return(out)
}