# Selection-aware likelihood.
#
# Variants are selected when min_k p_jk < selection_pval (S_j = 1). Conditional on selection (Bayes' theorem)
#     l(Y_j | Theta, S_j = 1) = l(Y_j | Theta) - log P(S_j = 1 | Theta),
# so the objective for the n selected variants is  ELBO - n * log P(S = 1 | Theta).
#
# Model: Y_j = fgbar a_j + e_j with a_jk ~ (1 - pi_k) delta_0 + pi_k N(0, sigma_k^2) and e_j ~ N(0, S_j R S_j).
# On the z-score scale z_j = Y_j / S_j,  z | pattern a ~ N(0, D fgbar diag(a sigma^2) fgbar' D + R),  D = diag(1 / Sbar),
# where Sbar is the median standard error of each trait (the selection probability is evaluated at one representative S).
# P(S = 0 | Theta) = sum_a w_a P(|z_k| <= thr for all k | a), w_a = prod_k pi_k^a_k (1 - pi_k)^(1 - a_k).
#   - exact (pattern enumeration + mvtnorm box probabilities) when there are few factors and traits,
#   - otherwise Monte Carlo with fixed draws and importance weights in pi, so it is smooth in Theta.
# The prior update of every factor then maximizes  sum_i log m(x_i | g_k) - n log P(S = 1 | Theta)  over g_k = (pi_k, sigma_k).
# The direct effects (B) are still updated by the unpenalized closed-form step (plug-in).

selection_setup <- function(dat, selection_pval, mc_draws = 20000, mc_smooth = 0.05, seed = 1) {
    stopifnot(is.numeric(selection_pval), length(selection_pval) %in% c(1, dat$p), all(selection_pval > 0 & selection_pval < 1))
    thr <- qnorm(1 - rep_len(selection_pval, dat$p) / 2)
    Sbar <- apply(dat$S, 2, stats::median, na.rm = TRUE)
    omega1 <- if (dat$s_equal) dat$omega else dat$omega[[1]]
    Rn <- stats::cov2cor(solve(omega1))
    K <- dat$k
    exact <- (dat$p <= 3) && (K <= 8)
    sel <- list(thr = thr, Sbar = Sbar, Rn = Rn, exact = exact, mc_smooth = mc_smooth, K = K, p = dat$p)
    if (exact) {
        sel$patterns <- as.matrix(expand.grid(rep(list(0:1), K)))
    } else {
        set.seed(seed)
        sel$u <- matrix(stats::runif(mc_draws * K), mc_draws, K) # pattern draws (importance reference below)
        sel$xi <- matrix(stats::rnorm(mc_draws * K), mc_draws, K)
        sel$noise <- matrix(stats::rnorm(mc_draws * dat$p), mc_draws, dat$p) %*% chol(Rn)
        sel$pi_ref <- rep(0.1, K) # reference spike weights for the importance weights
    }
    sel
}

# pi_k and sigma_k of a point-normal normalmix fitted by ebnm
g_pi_sigma <- function(g) {
    nz <- g$sd > 0
    c(pi = if (any(nz)) sum(g$pi[nz]) else 0, sigma = if (any(nz)) sqrt(sum(g$pi[nz] * g$sd[nz]^2) / sum(g$pi[nz])) else 0)
}

# log P(S = 1 | Theta). `pis`, `sigmas`: length K; `fgbar`: p x K.
selection_log_p1 <- function(sel, fgbar, pis, sigmas) {
    D <- diag(1 / sel$Sbar, sel$p)
    DF <- D %*% fgbar # p x K
    if (sel$exact) {
        A <- sel$patterns
        p0 <- 0
        for (r in seq_len(nrow(A))) {
            a <- A[r, ]
            w <- prod(ifelse(a == 1, pis, 1 - pis))
            if (w < 1e-14) next
            C <- DF %*% diag(a * sigmas^2, sel$K) %*% t(DF) + sel$Rn
            C <- (C + t(C)) / 2
            bp <- if (sel$p == 1) {
                diff(stats::pnorm(c(-sel$thr, sel$thr), 0, sqrt(C[1, 1])))
            } else {
                mvtnorm::pmvnorm(lower = -sel$thr, upper = sel$thr, sigma = C)[1]
            }
            p0 <- p0 + w * bp
        }
    } else {
        pr <- sel$pi_ref
        pis <- pmin(pmax(pis, 1e-6), 1 - 1e-6)
        a <- sweep(sel$u, 2, pr, "<") * 1 # patterns drawn from the reference pi
        w <- exp(rowSums(a * rep(log(pis / pr), each = nrow(a)) + (1 - a) * rep(log((1 - pis) / (1 - pr)), each = nrow(a))))
        w <- w / nrow(a)
        sig <- a * sweep(sel$xi, 2, sigmas, "*")
        z <- sig %*% t(DF) + sel$noise
        pn <- stats::pnorm((matrix(sel$thr, nrow(z), sel$p, byrow = TRUE) - abs(z)) / sel$mc_smooth)
        inside <- exp(rowSums(log(pmax(pn, 1e-300))))
        p0 <- sum(w * inside)
    }
    log(max(1 - p0, 1e-300))
}

# Penalty (to subtract from the ELBO): -n log P(S = 1 | Theta); NULL if some prior is not available yet
selection_penalty <- function(dat) {
    if (!isTRUE(dat$selection_adjust)) return(0)
    g <- dat$l$g_hat
    if (length(g) < dat$k || any(vapply(seq_len(dat$k), function(k) is.null(g[[k]]), logical(1)))) return(0)
    ps <- t(vapply(seq_len(dat$k), function(k) g_pi_sigma(g[[k]]), numeric(2)))
    -dat$n * selection_log_p1(dat$sel, dat$f$fgbar, ps[, "pi"], ps[, "sigma"])
}

# Refit the point-normal prior of factor j to maximize  sum_i log m(x_i | g) - n log P(S = 1 | Theta).
# Returns a normalmix to use as g_init (with fix_g = TRUE) in ebnm.
selection_refit_g <- function(dat, j, x, s, g0) {
    g <- dat$l$g_hat
    if (length(g) < dat$k || any(vapply(seq_len(dat$k), function(k) is.null(g[[k]]), logical(1)))) return(NULL)
    ps <- t(vapply(seq_len(dat$k), function(k) g_pi_sigma(g[[k]]), numeric(2)))
    pis <- ps[, "pi"]; sigmas <- ps[, "sigma"]
    start <- g_pi_sigma(g0)
    start["pi"] <- min(max(start["pi"], 1e-4), 1 - 1e-4)
    start["sigma"] <- max(start["sigma"], 1e-6)
    sel <- dat$sel
    fgbar <- dat$f$fgbar
    n <- dat$n
    obj <- function(th) {
        pi_j <- stats::plogis(th[1]); sg_j <- exp(th[2])
        ll <- sum(log(pi_j * stats::dnorm(x, 0, sqrt(sg_j^2 + s^2)) + (1 - pi_j) * stats::dnorm(x, 0, s) + 1e-300))
        pj <- pis; sj <- sigmas; pj[j] <- pi_j; sj[j] <- sg_j
        -(ll - n * selection_log_p1(sel, fgbar, pj, sj))
    }
    # Multiple starts: the unadjusted fit can sit on a flat boundary (pi ~ 1) that the selection penalty makes a poor optimum
    sg0 <- start["sigma"]
    starts <- list(
        c(stats::qlogis(start["pi"]), log(sg0)),
        c(stats::qlogis(0.05), log(2 * sg0)),
        c(stats::qlogis(0.2), log(sg0)),
        c(stats::qlogis(0.5), log(sg0))
    )
    best <- NULL
    for (st in starts) {
        o <- tryCatch(
            stats::optim(unname(st), obj, method = "Nelder-Mead", control = list(maxit = 200, reltol = 1e-8)),
            error = function(e) NULL)
        if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    o <- best
    if (is.null(o)) return(NULL)
    g0$pi <- c(1 - stats::plogis(o$par[1]), stats::plogis(o$par[1]))
    g0$sd <- c(0, exp(o$par[2]))
    g0$mean <- c(0, 0)
    g0
}
