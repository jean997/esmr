# TODO: Want to add option to re-compute ELBO for different beta
ELBO <- function(x, new_fgbar = NULL) {
    # TODO: Maybe pass in direct effects instead?
    if (is.null(new_fgbar)) {
        new_fgbar <- x$f$fgbar
    }
    # If we want to restrict a specific edge we need to remove it from fgbar
    # Unclear how to update this exactly if we are removing a direct effect edge
    ll <- calc_ell2(x$Y, x$l$abar, x$l$a2bar, new_fgbar, x$omega, x$omega_logdet, x$s_equal)
    neg_kl_l <- x$l$kl
    # Note: Will need to do something like this instead:
    # dat$beta$kl <- - kl_mvn(
    #    dat$beta$beta_m[kl_ix], dat$beta$V[kl_ix, kl_ix,drop=F], 0, prior_cov_mat)
    neg_kl_beta <- x$beta$kl

    ll + neg_kl_l + neg_kl_beta
}
