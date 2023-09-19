library(dplyr)
library(GRAPPLE)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
R <- readRDS(args[2])
p_thresh = as.numeric(args[3])
nb_files = args[-c(1:3)]




# Read in data
X <- map_dfr(nb_files, readRDS)

beta_hat <- X %>%
         select(ends_with(".beta"))

se <- X %>%
      select(ends_with(".se"))

p <- X %>%
      select(ends_with(".p"))
snps <- X$snp

nms <- names(X)[grep(".beta$", names(X))]
nms <- stringr::str_replace(nms, ".beta", "")

stopifnot(all(R$names %in% nms))
o <- match(R$names, nms)
beta_hat <- data.frame(beta_hat[, o])
se <- data.frame(se[, o])
Rcor <- cov2cor(R$R)

i <- ncol(beta_hat)
grapple_dat <- data.frame(cbind(beta_hat, se))
names(grapple_dat) <- c("gamma_out", paste0("gamma_exp", 1:(i-1)),
                        "se_out", paste0("se_exp", 1:(i-1)));
grapple_dat$selection_pvals <- apply(p[,-1, drop = F],1, min);

ix <- which(grapple_dat$selection_pvals < p_thresh)
if(length(ix) > 5000){
  ixn <- sample(ix, size = 5000, replace = FALSE)
  grapple_dat <- grapple_dat[ixn,]
}


t <- system.time(
        res <- grappleRobustEst(data = grapple_dat,
                                plot.it =FALSE,
                                p.thres = p_thresh,
                                cor.mat = Rcor,
                                niter = 100000));

res$time <- t
res$names <- nms

saveRDS(res, file=out)

