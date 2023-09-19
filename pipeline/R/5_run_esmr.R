library(dplyr)
library(purrr)
library(ebnm)
library(esmr)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
R <- readRDS(args[2])
ix1 = args[3]
nb_files = args[-c(1:3)]




# Read in data
X <- map_dfr(nb_files, readRDS)

ntrait <- X %>%
          dplyr::select(ends_with(".beta")) %>%
          ncol()

beta_hat <- X %>%
         dplyr::select(ends_with(".beta")) 

se <- X %>%
      dplyr::select(ends_with(".se")) 

snps <- X$snp

nms <- names(X)[grep(".beta$", names(X))]
nms <- stringr::str_replace(nms, ".beta", "")

stopifnot(all(R$names %in% nms))
o <- match(R$names, nms)
beta_hat <- data.frame(beta_hat[, o])
se <- data.frame(se[, o])
Rcor <- cov2cor(R$R)
p <- ncol(beta_hat)
if(ix1 == "NULL"){
    ix1 <- NULL
}
t <- system.time(
      fit <- esmr(beta_hat_Y <- beta_hat[,1], 
                  se_Y <- se[,1], 
                  beta_hat_X <- beta_hat[,2:p],
                  se_X <- se[, 2:p],
                  R = Rcor,
                  augment_G = TRUE,
                  g_type = "gfa",
                  ix1 = ix1,
                  ix0 = FALSE,
                  lfsr_thresh = 1))
fit$time <- t
fit$names <- nms
saveRDS(fit, file=out)

