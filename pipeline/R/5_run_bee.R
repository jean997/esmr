library(dplyr)
library(MRBEE)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
R <- readRDS(args[2])
p_thresh = as.numeric(args[3])
pleio_p_thresh = as.numeric(args[4])
nb_files = args[-c(1:4)]


# Read in data
X <- map_dfr(nb_files, readRDS)

ntrait <- X %>%
          select(ends_with(".beta")) %>%
          ncol()

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
pmin <- apply(p[,-1, drop  = F], 1, min)
ix <- which(pmin < p_thresh)
bT <- list(R = Rcor, Ncor = Inf,
           EstHarm = beta_hat[ix,], 
           SEHarm =  se[ix,])
pD <- prepData(bT)

t1 <- system.time(fit <- MRBEE.IMRP(pD, PleioPThreshold = pleio_p_thresh))

fit$time <- t1
fit$names <- nms

saveRDS(fit, file=out)

