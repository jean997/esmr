library(dplyr)
library(TwoSampleMR)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)

out = args[1]
p_thresh = as.numeric(args[2])
nb_files = args[-c(1, 2)]




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


i <- ncol(beta_hat)
beta_hat <- data.frame(beta_hat)
se <- data.frame(se)
p <- data.frame(p)
pmin <- apply(p[,-1, drop = F], 1, min)
ix <- which(pmin < p_thresh)

exp <- as.matrix(beta_hat[ix, 2:i])
colnames(exp) <- nms[-1]
hdat <-  list(exposure_beta = exp,
              exposure_pval = as.matrix(p[ix, 2:i]),
              exposure_se = as.matrix(se[ix,2:i]),
              outcome_beta = beta_hat[ix,1],
              outcome_pval = p[ix,1],
              outcome_se = se[ix,1],
              expname = data.frame(id.exposure = nms, exposure = nms),
              outname = data.frame(id.outcome = nms[1], outcome = nms[1]))

t1 <- system.time(res1 <- mv_multiple(hdat, 
                      instrument_specific = FALSE));

t2 <- system.time(res2 <- mv_multiple(hdat, 
                      instrument_specific = TRUE));

res1$time <- t1
res2$time <- t2
res <- list(res1, res2)

saveRDS(res, file=out)

