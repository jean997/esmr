library(dplyr)
library(purrr)
library(readr)
library(bigsnpr)
library(stringr)


args <- commandArgs(trailingOnly=TRUE)
l2_dir <- args[1]
gwas_info <- read_csv(args[2])
out <- args[3]
inp_files <- args[-c(1:3)]


ld <- purrr::map_dfr(1:22, function(c){
  read_table(paste0(l2_dir, c, ".l2.ldscore.gz"))
})

M <- purrr:::map(1:22, function(c){
  read_lines(paste0(l2_dir, c, ".l2.M_5_50"))
}) %>% unlist() %>% as.numeric() %>% sum()


dat <- map_dfr(inp_files, function(x){readRDS(x)})
dat <- filter(dat, snp %in% ld$SNP)

Z <- dat %>% dplyr::select(ends_with(".beta"))
n <- names(Z)
s <- dat %>% dplyr::select(ends_with(".se"))
Z <- Z/s
n <- names(Z)
n <- stringr::str_replace(n, ".beta", "")
names(Z) <- n
Z$SNP <- dat$snp

full_dat <- inner_join(Z, ld) 

h2 <- lapply(n, function(nn){
                 ss <- gwas_info$pub_sample_size[gwas_info$name == nn]
                 i <- which(!is.na(full_dat[[nn]]))
                 snp_ldsc(ld_score = full_dat$L2[i],
                    ld_size = M,
                    chi2 = (full_dat[[nn]][i])^2,
                    sample_size = ss,
                    blocks = NULL)
})


name_pairs <- expand.grid(t1 = seq_along(n), t2 = seq_along(n)) %>% 
                filter(t1 < t2) 
name_pairs$t1 <- n[name_pairs$t1]
name_pairs$t2 <- n[name_pairs$t2]

np_res <- lapply(seq(nrow(name_pairs)), function(ii){
                     cat(ii, " ")
                     nn1 <- name_pairs$t1[ii]
                     nn2 <- name_pairs$t2[ii]
                     ss1 <- gwas_info$pub_sample_size[gwas_info$name == nn1]
                     ss2<- gwas_info$pub_sample_size[gwas_info$name == nn2]
                     i <- which(!is.na(full_dat[[nn1]]) & !is.na(full_dat[[nn2]]))
                    snp_ldsc_rg(ld_score = full_dat$L2[i],
                                ld_size = M,
                                sample_size_1 = ss1,
                                sample_size_2 = ss2,
                                z1 = full_dat[[nn1]][i],
                                z2 = full_dat[[nn2]][i],
                                blocks = NULL, h2_se = FALSE)
})
name_pairs$cov <- np_res %>%map("int") %>% unlist()
x <- name_pairs
names(x) <- c("t2", "t1", "cov") 
name_pairs <- bind_rows(name_pairs, data.frame(t1 = n, t2 = n, cov = h2 %>% map("int") %>% unlist()))
name_pairs <- bind_rows(name_pairs, x)
cov_mat <- reshape2::dcast(name_pairs, t1 ~ t2)

nms <- paste0(as.vector(cov_mat$t1))
R <- as.matrix(cov_mat[,-1])
R <- Matrix::nearPD(R,  posd.tol = 1e-3) %>% with(., as.matrix(mat));

o <- match(gwas_info$name, nms)
R <- R[o, o]
eS <- eigen(R)

ret <- list(R = R, names = nms[o], eS = eS)

saveRDS(ret, file=out)
             

