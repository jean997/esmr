
#'@export
retrieve_traits <- function(id_x, pval_x, pval_z, pop, batch,
                            r2 = 0.001, kb = 10000,
                            access_token = check_access_token(),
                            min_snps = 3){
  top_hits <- tophits(id = id_x, pval = pval_x,
                      r2 = r2, kb = kb, pop = pop,
                      access_token =  access_token)
  cat("Retrieved", nrow(top_hits), "instruments for", id_x, "\n")
  phe <- phewas(variants = top_hits$rsid, pval = pval_z,
                batch = batch, access_token = access_token)
  cat("Retrieved", nrow(phe), "associations with", length(unique(phe$id)), "traits", "\n")
  x <- phe %>% group_by(id) %>% summarize(n = length(unique(rsid)))
  cat(sum(x$n >= min_snps), "traits have at least", min_snps, "shared variants with", id_x, "\n")
  ids <- x$id[x$n >= min_snps]
  phe <- filter(phe, id %in% ids)
  ret <- list("topx" = top_hits, "phe" = phe)
  return(ret)
}


#'@export
mr_pairs <- function(ids1, ids2, inst_pval= 5e-8, method_list = c("mr_ivw")){
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids1, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids2)
  dat_1_2 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m_1_2 <- mr(dat_1_2, method_list = method_list) # Replace for grapple or esmr

  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ids2, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids1)
  dat_2_1 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m_2_1 <- mr(dat_2_1, method_list = method_list)


  cor_vals <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = FALSE) %>%
              filter(id1 != id2)
  cor_vals$cor <- map2(cor_vals$id1, cor_vals$id2, function(x, y){
    X_1_2 <- filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
             rename(beta1 = beta.exposure, beta2 = beta.outcome) %>%
             select(SNP, beta1, beta2)
    X_2_1 <- filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
      rename(beta2 = beta.exposure, beta1 = beta.outcome) %>%
      select(SNP, beta1, beta2) %>%
      filter(!SNP %in% X_1_2$SNP)
    X <- bind_rows(X_1_2, X_2_1)
    with(X, cor(beta1, beta2))
  }) %>% unlist()

  return(list(mr12 = m_1_2, mr21 = m_2_1, cor = cor_vals))
}

