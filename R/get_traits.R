
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
mr_pairs <- function(ex_ids, out_ids, inst_pval= 5e-8, method_list = c("mr_ivw")){
  ex_dat <- TwoSampleMR::extract_instruments(outcomes = ex_ids, p1 = inst_pval)
  out_dat <- TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = out_ids)
  dat <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  m <- mr(dat, method_list = method_list)
  return(m)
}

sloppy_cor <- function(id1, id2, inst_pval= 5e-8){
  ex_dat1 <- TwoSampleMR::extract_instruments(outcomes = id1)
  ex_dat1$id <- str_split(ex_dat1$exposure, "id:") %>% map(., 2) %>% unlist()
  ex_dat2 <- TwoSampleMR::extract_instruments(outcomes = id2)
  ex_dat2$id <- str_split(ex_dat2$exposure, "id:") %>% map(., 2) %>% unlist()
  cor_vals <- expand.grid(id1 = id1, id2 = id2, stringsAsFactors = FALSE)
  c <- map2(cor_vals$id1, cor_vals$ids, function(x, y){

  })
  return(m)
}


filter_traits <- function(id_x, id_y, candidate_traits, cor_thresh = 0.9,
                          cor_func = sloppy_cor){

}

sloppy_cor <- function(idlist, r2, kb, pop, access_token){
  top_hits_all <- purrr::map(idlist, function(id){
                      tophits(id = id, pval = pval,
                        r2 = r2, kb = kb, pop = pop,
                        access_token =  access_token)
  })
}
