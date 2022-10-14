
#'@export
retrieve_traits <- function(id_x, pval_x, pval_z, pop,
                            batch,
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

my_steiger_filtering <- function (dat){
  plyr::ddply(dat, c("id.exposure", "id.outcome"), function(d){
    dd <- try(TwoSampleMR:::steiger_filtering_internal(d), silent = TRUE)
    if(class(dd) == "try-error"){
      dd <- d
      dd$rsq.exposure <- NA
      dd$effective_n.exposure <- NA
      dd$rsq.outcome <- NA
      dd$effective_n.outcome <- NA
      dd$steiger_dir <- NA
      dd$steiger_pval <- NA
    }
    return(dd)
  })
}

#'@export
mr_one_d <- function(ids1, ids2, inst_pval= 4.9e-8, method_list = c("mr_ivw"),
                     steiger_filter = TRUE, Qpval_thresh = 0.05, ntries = 3,
                     groups = 10){
  anal <- expand.grid(st =steiger_filter, q = Qpval_thresh, m = method_list, stringsAsFactors = F)

  if(length(ids1) > groups){
    start <- seq(1, length(ids1), by = groups)
    stop <- pmin(start + groups -1, length(ids1))
  }else{
    start <- 1
    stop <- length(ids1)
  }

  ex_dat <- map_dfr(seq_along(start), function(i){
              d <- try(TwoSampleMR::extract_instruments(outcomes = ids1[start[i]:stop[i]], p1 = inst_pval), silent = TRUE)
              cat("done ", i, "\n")
              return(d)
  })

  out_dat <- try(TwoSampleMR::extract_outcome_data(snps = ex_dat$SNP, outcomes = ids2), silent = TRUE)

  dat_1_2 <- TwoSampleMR::harmonise_data(ex_dat, out_dat)
  dat_1_2 <- add_metadata(dat_1_2)
  st <- my_steiger_filtering(dat_1_2)
  dat_1_2 <- st

  ids <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = F)
  m1 <- map2_dfr(ids$id1, ids$id2, function(x, y){
    d <- filter(dat_1_2, id.exposure == x & id.outcome == y)
    m <-  try(TwoSampleMR:::mr_mean_ivw(d), silent = TRUE)
    if(!class(m) == "try-error"){
      d <- left_join(d,m$outliers )
    }else if(nrow(d) == 0){
      d$Qpval <- NULL
    }else{
      cat(x, " ", y, "\n")
      d$Qpval <- NA
    }
    return(d)
  })
  dat_1_2 <- m1

  m_1_2 <- map_dfr(seq(nrow(anal)), function(i){
    d <- dat_1_2
    if(anal$st[i]){
      d <- filter(dat_1_2, steiger_dir == TRUE)
    }
    if(anal$q[i] > 0){
      d <- filter(d, Qpval > anal$q[i] & !is.na(Qpval) )
    }
    r <- try(TwoSampleMR::mr(d, method_list = anal$m[i]) , silent = TRUE)# Replace for grapple or esmr
    if(!class(r) == "try-error"){
      r$steiger_filtering <- anal$st[i]
      r$Qpval <- anal$q[i]
    }else{
      cat("error: ", i, "\n")
      r <- data.frame(id.exposure = NA, id.outcome = NA, outcome = NA, exposure = NA, method = anal$m[i],
                      snp = NA, b = NA, se = NA, pval = NA, steiger_filtering = anal$st[i], Qpval = anal$q[i])
    }
    return(r)
  })
  return(list(m = m_1_2, dat = dat_1_2))
}

#'@export
mr_pairs <- function(ids1, ids2, inst_pval= 4.9e-8, method_list = c("mr_ivw"), steiger_filter = TRUE, Qpval_thresh = 0.05){

  m_2_1 <- mr_one_d(ids2, ids1,inst_pval, method_list, steiger_filter, Qpval_thresh)
  dat_2_1 <- m_2_1$dat
  m_2_1 <- m_2_1$m


  m_1_2 <- mr_one_d(ids1, ids2,inst_pval, method_list, steiger_filter, Qpval_thresh)
  dat_1_2 <- m_1_2$dat
  m_1_2 <- m_1_2$m

  cor_vals <- expand.grid(id1 = ids1, id2 = ids2, stringsAsFactors = FALSE) %>%
              filter(id1 != id2)

  c2 <- map2(cor_vals$id1, cor_vals$id2, function(x, y){
    X_1_2 <- filter(dat_1_2, id.exposure == x & id.outcome == y) %>%
             rename(beta1 = beta.exposure, beta2 = beta.outcome,
                    maf1 = eaf.exposure, maf2 = eaf.outcome,
                    se1 = se.exposure, se2 = se.outcome) %>%
             select(SNP, beta1, beta2, maf1, maf2, se1, se2)
    X_2_1 <- filter(dat_2_1, id.exposure == y & id.outcome == x) %>%
      rename(beta2 = beta.exposure, beta1 = beta.outcome,
             maf2 = eaf.exposure, maf1 = eaf.outcome,
             se2 = se.exposure, se1 = se.outcome) %>%
      select(SNP, beta1, beta2, maf1, maf2, se1, se2) %>%
      filter(!SNP %in% X_1_2$SNP)
    X <- bind_rows(X_1_2, X_2_1)
    r <- with(X, c(cor(beta1, beta2), cor(maf1, maf2), cor(1/se1, 1/se2, method = "spearman")))
    return(r)
  })
  cor_vals$cor <- map(c2, 1) %>% unlist()
  cor_vals$maf_cor <- map(c2, 2) %>% unlist()
  cor_vals$seinv_cor <- map(c2, 3) %>% unlist()
  return(list(mr12 = m_1_2, mr21 = m_2_1, cor = cor_vals))
}


cor_groups <- function(cor, thresh){
  cor.save <- cor
  cor.save$g1 <- 0
  cor.save$g2 <- 0
  grps <- list()
  j <- 1
  while(sum(abs(cor$cor) > thresh) > 0){
    ii <- which(abs(cor$cor) > thresh)[1]
    ids <- c(cor$id1[ii], cor$id2[ii])
    done <- FALSE
    while(!done){
      z <- filter(cor, abs(cor) > thresh & (id1 %in% ids | id2 %in% ids ) )
      if(nrow(z) == 0){
        done <- TRUE
      }else{
        new_ids <- unique(z$id1, z$id2)
        ids <- unique(c(ids, new_ids))
        cor <- filter(cor, !(abs(cor) > thresh & (id1 %in% ids | id2 %in% ids ) ))
      }
    }
    grps[[j]] <- ids
    cor.save$g1[cor.save$id1 %in% ids] <- j
    cor.save$g2[cor.save$id2 %in% ids] <- j
    j <- j + 1
  }
  return(list(g = grps, cor = cor.save))
}
