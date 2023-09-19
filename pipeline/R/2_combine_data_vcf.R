library(VariantAnnotation)
library(gwasvcf)
library(dplyr)
library(rlang)
library(readr)
library(purrr)
library(stringr)
library(ebnm)

args <- commandArgs(trailingOnly=TRUE)
c <- args[1]
gwas_info_file <- args[2]
dir <- args[3]
nmiss_thresh <- args[4]
out <- args[5]
nm_out <- args[6]

info <- read_csv(gwas_info_file)

fulldat <- map(seq(nrow(info)),   function(i){
                        f <- paste0(dir, info$name[i], ".vcf.bgz")
                        n <- info$name[i]
                        v <- query_chrompos_file(paste0(c, ":1-536870911"), f)
                        pos_name <- as_name(paste0(n, ".pos"))
                        beta_name <- as_name(paste0(n, ".beta"))
                        se_name <- as_name(paste0(n, ".se"))
                        p_name <- as_name(paste0(n, ".p"))
                        #wp_name <- as_name(paste0(n, ".wp"))
                        dat <- vcf_to_tibble(v) %>%
                                dplyr::rename(snp=rsid) %>%
                                dplyr::mutate(Z  = ES/SE,
                                              P  = 2*pnorm(-abs(Z))) 
                        ii <- which(!is.na(dat$ES))
                        f <- ebnm(x = dat$ES[ii], s = dat$SE[ii], 
                                  prior_family = "point_normal", output = ebnm::output_all());
                        pi0 <- f$fitted_g$pi[1];
                        mu <- f$fitted_g$mean[2];
                        s2 <- f$fitted_g$sd[2]^2;
                        w <- 1-pi0;
                        a <- 1/s2;
                        #wpost <- ebnm:::wpost_normal(x=dat$ES[ii], s=dat$SE[ii], w, a, mu);
                        #dat$wpost <- NA
                        #dat$wpost[ii] <- wpost
                        dat <- dat %>% 
                                dplyr::select(chr:=seqnames,  snp, REF, ALT, #start, Z, SS)
                                              !!pos_name := start,
                                              !!beta_name := ES,
                                              !!se_name := SE, 
                                              !!p_name := P) 
                                              #!!wp_name := wpost)
                        return(dat)
                 }) %>%
       purrr::reduce(full_join, by = c("chr", "snp", "REF", "ALT"))

dup_snps <- fulldat$snp[duplicated(fulldat$snp)]
if(length(dup_snps) > 0){
    fulldat <- filter(fulldat, !snp %in% dup_snps)
}



# Save table of how traits are missing each SNP for LD clumping
miss <- fulldat %>%
        select(ends_with(".beta")) 
    
miss <- rowSums(is.na(miss))

pmin <- fulldat %>%
        select(ends_with(".p")) 
pmin <- apply(pmin[,-1], 1, function(x){min(x, na.rm=TRUE)})
        
#wpmax <- fulldat %>%
#        select(ends_with(".wp")) 
#wpmax <- apply(wpmax[,-1], 1, function(x){max(x, na.rm=TRUE)})

#df <- data.frame(snp = fulldat$snp, pmin = pmin, wpmax = wpmax, miss = miss)
df <- data.frame(snp = fulldat$snp, pmin = pmin, miss = miss)

ix <- which(miss <= nmiss_thresh)

saveRDS(fulldat[ix,], file=out)

saveRDS(df[ix,], file=nm_out)

