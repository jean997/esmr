library(dplyr)
library(ieugwasr)

args <- commandArgs(trailingOnly=TRUE)
info <- readRDS(args[1])
r2_thresh <- as.numeric(args[2])
clump_kb <- args[3]
ref_path  <- args[4]
nm_thresh <- as.numeric(args[5])
thresh <- as.numeric(args[6])
thresh_type <- args[7] # can be pvalue or none
seed <- args[8]
out <- args[9]

info <- info %>%
         rename(rsid = snp) 
     
#if(thresh_type == "wpost" & thresh >=0){
#    info <- info %>% 
#         filter(miss <= nm_thresh & wpmax >= thresh)
#    info$pval <- 1-info$wpmax
#}else 
if(thresh_type == "pvalue"){
    info <- info %>% 
         filter(miss <= nm_thresh & pmin <= thresh)
    info$pval <- info$pmin
}else if(thresh_type == "none"){
    set.seed(seed)
    info$pval <- runif(n = nrow(info))
}else{
    stop("Invalid thresh_type")
}

my_clump <- ld_clump(dat = info,
                     clump_r2 = r2_thresh,
                     clump_p = 1,
                     clump_kb = clump_kb,
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = ref_path)

my_clump <- rename(my_clump, snp=rsid)

saveRDS(my_clump, file=out)

