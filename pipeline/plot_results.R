library(ggplot2)
library(dplyr)
library(purrr)
library(gridExtra)


bee_types <- c("1e-5_0.05", "1e-6_0.05", "5e-8_0.05",
               "1e-5_0", "1e-6_0", "5e-8_0")
esmr_types <- c("NULL", "pval-0.00001", "zl-2", "zl-3", "zl-4", "pval-5e-8")
grapple_types <- c("1e-3", "1e-5", "1e-6")


esmr_res <- map_dfr(esmr_types, function(t){
                        r <- readRDS(paste0("results/cad_v3_esmr_", t, ".ldpruned_r20.01_kb10000_seed0.R_ldsc_full.RDS"))
                        x <- data.frame(exposure = r$names[-1], b = r$beta$beta_m, s = r$beta$beta_s, type = paste0("esmrN_", t))
                        x$p <- with(x, 2*pnorm(-abs(b/s)))
                        return(x)
                  })

bee_res <- map_dfr(bee_types, function(t){
                        r <- readRDS(paste0("results/cad_v3_bee_", t, ".ldpruned_r20.01_kb10000_seed0.R_ldsc_full.RDS"))
                        x <- data.frame(exposure = r$names[-1], b = r$CausalEstimates[-1], s = sqrt(diag(r$VCovCausalEstimates))[-1], type = paste0("bee_", t))
                        x$p <- with(x, 2*pnorm(-abs(b/s)))
                        return(x)
                  })

rivw <- readRDS("results/cad_v3_ivw_5e-8.ldpruned_r20.01_kb10000_seed0.RDS")
ivw_res <- rivw[[1]]$result %>% select(id.exposure, b, se) %>% rename(exposure = id.exposure, b = b, s = se) %>% mutate( type = "MV-IVW")
ivw_res2 <- rivw[[2]]$result %>% select(id.exposure, b, se) %>% rename(exposure = id.exposure, b = b, s = se) %>% mutate( type = "MV-IVW Instrument Specific")
ivw_res <- bind_rows(ivw_res, ivw_res2)
ivw_res$p <- with(ivw_res, 2*pnorm(-abs(b/s)))


restot <- bind_rows(esmr_res, bee_res, ivw_res)


namedir <- data.frame(exposure = esmr_res$exposure[1:9],
                      Exposure = c("BMI", "HDL-C", "LDL-C", "Trig", "Fast. Glu.", "DBP", "SBP", "W2HR", "CRP"))
restot <- full_join(restot, namedir)

restot$b[restot$Exposure %in% c("SBP", "DBP")] <- 10*restot$b[restot$Exposure %in% c("SBP", "DBP")]
restot$s[restot$Exposure %in% c("SBP", "DBP")] <- sqrt(10)*restot$s[restot$Exposure %in% c("SBP", "DBP")]

restot$Exposure = factor(restot$Exposure, levels = c("BMI", "W2HR","CRP", "SBP", "DBP", "Fast. Glu.", "Trig", "HDL-C", "LDL-C"))


plt <- restot %>%
       filter(! type %in% c("bee_1e-5_0.05", "bee_1e-6_0.05", "bee_5e-8_0.05")) %>%
       ggplot() +
       geom_vline(xintercept = 0) +
       geom_point(aes(y = Exposure, x = b, color = type,  group = type), position=position_dodge(width = 0.5), size = 3) +
       geom_errorbar(aes(y = Exposure, xmin =b-qnorm(0.975)*s, xmax = b + qnorm(0.975)*s, color = type), position=position_dodge(width = 0.5)) +
       xlab("Estimate") + coord_flip() +
       theme_bw() +
       theme(axis.text.y = element_text(size = 16),
             axis.text.x = element_text(size = 16, angle = 0),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        legend.position = "bottom")



