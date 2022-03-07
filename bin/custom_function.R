#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    require(broom, quietly = TRUE)
    require(dplyr, quietly = TRUE)
})


custom_fit    = function (x_df) {
    fit       = lm(data = x_df, 
                   PHENOTYPE ~ DOSAGE + pc1 + pc2 + pc3 + pc4 + pc5,
                   weights = WEIGHT)
    fit_df    = tidy(fit)
    estimate  = fit_df %>% filter(term == "DOSAGE") %>% select(estimate)
    se        = fit_df %>% filter(term == "DOSAGE") %>% select(std.error)
    p         = fit_df %>% filter(term == "DOSAGE") %>% select(p.value)
    result_df = data.frame(estimate, se, p)
    colnames(result_df) = c("estimate", "se", "p")
    return(result_df)
}