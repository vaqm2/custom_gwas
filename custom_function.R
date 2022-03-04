#!/usr/bin/env Rscript

require(broom)
require(dplyr)

custom_fit   = function (x_df) {
    fit      = lm(data = x_df, 
                  PHENOTYPE ~ DOSAGE + pc1 + pc2 + pc3 + pc4 + pc5,
                  weights = WEIGHT)
    fit_df   = tidy(fit)
    estimate = fit_df %>% filter(term == "DOSAGE") %>% select(estimate)
    SE       = fit_df %>% filter(term == "DOSAGE") %>% select(std.error)
    p        = fit_df %>% filter(term == "DOSAGE") %>% select(p.value)
    
    return(estimate, SE, p)
}