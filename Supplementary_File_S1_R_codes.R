
##########################################################################################
					# R CODES USED IN THE ANALYIS#
##########################################################################################

## Load the required packages for the analysis

library(WGCNA)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(readr)
library(readxl)
library(gtsummary)
library(gt)
library(webshot2)
library(webshot)
library(plm)
library(lcmm)
library(rstatix)
library(lme4)
library(memisc)
library(kableExtra)
library(mclogit)
library(pbkrtest)
library(lavaan)
library(lavaanPlot)
library(janitor)

##########################################################################################
		# TABLE 1
##########################################################################################

### The code below generates Table 1 showing baseline characteristics of the study children

HIV_ncc_growth_Table1_recodID %>% 
  tbl_summary(by = hiv_status,
              type = list(#age_adm ~ 'continuous2',
                          discharge_age ~ 'continuous2',
                          muac_disch ~ 'continuous2',
                          haz_disch ~ 'continuous2',
                          waz_disch ~ 'continuous2',
                          whz_disch ~ 'continuous2'),
              statistic = list(all_continuous() ~ c("{median} ({p25} - {p75})")),
              missing =  "no",
              label = list(sex_adm ~ "Males (%)",
                           discharge_age ~ " Discharge age, months",
                           site ~ "Site",
                           nutritional_category ~ "Nutritional status at discharge",
                           diarrhoea ~ "Diarrhoea",
                           muac_disch ~ "MUAC (cm)",
                           haz_disch ~ "HAZ score",
                           waz_disch ~ "WAZ score",
                           whz_disch ~ "WHZ score",
                           oedema_disch ~ "Oedema - Yes",
                           infec_diag_conrmedmalaria_disch_processed ~  "Malaria Positive (RDT)",
                           infec_diag_measles_disch_processed ~ "Measles",
                           infec_diag_sepsis_disch_processed ~ "Sepsis",
                           resp_diag_pulmonarytb_disch ~ "Respiratory Pulmonary TB",
                           resp_diag_lrti_disch_processed ~ "Lower Respiratory Tract Infection",
                           severe_pneumonia ~ "Pneumonia"),
              missing_text = "Missing",
              value = list(oedema_disch ~ "Yes",
                           sex_adm ~ "Male")) %>%
   
  modify_spanning_header(c("stat_1", "stat_2") ~ "**HIV Status**") %>%
  modify_header(label = "**Participant characteristics at discharge**",
                all_stat_cols() ~ "**{level}**\n(n={n})") %>%
  
  add_overall(last = T, col_label = "**Total**\n(n={N})") %>%
  
  modify_footnote(update = everything() ~ NA) %>%
  gtsummary::as_tibble() %>%
  writexl::write_xlsx(., "../Proteomics/Somalogic-2/Significant_modules/NCC_HIV_substudy/Tables/HIV_ncc_table_v1.xlsx")
  
 
##########################################################################################  
	# Growth from discharge through 6 months of post-discharge follow-up
##########################################################################################


## 1. MUAC

ggboxplot(NCC_HIV_substudy_growth_analysis_lonFormat_box_plot,
          x = "TimePoint", y = "muac", fill = "hiv_status", 
          bxp.errorbar = TRUE, bxp.errorbar.width = 0.2, size = 0.25) +
  geom_hline(yintercept = 11.5, color = "darkgreen", linetype = 2, linewidth = 2.5) +
  geom_smooth(aes(group=hiv_status, color = hiv_status), method = "loess", se = TRUE)+
  labs(x = "Time point", y = "MUAC (cm)") +
  scale_fill_startrek(alpha = 1) +
  scale_fill_startrek(alpha = 1) +
  theme(legend.position = "none",
        text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 360))
        


## 2. HAZ

NCC_HIV_substudy_growth_analysis_lonFormat_box_plot %>%
  dplyr::filter(!is.na(haz)) %>%
  ggplot(aes(TimePoint, haz, group=record_id)) + 
  geom_point(size = 0.5) +
  geom_line(aes(color = hiv_status)) + 
  geom_smooth(aes(group=hiv_status), method = "loess", se = TRUE)+
  geom_hline(yintercept = -3, color = "darkgreen", linetype = 2, linewidth = 2.5) +
  labs(x="Time point", y="HAZ score", color = "HIV status") +
  theme_classic() +
  scale_color_startrek()+
   theme(legend.position = "none",
        text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 360)) +
  facet_wrap(~ hiv_status)



## 3. WAZ


NCC_HIV_substudy_growth_analysis_lonFormat_box_plot %>%
  dplyr::filter(!is.na(waz)) %>%
  ggplot(aes(TimePoint, waz, group=record_id)) + 
  geom_point(size = 0.5) +
  geom_line(aes(color = hiv_status)) + 
  geom_smooth(aes(group=hiv_status), method = "loess", se = TRUE)+
  geom_hline(yintercept = -3, color = "darkgreen", linetype = 2, linewidth = 2.5) +
  labs(x="Time point", y="WAZ score", color = "HIV status") +
  theme_classic() +
  scale_color_startrek()+
   theme(legend.position = "none",
        text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 360)) +
  facet_wrap(~ hiv_status)
  
  
## 4. WHZ

NCC_HIV_substudy_growth_analysis_lonFormat_box_plot %>%
  dplyr::filter(!is.na(whz)) %>%
  ggplot(aes(TimePoint, whz, group=record_id)) + 
  geom_point(size = 0.5) +
  geom_line(aes(color = hiv_status)) + 
  geom_smooth(aes(group=hiv_status), method = "loess", se = TRUE)+
  geom_hline(yintercept = -3, color = "darkgreen", linetype = 2, linewidth = 2.5) +
  labs(x="Time point", y="WHZ score", color = "HIV status") +
  theme_classic() +
  scale_color_startrek()+
   theme(legend.position = "none",
        text = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 360)) +
  facet_wrap(~ hiv_status)
  

##########################################################################################
	# Fixed-effects panel model specification analysis
##########################################################################################


muac_plm_table <- tbl_regression(plm(muac ~ TimePoint + muac_disch + TimePoint*hiv_status + site + sex_adm,  
                                     exponentiate = FALSE,
    data = NCC_HIV_substudy_growth_analysis_lonFormat_plm,
    model = "within"),
    conf.level = 0.95) %>%
  bold_p()
  
  
haz_plm_table <- tbl_regression(plm(haz ~ TimePoint + haz_disch + TimePoint*hiv_status + site + sex_adm, 
                                    exponentiate = FALSE,
    data = NCC_HIV_substudy_growth_analysis_lonFormat_plm,
    model = "within"),
    conf.level = 0.95) %>%
  bold_p()
  

waz_plm_table <- tbl_regression(plm(waz ~ TimePoint + waz_disch + TimePoint*hiv_status + site + sex_adm,  exponentiate = FALSE,
    data = NCC_HIV_substudy_growth_analysis_lonFormat_plm,
    model = "within"),
    conf.level = 0.95) %>%
  bold_p()
  

whz_plm_table <- tbl_regression(plm(whz ~ TimePoint + whz_disch + TimePoint*hiv_status + site + sex_adm,  exponentiate = FALSE,
    data = NCC_HIV_substudy_growth_analysis_lonFormat_plm,
    model = "within"),
    conf.level = 0.95) %>%
  bold_p()
  

library(webshot)
  tbl_merge(
    tbls = list(muac_plm_table, waz_plm_table, 
                whz_plm_table, haz_plm_table),
    tab_spanner = c("MUAC", "WAZ", "WHZ", "HAZ")) %>%
  as_gt() %>%
  gt::tab_options(table.font.names = "Calibri",
                  table.font.size = 12,
                  table.background.color = "lightgrey")
                  

##########################################################################################
	# Inverse Probability weighting  (IPW) to balance baseline characteristics
##########################################################################################

# Using age, sex, enrolment site, nutritional strata and clinical syndromes (malaria, diarrhoea and pneumonia) generate the weights to be used for downstream analysis.

# This chunk is to assist in generating the weights. Use the final output file for subsequent analysis. The generated weight added to the Proteins_discharge_modules dataset.


NCC_HIV_substudy_growth_analysis_weights <- NCC_HIV_substudy_growth_analysis %>% 

  dplyr::mutate(hiv_status_recored = dplyr::recode(hiv_status, "negative" = 0, "positive" = 1)) 
  
NCC_HIV_substudy_growth_analysis_weights$sex_adm <- factor(NCC_HIV_substudy_growth_analysis_weights$sex_adm, 
                                                           levels = c("Female", "Male"))

NCC_HIV_substudy_growth_analysis_weights$infec_diag_conrmedmalaria_disch_processed <- factor(NCC_HIV_substudy_growth_analysis_weights$infec_diag_conrmedmalaria_disch_processed, levels = c("No", "Yes"))


NCC_HIV_substudy_growth_analysis_weights$diarrhoea <- factor(NCC_HIV_substudy_growth_analysis_weights$diarrhoea, 
                                                             levels = c("No", "Yes")) 

NCC_HIV_substudy_growth_analysis_weights$nutritional_category <- factor(NCC_HIV_substudy_growth_analysis_weights$nutritional_category, levels = c("NW", "MW", "SW"))


NCC_HIV_substudy_growth_analysis_weights <- NCC_HIV_substudy_growth_analysis_weights %>% mutate(discharge_age = coalesce(NCC_HIV_substudy_growth_analysis_weights$discharge_age,
                                  NCC_HIV_substudy_growth_analysis_weights$age_adm))


NCC_HIV_substudy_growth_analysis_weights$infec_diag_conrmedmalaria_disch_processed <- NCC_HIV_substudy_growth_analysis_weights$infec_diag_conrmedmalaria_disch_processed %>% replace_na("No")

NCC_HIV_substudy_growth_analysis_weights$diarrhoea <- NCC_HIV_substudy_growth_analysis_weights$diarrhoea %>% replace_na("No")


NCC_HIV_substudy_growth_analysis_weights$severe_pneumonia <- factor(NCC_HIV_substudy_growth_analysis_weights$severe_pneumonia, levels = c("No", "Yes"))


missing_disch_age <- NCC_HIV_substudy_growth_analysis_weights %>% 
  dplyr::filter(is.na(discharge_age)) %>% dplyr::select(record_id) 

missing_disch_age <- missing_disch_age %>% dplyr::left_join(., anthropometry_chain, 
                                                            by = "record_id")


hiv_ncc_logit <- glmer(hiv_status_recored ~ discharge_age + sex_adm  +  
                         nutritional_category +
                         infec_diag_conrmedmalaria_disch_processed +
                         diarrhoea +
                         severe_pneumonia + (1|site),
                       data = NCC_HIV_substudy_growth_analysis_weights, 
                       nAGQ = 100,
                       family = binomial(link = "logit"))

weights_hiv_ncc <- predict(hiv_ncc_logit,type = "response", na.action = na.exclude)

NCC_HIV_substudy_growth_analysis_weights$weights <- predict(hiv_ncc_logit, 
                                                            type = "response", 
                                                            na.action = na.exclude)


NCC_HIV_substudy_growth_analysis_weights <- NCC_HIV_substudy_growth_analysis_weights %>% 
  dplyr::rename(malaria = infec_diag_conrmedmalaria_disch_processed, gastroenteritis = infec_diag_gastroenteritis_disch_processed )


NCC_HIV_substudy_growth_analysis_weights <- NCC_HIV_substudy_growth_analysis_weights %>% # generate inverse probability weights
  dplyr::mutate(como_ipw = (hiv_status_recored/weights) + ((1 - hiv_status_recored) / (1 - weights))) 



NCC_HIV_substudy_growth_analysis_ipw <- NCC_HIV_substudy_growth_analysis_weights %>%  # generate inverse probability weights 
  dplyr::mutate(como_ipw_B = if_else(hiv_status_recored == 1, (1/weights), (1/1-weights)))




NCC_HIV_substudy_growth_analysis_ipw <- NCC_HIV_substudy_growth_analysis_ipw %>% 
  dplyr::select(!starts_with(c("heent_", "skin_rash_", "rash_site_", "cns_", 
                               "infec_diag_", "breathing_", "resp_", "gen_diag",
                               "abdom_", "oxysat_", "relac", "bfeed", "againstadvice", "ricket")))
                               
                               

##########################################################################################
# Weighted Correlation Network analysis Using WGCNA
##########################################################################################

powers <- c(c(1:10), seq(from = 10, to = 50, by = 2 ))


## Generate soft thresholding powers 

SoftThreshold <- pickSoftThreshold(protein_data_human_D0_wide_log_auto_scale,
                  powerVector = powers,
                  networkType = "signed",
                  corFnc = "bicor",
                  corOptions = list(maxPOutliers=0.1),
                  verbose = 5)
                  

## Plotting to identify optimal thresholding power                  

SoftThreshold_data <- SoftThreshold$fitIndices

rsqurd <- ggplot(SoftThreshold_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point(size = 0.35) +
  geom_text(nudge_y = 0.05, size = 3.5) +
  geom_hline(yintercept = 0.85, color = "red")+
  #geom_vline(xintercept = 9 , color = "orange")+
  labs(x = "Soft Threshold (power)", 
       y = "Scale free topology model fit, signed R^2",
       title = "Soft thresholding power criterion") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        text = element_text(size = 12, face = "bold"))


mean_connectivity <- ggplot(SoftThreshold_data, 
                            aes(Power, mean.k., label = Power)) +
  geom_point(size = 0.35) +
  geom_text(nudge_y = 0.05, size = 3.5) +
  #geom_vline(xintercept = 9 , color = "orange")+
  labs(x = "Soft Threshold (power)", y = "Mean Connectivity") + 
  theme_classic() +
  theme(text = element_text(size = 12, face = "bold"))

plt_combined_D0 <- ggarrange(rsqurd, 
                          mean_connectivity,
                          nrow = 2)

plt_combined_D0



soft_power = 9


## Using Automatic single-step network construction generate protein network 

bwnt_D0_bicor <- blockwiseModules(protein_data_human_D0_wide_log_auto_scale,
                 maxBlockSize = 20000,
                 TOMType = "signed",
                 corType = "bicor",
                 networkType = "signed",
                 maxPOutliers = 0.1,
                 minModuleSize = 10,
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 deepSplit = 4,
                 numericLabels = FALSE,
                 saveTOMs = TRUE,
                 saveTOMFileBase = "protein_discharge_blckwise_bicor",
                 randomSeed = 12345,
                 verbose = 1)

##########################################################################################
	# Association between HIV status and Protein Modules
##########################################################################################


### The code below iterates through all the 39 proteins modules to estimate the associate between HIV status and the modules

Adjustedmod_summaries <- list()  

Adjustedmodel_tidy<- list() 

AdjustedModel_PValues <- list()

for (col in colnames(module_cols_cov_adj_data)) {
  var = colnames(Proteins_discharge_modules[,col])
  modules = Proteins_discharge_modules[,col]
 
  Adjustedmod_summaries[[col]] <- lmerTest::lmer(modules[[col]] ~ hiv_status_recored  + (1|site),
                                                 weights = como_ipw_B,
                                                 data = Proteins_discharge_modules)
  
}

coef_estimates <- list()



for (col in colnames(module_cols_cov_adj_data)) {
  
  fixed_effects <- fixef(Adjustedmod_summaries[[col]]) # Extract the fixed effects coefficients

  SE <- sqrt(diag(vcov(Adjustedmod_summaries[[col]]))) # Compute the standard errors

  coef_df <- data.frame(    # Combine the coefficient and standard error data into a data frame
    Protein = names(fixed_effects),
    Estimate = fixed_effects,
    SE = SE
  )

  coef_summary <- broom.mixed::tidy(Adjustedmod_summaries[[col]])     # Use broom.mixed to extract p-values
  p_values <- coef_summary$P.value


  coef_df$p_value <- p_values      # Add p-values to the coefficient data frame

 
  coef_estimates[[col]] <- coef_df     # Append the coefficient data frame to the list
}

adj_coefs_df1 <- do.call(rbind, coef_estimates)   # Combine the coefficient data frames into a single data frame



exclude_values <- c("discharge_age", "sexMale", "(Intercept)", "sex_admMale", 
                    "nutritional_categoryMW", "nutritional_categorySW")  # remove unwanted covariates and their associated coefficients
adj_data <- subset(adj_coefs_df1, !Protein %in% exclude_values)


adj_data$Modules <- rownames(adj_data)

adj_data <- adj_data %>%
  separate_wider_delim(Modules, delim = ".", names = c("Modules", "Status")) %>%
  dplyr::select(-c(Protein, Status))
  
  
  
adj_filtered <- adj_data %>%
  group_by(Modules) %>%
  mutate(
    Lower = Estimate - qnorm(0.975) * SE,
    Upper = Estimate + qnorm(0.975) * SE
  ) %>%
  ungroup()  # Calculate the upper and lower error bars
  
#### Adding a significance column

adj_filtered$significance <- ""
adj_filtered[which(adj_filtered$Lower<0 | adj_filtered$Lower>0 & adj_filtered$Upper<0 |adj_filtered$Upper>0),"significance"] <- "Yes"

 

adj_filtered[which(adj_filtered$Lower<0 & adj_filtered$Upper>0),
"significance"] <- "No" 


significantly_associated_hiv_with_baseline_weights_only <- adj_filtered %>% dplyr::filter(significance == "Yes")


### Wrangle the generated file

adj_filtered_arranged <- adj_filtered 

adj_filtered_arranged <- adj_filtered_arranged %>%  
  dplyr::mutate(Modules = substring(Modules, 3))

adj_filtered_arranged <- adj_filtered_arranged %>% 
  dplyr::inner_join(., proteinModuleSizes_recoded, by = c("Modules" = "module"))

  
adj_filtered_arranged_ordered <- adj_filtered_arranged %>% 
  arrange(Estimate)

level_order <- c(paste(adj_filtered_arranged_ordered$module_recoded_PM, sep = ","))


adj_filtered_arranged_ordered$module_recoded_PM <- factor(adj_filtered_arranged_ordered$module_recoded_PM,
                                                 levels = level_order)
## Generate Forest Plot


adj_filtered_arranged_plt <- 
  adj_filtered_arranged_ordered  %>% 
  ggplot(aes(y = module_recoded_PM, col = significance)) + 
  geom_point(aes(x=Estimate), size= 0.35) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.35, linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype="dashed",  linewidth = 0.35 ) +

  theme_light()+
  labs(x="Estimates (95% CI)", y="Modules", title = "HIV (-)   HIV (+)", color = "BP Signf status") +
  theme(plot.title = element_text(hjust = 0.53),
        legend.position = c(0.95, 0.905),
        text = element_text(size = 9.5),
        axis.text = element_text(size = 9.5),
        legend.title  =  element_text(size=5),
        axis.ticks.length=unit(0.05,"cm"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                colour = "white"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                colour = "grey95")) +
  ggsci::scale_color_lancet() +
  theme(legend.position = "none") +
  guides(y.sec=ggh4x::guide_axis_manual(title = element_blank(), 
                                        breaks = adj_filtered_arranged_ordered$module_recoded_PM, 
                                        labels = paste0("n=", adj_filtered_arranged_ordered$module_size))) 
  
adj_filtered_arranged_plt


## Eigenprotein adjacency heatmap

sizeGrWindow(12.5,8.5);
par(cex = 0.85)
plotEigengeneNetworks(Eigenprotein_adjacency, "Eigenprotein adjacency heatmap", 
                      marDendro = c(0,6,2,2), marHeatmap = c(3,3,3,3), 
                      cex.lab = 0.7, xLabelsAngle = 90,
                      cex.lab.x = 0.7,
                      signed = TRUE,
                      plotDendrograms = FALSE, excludeGrey = FALSE,
                      yLabels = names(MET_recoded_transposed_df),
                      xLabels = names(MET_recoded_transposed_df))
                      

