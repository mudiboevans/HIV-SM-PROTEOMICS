
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
library(VGAM)
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
                      
                      
##########################################################################################
	# Correlation between absolute effects sizes (Protein-HIV Significance) and intra-module network connectivity
##########################################################################################

### PM37

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(230)

brown_ki_cmbd_protein <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM37")

hub_prot_37 <- brown_ki_cmbd_protein %>% dplyr::filter(Target == "DNS2B")

ggpubr::ggscatter(brown_ki_cmbd_protein, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.45, label.y = 0.22, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM37", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(brown_ki_cmbd_protein, PS.hiv_status_abs > 0.13 & kWithin > 0.95),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
    geom_point(data = hub_prot_37,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.8, fig.height=3.6, warning=FALSE}


brown <- protein_annot_file %>% 
  dplyr::filter(Modules == "brown") 

brown_boxplot <- brown %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


brown_boxplot_SPEF1 <- brown_boxplot %>% dplyr::filter(Target == "DNS2B")

brown_boxplot_SPEF1$hiv_status_recoded <- factor(brown_boxplot_SPEF1$hiv_status_recoded, 
                                                 levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

brown_boxplot_SPEF1$hiv_status_recoded <- as.factor(brown_boxplot_SPEF1$hiv_status_recoded)

t_test_p37 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = brown_boxplot_SPEF1)

t_test_p37 <- t_test_p37 %>% mutate(y.position = 10.1)


ggboxplot(brown_boxplot_SPEF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "DNS2B levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p37, label = "p.adj") +
  scale_color_startrek()+
  scale_fill_startrek()
  
```

### PM31

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(231)

p31 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM31")

hub_prot_31 <- p31 %>% dplyr::filter(Target == "SEM4F")


ggpubr::ggscatter(p31, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.45, label.y = 0.18, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM31", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p31, PS.hiv_status_abs < 0.05 & kWithin > 0.78),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  #geom_point(color = ifelse(p31$Target == "SEM4F", "red", "black"), size = 1.5) +
  geom_point(data = hub_prot_31,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


magenta <- protein_annot_file %>% 
  dplyr::filter(Modules == "magenta") 

magenta_boxplot <- magenta %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


magenta_boxplot_LRRT4 <- magenta_boxplot %>% dplyr::filter(Target == "SEM4F")


magenta_boxplot_LRRT4$hiv_status_recoded <- factor(magenta_boxplot_LRRT4$hiv_status_recoded, 
                                                 levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

magenta_boxplot_LRRT4$hiv_status_recoded <- as.factor(magenta_boxplot_LRRT4$hiv_status_recoded)

t_test_p31 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = magenta_boxplot_LRRT4)

t_test_p31 <- t_test_p31 %>% mutate(y.position = 9.3)


ggboxplot(magenta_boxplot_LRRT4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "SEM4F levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p31, label = "p.adj") +
  scale_fill_startrek()+
  scale_color_startrek()

```

### PM32

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(231)

p32 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM32")


hub_prot_32 <- p32 %>% dplyr::filter(Target == "FGF-16")

ggpubr::ggscatter(p32, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.15, label.y = 0.165, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM32", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p32, PS.hiv_status_abs > 0.11 & kWithin > 0.9),
                           aes(label = Target),
                           size = 3,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.125, nudge_y = -0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 # geom_point(color = ifelse(p32$Target == "FGF-16", "red", "black"), size = 1.5) +
  
   geom_point(data = hub_prot_32,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)
  
```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


pink <- protein_annot_file %>% 
  dplyr::filter(Modules == "pink") 

pink_boxplot <- pink %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


pink_boxplot_FGF16 <- pink_boxplot %>% dplyr::filter(Target == "FGF-16")


pink_boxplot_FGF16$hiv_status_recoded <- factor(pink_boxplot_FGF16$hiv_status_recoded, 
                                                 levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

pink_boxplot_FGF16$hiv_status_recoded <- as.factor(pink_boxplot_FGF16$hiv_status_recoded)

t_test_p32 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = pink_boxplot_FGF16)

t_test_p32 <- t_test_p32 %>% mutate(y.position = 5.15)


ggboxplot(pink_boxplot_FGF16,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "FGF-16 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p32, label = "p.adj") +
  scale_fill_startrek()+
  scale_color_startrek()

```

### PM5

```{r message=FALSE, fig.width=2.55, fig.height=3.4}

set.seed(232)

p5 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM5")

hub_prot_5 <- p5 %>% dplyr::filter(Target == "TPPC3")

ggpubr::ggscatter(p5, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.25, label.y = 0.14, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM5", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p5, PS.hiv_status_abs > 0.1 & kWithin > 0.75),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.05, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  #geom_point(color = ifelse(p5$Target == "TPPC3", "red", "black"), size = 1.5) +
     geom_point(data = hub_prot_5,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.55, fig.height=3.6, warning=FALSE}


sienna3 <- protein_annot_file %>% 
  dplyr::filter(Modules == "sienna3") 

sienna3_boxplot <- sienna3 %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


sienna3_boxplot_HNF4A <- sienna3_boxplot %>% dplyr::filter(Target == "TPPC3")


sienna3_boxplot_HNF4A$hiv_status_recoded <- factor(sienna3_boxplot_HNF4A$hiv_status_recoded, 
                                                  levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

sienna3_boxplot_HNF4A$hiv_status_recoded <- as.factor(sienna3_boxplot_HNF4A$hiv_status_recoded)

t_test_p5 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = sienna3_boxplot_HNF4A)

t_test_p5 <- t_test_p5 %>% mutate(y.position = 5.5)


ggboxplot(sienna3_boxplot_HNF4A,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "TPPC3 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p5, label = "p.adj") +
  scale_fill_startrek() 

```

### PM16

```{r message=FALSE, fig.width=2.65, fig.height=3.4}

set.seed(233)

p16 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM16")

p16 <- p16 %>% dplyr::mutate(Target = dplyr::recode(Target, "FGF-8B" = "FGF8B"))

hub_prot_16 <- p16 %>% dplyr::filter(Target == "FGF8B")

ggpubr::ggscatter(p16, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  
  stat_cor(method  = "pearson", label.x = 0.4, label.y = 0.175, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM16", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
  ggrepel::geom_text_repel(data=subset(p16, PS.hiv_status_abs > 0.1 & kWithin == 1),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.13, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
  #geom_point(color = ifelse(p16$Target == "FGF8B", "red", "black"), size = 1.5) +
  geom_point(data = hub_prot_16,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


darkgrey <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkgrey") 

darkgrey_boxplot <- darkgrey %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkgrey_boxplot_FGF8B <- darkgrey_boxplot %>% dplyr::filter(Target == "FGF-8B")


darkgrey_boxplot_FGF8B$hiv_status_recoded <- factor(darkgrey_boxplot_FGF8B$hiv_status_recoded, 
                                                  levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkgrey_boxplot_FGF8B$hiv_status_recoded <- as.factor(darkgrey_boxplot_FGF8B$hiv_status_recoded)

t_test_p16 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkgrey_boxplot_FGF8B)

t_test_p16 <- t_test_p16 %>% mutate(y.position = 7)


ggboxplot(darkgrey_boxplot_FGF8B,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "FGF-8B levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p16, label = "p.adj") +
  scale_fill_startrek() 

```

### PM6

```{r message=FALSE, fig.width=2.65, fig.height=3.4}

set.seed(234)

p6 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM6")

p6 <- p6 %>% dplyr::mutate(Target = dplyr::recode(Target, "OA synthetase" = "OAS"))

hub_prot_6 <- p6 %>% dplyr::filter(Target == "SC5A8")

ggpubr::ggscatter(p6, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.5, label.y = 0.065, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM6", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
  ggrepel::geom_text_repel(data=subset(p6, PS.hiv_status_abs < 0.025 & kWithin > 0.75),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.000145, nudge_y = 0.0055,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  #geom_point(color = ifelse(p6$Target == "SC5A8", "red", "black"), size = 1.5) +
      geom_point(data = hub_prot_6,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.61, fig.height=3.6, warning=FALSE}


darkmagenta <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkmagenta") 

darkmagenta_boxplot <- darkmagenta %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkmagenta_boxplot_SC5A8 <- darkmagenta_boxplot %>% dplyr::filter(Target == "SC5A8")


darkmagenta_boxplot_SC5A8$hiv_status_recoded <- factor(darkmagenta_boxplot_SC5A8$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkmagenta_boxplot_SC5A8$hiv_status_recoded <- as.factor(darkmagenta_boxplot_SC5A8$hiv_status_recoded)

t_test_p6 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkmagenta_boxplot_SC5A8)

t_test_p6 <- t_test_p6 %>% mutate(y.position = 4.2)


ggboxplot(darkmagenta_boxplot_SC5A8,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "SC5A8 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p6, label = "p.adj") +
  scale_fill_startrek() 

```

### PM39

```{r message=FALSE, fig.width=2.70, fig.height=3.4}

set.seed(234)

p39 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM39")

hub_prot_39 <- p39 %>% dplyr::filter(AptName == "seq.20139.57")


ggpubr::ggscatter(p39, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95, 
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.3, label.y = 0.18, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM39", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 11, face = "bold"),
        axis.text = element_text(face = "plain")) +
  ggrepel::geom_text_repel(data=subset(p39, AptName == "seq.20139.57"),
                           aes(label = Target),
                           size = 3,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = -0.08, nudge_y = - 0.0085,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  geom_point(data = hub_prot_39,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.61, fig.height=3.6, warning=FALSE}


turquoise <- protein_annot_file %>% 
  dplyr::filter(Modules == "turquoise") 

turquoise_boxplot <- turquoise %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


turquoise_boxplot_SH3G2 <- turquoise_boxplot %>% dplyr::filter(Target == "SH3G2")


turquoise_boxplot_SH3G2$hiv_status_recoded <- factor(turquoise_boxplot_SH3G2$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

turquoise_boxplot_SH3G2$hiv_status_recoded <- as.factor(turquoise_boxplot_SH3G2$hiv_status_recoded)

t_test_p39 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = turquoise_boxplot_SH3G2)

t_test_p39 <- t_test_p39 %>% mutate(y.position = 3.5)


ggboxplot(turquoise_boxplot_SH3G2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "SH3G2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p39, label = "p.adj") +
  scale_fill_startrek() 

```

### PM8

```{r message=FALSE, fig.width=2.5, fig.height=3.4}

set.seed(235)

p8 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM8")


p8 <- p8 %>% dplyr::mutate(Target = dplyr::recode(Target, "complement factor H-related 5" = "CFHR5"))

hub_prot_8 <- p8 %>% dplyr::filter(AptName == "seq.16055.3")

ggpubr::ggscatter(p8, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.5, label.y = 0.129, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM8", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p8, PS.hiv_status_abs > 0.1 & kWithin > 0.9),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.145, nudge_y = -0.00195,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p8$AptName == "seq.16055.3", "red", "black"), size = 1.8) +
  
      geom_point(data = hub_prot_8,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


paleturquoise <- protein_annot_file %>% 
  dplyr::filter(Modules == "paleturquoise") 

paleturquoise_boxplot <- paleturquoise %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


paleturquoise_boxplot_CFHR5 <- paleturquoise_boxplot %>% dplyr::filter(AptName == "seq.16055.3")

paleturquoise_boxplot_CFHR5$hiv_status_recoded <- factor(paleturquoise_boxplot_CFHR5$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

paleturquoise_boxplot_CFHR5$hiv_status_recoded <- as.factor(paleturquoise_boxplot_CFHR5$hiv_status_recoded)

t_test_p8 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = paleturquoise_boxplot_CFHR5)

t_test_p8 <- t_test_p8 %>% mutate(y.position = 2.2)


ggboxplot(paleturquoise_boxplot_CFHR5,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "CFHR5 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p8, label = "p.adj") +
  scale_fill_startrek() 


```

### PM2

```{r message=FALSE, fig.width=2.6, fig.height=3.0}

set.seed(236)

p2 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM2")

hub_prot_2 <- p2 %>% dplyr::filter(Target == "HCE000104")

ggpubr::ggscatter(p2, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.25, label.y = 0.1, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM2", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p2, PS.hiv_status_abs > 0.11 & kWithin > 0.95),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = -0.185, nudge_y = -0.00195,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p2$Target == "HCE000104", "red", "black"), size = 1.8) +
     geom_point(data = hub_prot_2,
             aes(x = kWithin, y = PS.hiv_status_abs),
             color = "red", size = 1.25)
             
```

#### Box plot

```{r fig.width=1.58, fig.height=3.6, warning=FALSE}


plum1 <- protein_annot_file %>% 
  dplyr::filter(Modules == "plum1") 

plum1_boxplot <- plum1 %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


plum1_boxplot_HCE000104 <- plum1_boxplot %>% dplyr::filter(Target == "HCE000104")


plum1_boxplot_HCE000104$hiv_status_recoded <- factor(plum1_boxplot_HCE000104$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

plum1_boxplot_HCE000104$hiv_status_recoded <- as.factor(plum1_boxplot_HCE000104$hiv_status_recoded)

t_test_p2 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = plum1_boxplot_HCE000104)

t_test_p2 <- t_test_p2 %>% mutate(y.position = 2.5)


ggboxplot(plum1_boxplot_HCE000104,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "HCE000104 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p2, label = "p.adj") +
  scale_fill_startrek() 

```

### PM30

```{r message=FALSE, fig.width=2.65, fig.height=3.1}

set.seed(236)

p30 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM30")

#p30 <- p30 %>% dplyr::mutate(Target = dplyr::recode(Target, "Lectin, mannose-binding 2" = "MBL2"))

hub_prot_30 <- p30 %>% dplyr::filter(Target == "OPT")

ggpubr::ggscatter(p30, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0., label.y = 0.155, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM30", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) + 
    ggrepel::geom_text_repel(data=subset(p30, PS.hiv_status_abs > 0.05 & kWithin > 0.9),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.0185, nudge_y = -0.0095,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p30$Target == "OPT", "red", "black"), size = 1.5) +
  geom_point(data = hub_prot_30,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


purple <- protein_annot_file %>% 
  dplyr::filter(Modules == "purple") 

purple_boxplot <- purple %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


purple_boxplot_OPT <- purple_boxplot %>% dplyr::filter(Target == "OPT")


purple_boxplot_OPT$hiv_status_recoded <- factor(purple_boxplot_OPT$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

purple_boxplot_OPT$hiv_status_recoded <- as.factor(purple_boxplot_OPT$hiv_status_recoded)

t_test_p30 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = purple_boxplot_OPT)

t_test_p30 <- t_test_p30 %>% mutate(y.position = 7)


ggboxplot(purple_boxplot_OPT,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "OPT levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p30, label = "p.adj") +
  scale_fill_startrek() 

```

### PM10

```{r message=FALSE, fig.width=2.65, fig.height=3.1}

set.seed(236)

PM10 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM10")

hub_prot_10 <- PM10 %>% dplyr::filter(Target == "CCD43")

ggpubr::ggscatter(PM10, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.3, label.y = 0.1, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM10", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(PM10, AptName == "seq.21595.8"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = -0.085, nudge_y = 0.0055,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(PM10$Target == "CCD43", "red", "black"), size = 1.5) +
  
  geom_point(data = hub_prot_10,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


violet <- protein_annot_file %>% 
  dplyr::filter(Modules == "violet") 

violet_boxplot <- violet %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


violet_boxplot_CCD43 <- violet_boxplot %>% dplyr::filter(Target == "CCD43")


violet_boxplot_CCD43$hiv_status_recoded <- factor(violet_boxplot_CCD43$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

violet_boxplot_CCD43$hiv_status_recoded <- as.factor(violet_boxplot_CCD43$hiv_status_recoded)

t_test_p10 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = violet_boxplot_CCD43)

t_test_p10 <- t_test_p10 %>% mutate(y.position = 8.3)


ggboxplot(violet_boxplot_CCD43,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "CCD43 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p10, label = "p.adj") +
  scale_fill_startrek() 

```

### PM25

```{r message=FALSE, fig.width=2.65, fig.height=3.1}

set.seed(236)

p25 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM25")

hub_prot_25 <- p25 %>% dplyr::filter(Target == "CQ10A")

ggpubr::ggscatter(p25, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.5, label.y = 0.13, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM25", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p25, PS.hiv_status_abs < 0.02 & kWithin > 0.77),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0055,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p25$Target == "CQ10A", "red", "black"), size = 1.8) +
  
    geom_point(data = hub_prot_25,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.65, fig.height=3.6, warning=FALSE}


midnightblue <- protein_annot_file %>% 
  dplyr::filter(Modules == "midnightblue") 

midnightblue_boxplot <- midnightblue %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


midnightblue_boxplot_CQ10A <- midnightblue_boxplot %>% dplyr::filter(Target == "CQ10A")


midnightblue_boxplot_CQ10A$hiv_status_recoded <- factor(midnightblue_boxplot_CQ10A$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

midnightblue_boxplot_CQ10A$hiv_status_recoded <- as.factor(midnightblue_boxplot_CQ10A$hiv_status_recoded)

t_test_p25 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = midnightblue_boxplot_CQ10A)

t_test_p25 <- t_test_p25 %>% mutate(y.position = 7)


ggboxplot(midnightblue_boxplot_CQ10A,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "CQ10A levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p25, label = "p.adj") +
  scale_fill_startrek() 

```

### PM33

```{r message=FALSE, fig.width=2.65, fig.height=3.01}

set.seed(186)

p33 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM33")

hub_prot_33 <- p33 %>% dplyr::filter(Target == "STMN4")
ggpubr::ggscatter(p33, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.85,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.7)) +

  stat_cor(method  = "pearson", label.x = 0.5, label.y = 0.165, size = 4.5,
           label.sep='\n') +
  labs(x = "Connectivity within PM33", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p33, PS.hiv_status_abs > 0.1 & kWithin > 0.9),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.125, nudge_y = 0.0095,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p33$Target == "STMN4", "red", "black"), size = 1.8) +
  
    geom_point(data = hub_prot_33,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


black <- protein_annot_file %>% 
  dplyr::filter(Modules == "black") 

black_boxplot <- black %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


black_boxplot_STMN4 <- black_boxplot %>% dplyr::filter(Target == "STMN4")


black_boxplot_STMN4$hiv_status_recoded <- factor(black_boxplot_STMN4$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

black_boxplot_STMN4$hiv_status_recoded <- as.factor(black_boxplot_STMN4$hiv_status_recoded)

t_test_p33 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = black_boxplot_STMN4)

t_test_p33 <- t_test_p33 %>% mutate(y.position = 4.05)


ggboxplot(black_boxplot_STMN4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "STMN4 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p33, label = "p.adj") +
  scale_fill_startrek() 

```

### PM34

```{r message=FALSE, fig.width=2.65, fig.height=3.01}

set.seed(186)

p34 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM34")

p34 <- p34 %>% dplyr::mutate(Target = dplyr::recode(Target, "nectin-1 gamma:CD" = "NECTIN1"))

hub_prot_34 <- p34 %>% dplyr::filter(Target == "SATB1")

ggpubr::ggscatter(p34, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.85,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.7)) +

  stat_cor(method  = "pearson", label.x = 0.35, label.y = 0.165, size = 4.5,
           label.sep='\n') +
  labs(x = "Connectivity within PM34", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p34, Target == "SATB1"),
                           aes(label = Target),
                           size = 3.5,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.125, nudge_y = -0.0095,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p34$Target == "SATB1", "red", "black"), size = 1.8) +
  
   geom_point(data = hub_prot_34,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


red <- protein_annot_file %>% 
  dplyr::filter(Modules == "red") 

red_boxplot <- red %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


red_boxplot_SATB1 <- red_boxplot %>% dplyr::filter(Target == "SATB1")


red_boxplot_SATB1$hiv_status_recoded <- factor(red_boxplot_SATB1$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

red_boxplot_SATB1$hiv_status_recoded <- as.factor(red_boxplot_SATB1$hiv_status_recoded)

t_test_p34 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = red_boxplot_SATB1)

t_test_p34 <- t_test_p34 %>% mutate(y.position = 9.55)


ggboxplot(red_boxplot_SATB1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "SATB1 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p34, label = "p.adj") +
  scale_fill_startrek() 

```

### PM17

```{r message=FALSE, fig.width=2.6, fig.height=3.0}

set.seed(236)

p17 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM17")

p17 <- p17 %>% dplyr::mutate(Target = dplyr::recode(Target, "Lectin, mannose-binding 2" = "MBL2"))

hub_prot_17 <- p17 %>% dplyr::filter(Target == "MBL2")

ggpubr::ggscatter(p17, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.5, label.y = 0.13, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM17", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p17, PS.hiv_status_abs < 0.11 & kWithin > 0.99),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0055,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p17$Target == "MBL2", "red", "black"), size = 1.8) +
  
   geom_point(data = hub_prot_17,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


darkturquoise <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkturquoise") 

darkturquoise_boxplot <- darkturquoise %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkturquoise_boxplot_Lectin <- darkturquoise_boxplot %>% dplyr::filter(Target == "Lectin, mannose-binding 2")


darkturquoise_boxplot_Lectin$hiv_status_recoded <- factor(darkturquoise_boxplot_Lectin$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkturquoise_boxplot_Lectin$hiv_status_recoded <- as.factor(darkturquoise_boxplot_Lectin$hiv_status_recoded)

t_test_p17 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkturquoise_boxplot_Lectin)

t_test_p17 <- t_test_p17 %>% mutate(y.position = 5)


ggboxplot(darkturquoise_boxplot_Lectin,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "MBL2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p17, label = "p.adj") +
  scale_fill_startrek() 

```

### PM38

```{r message=FALSE, fig.width=2.5, fig.height=3.0}

set.seed(237)

p38 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM38")

hub_prot_38 <- p38 %>% dplyr::filter(AptName == "seq.7190.50")

ggpubr::ggscatter(p38, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.4, label.y = 0.25, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM38", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p38, PS.hiv_status_abs > 0.11 & kWithin > 0.99),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0095,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p38$AptName == "seq.7190.50", "red", "black"), size = 1.8) +
  
   geom_point(data = hub_prot_38,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


blue <- protein_annot_file %>% 
  dplyr::filter(Modules == "blue") 

blue_boxplot <- blue %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


blue_boxplot_EPHA4 <- blue_boxplot %>% dplyr::filter(AptName ==  "seq.7190.50")


blue_boxplot_EPHA4$hiv_status_recoded <- factor(blue_boxplot_EPHA4$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

blue_boxplot_EPHA4$hiv_status_recoded <- as.factor(blue_boxplot_EPHA4$hiv_status_recoded)

t_test_p38 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = blue_boxplot_EPHA4)

t_test_p38 <- t_test_p38 %>% mutate(y.position = 3)


ggboxplot(blue_boxplot_EPHA4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "EPHA4 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p38, label = "p.adj") +
  scale_fill_startrek() 


```

### PM14

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(238)

p14 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM14")

hub_prot_14 <- p14 %>% dplyr::filter(Target == "DPH5")

ggpubr::ggscatter(p14, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.55, label.y = 0.15, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM14", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p14, PS.hiv_status_abs > 0.09 & kWithin > 0.99),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.105, nudge_y = 0.0095,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p14$Target == "DPH5", "red", "black"), size = 1.8) +
    
   geom_point(data = hub_prot_14,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

### Box plot

```{r fig.width=1.65, fig.height=3.6, warning=FALSE}


darkorange <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkorange") 

darkorange_boxplot <- darkorange %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkorange_boxplot_DPH5 <- darkorange_boxplot %>% dplyr::filter(Target ==  "DPH5")


darkorange_boxplot_DPH5$hiv_status_recoded <- factor(darkorange_boxplot_DPH5$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkorange_boxplot_DPH5$hiv_status_recoded <- as.factor(darkorange_boxplot_DPH5$hiv_status_recoded)

t_test_p14 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkorange_boxplot_DPH5)

t_test_p14 <- t_test_p14 %>% mutate(y.position = 6.35)


ggboxplot(darkorange_boxplot_DPH5,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "DPH5 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p14, label = "p.adj") +
  scale_fill_startrek() 

```

### PM1

```{r message=FALSE, fig.width=2.5, fig.height=3.2}

set.seed(239)

p1 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM1")

hub_prot_1 <- p1 %>% dplyr::filter(Target == "Transferrin")

ggpubr::ggscatter(p1, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.25, label.y = 0.12, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM1", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p1, PS.hiv_status_abs > 0.1 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.25, nudge_y = -0.00395,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p1$Target == "Transferrin", "red", "black"), size = 1.8) +
  
     geom_point(data = hub_prot_1,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)
              
```

#### Box plot

```{r fig.width=1.58, fig.height=3.6, warning=FALSE}


orangered4 <- protein_annot_file %>% 
  dplyr::filter(Modules == "orangered4") 

orangered4_boxplot <- orangered4 %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkorange_boxplot_Transferrin <- orangered4_boxplot %>% dplyr::filter(Target ==  "Transferrin")



darkorange_boxplot_Transferrin$hiv_status_recoded <- factor(darkorange_boxplot_Transferrin$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkorange_boxplot_Transferrin$hiv_status_recoded <- as.factor(darkorange_boxplot_Transferrin$hiv_status_recoded)

t_test_p1 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkorange_boxplot_Transferrin)

t_test_p1 <- t_test_p1 %>% mutate(y.position = 4.4)


ggboxplot(darkorange_boxplot_Transferrin,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "Transferrin levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p1, label = "p.adj") +
  scale_fill_startrek() 

```

### PM12

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(240)

p12 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM12")


hub_prot_12 <- p12 %>% dplyr::filter(Target == "MOB3B")

ggpubr::ggscatter(p12, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.25, label.y = 0.19, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM12", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p12, PS.hiv_status_abs > 0.155 & kWithin > 0.9),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = -0.01225,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p12$Target == "MOB3B", "red", "black"), size = 1.8) +
     geom_point(data = hub_prot_12,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


skyblue <- protein_annot_file %>% 
  dplyr::filter(Modules == "skyblue") 

skyblue_boxplot <- skyblue %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


skyblue_boxplot_MOB3B <- skyblue_boxplot %>% dplyr::filter(Target ==  "MOB3B")


skyblue_boxplot_MOB3B$hiv_status_recoded <- factor(skyblue_boxplot_MOB3B$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

skyblue_boxplot_MOB3B$hiv_status_recoded <- as.factor(skyblue_boxplot_MOB3B$hiv_status_recoded)

t_test_p12 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = skyblue_boxplot_MOB3B)

t_test_p12 <- t_test_p12 %>% mutate(y.position = 4.2)


ggboxplot(skyblue_boxplot_MOB3B,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "MOB3B levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p12, label = "p.adj") +
  scale_fill_startrek() 

```

### PM18

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(241)

p18 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM18")

p18 <- p18 %>% dplyr::filter(Target != "No Protein")


p18 <- p18 %>% dplyr::mutate(Target = dplyr::recode(Target, "Carbonic anhydrase VII" = "CA7"))

hub_prot_18 <- p18 %>% dplyr::filter(Target == "CA7")

ggpubr::ggscatter(p18, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.38, label.y = 0.145, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM18", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p18, PS.hiv_status_abs > 0.14 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = -0.010,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p18$Target == "CA7", "red", "black"), size = 1.8) +
     geom_point(data = hub_prot_18,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)
     
```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}


darkgreen <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkgreen") 

darkgreen_boxplot <- darkgreen %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkgreen_boxplot_CA7 <- darkgreen_boxplot %>% dplyr::filter(Target ==  "Carbonic anhydrase VII")

darkgreen_boxplot_CA7$hiv_status_recoded <- factor(darkgreen_boxplot_CA7$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkgreen_boxplot_CA7$hiv_status_recoded <- as.factor(darkgreen_boxplot_CA7$hiv_status_recoded)

t_test_p18 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkgreen_boxplot_CA7)

t_test_p18 <- t_test_p18 %>% mutate(y.position = 8.8)


ggboxplot(darkgreen_boxplot_CA7,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "CA7 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p18, label = "p.adj") +
  scale_fill_startrek() 

```

### PM27

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(242)

p27 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM27")

hub_prot_27 <- p27 %>% dplyr::filter(Target == "FR1OP")

ggpubr::ggscatter(p27, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.099, label.y = 0.165, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM27", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p27, PS.hiv_status_abs > 0.15 & kWithin > 0.75),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = -0.010,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p27$Target == "FR1OP", "red", "black"), size = 1.8) +
  
       geom_point(data = hub_prot_27,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}


salmon <- protein_annot_file %>% 
  dplyr::filter(Modules == "salmon") 

salmon_boxplot <- salmon %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


salmon_boxplot_FR1OP <- salmon_boxplot %>% dplyr::filter(Target ==  "FR1OP")


salmon_boxplot_FR1OP$hiv_status_recoded <- factor(salmon_boxplot_FR1OP$hiv_status_recoded, 
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

salmon_boxplot_FR1OP$hiv_status_recoded <- as.factor(salmon_boxplot_FR1OP$hiv_status_recoded)

t_test_p27 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = salmon_boxplot_FR1OP)

t_test_p27<- t_test_p27 %>% mutate(y.position = 10.4)


ggboxplot(salmon_boxplot_FR1OP,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "FR1OP levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p27, label = "p.adj") +
  scale_fill_startrek() 

```

## PM3

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(243)

p3 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM3")


p3 <- p3 %>% dplyr::mutate(Target = dplyr::recode(Target, "Kininogen, HMW, Two Chain" = "KNG2"))

hub_prot_3 <- p3 %>% dplyr::filter(Target == "KNG2")

ggpubr::ggscatter(p3, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.4, label.y = 0.165, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM3", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p3, Target == "KNG2"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = 0.00495,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p3$Target == "KNG2", "red", "black"), size = 1.8) +
  
  
       geom_point(data = hub_prot_3,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


skyblue3 <- protein_annot_file %>% 
  dplyr::filter(Modules == "skyblue3") 

skyblue3_boxplot <- skyblue3 %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


skyblue3_boxplot_KNG2 <- skyblue3_boxplot %>% dplyr::filter(Target ==  "Kininogen, HMW, Two Chain")


skyblue3_boxplot_KNG2$hiv_status_recoded <- factor(skyblue3_boxplot_KNG2$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))


my_comparisons <- list( c("HIV-", "HIV+"))

skyblue3_boxplot_KNG2$hiv_status_recoded <- as.factor(skyblue3_boxplot_KNG2$hiv_status_recoded)

t_test_p3 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = skyblue3_boxplot_KNG2)

t_test_p3<- t_test_p3 %>% mutate(y.position = 3.4)


ggboxplot(skyblue3_boxplot_KNG2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "KNG2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p3, label = "p.adj") +
  scale_fill_startrek() 

```

### PM7

```{r message=FALSE, fig.width=2.5, fig.height=3.2}

set.seed(244)

p7 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM7")


hub_prot_7 <- p7 %>% dplyr::filter(Target == "IGFBP-3")

ggpubr::ggscatter(p7, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.1, label.y = 0.185, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM7", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p7, Target == "IGFBP-3"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.19, nudge_y = 0.00895,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p7$Target == "IGFBP-3", "red", "black"), size = 1.8) +
    
       geom_point(data = hub_prot_7,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}

darkolivegreen <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkolivegreen") 

darkolivegreen_boxplot <- darkolivegreen %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkolivegreen_boxplot_IGF1 <- darkolivegreen_boxplot %>% dplyr::filter(Target ==  "IGFBP-3")


darkolivegreen_boxplot_IGF1$hiv_status_recoded <- factor(darkolivegreen_boxplot_IGF1$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))


my_comparisons <- list( c("HIV-", "HIV+"))

darkolivegreen_boxplot_IGF1$hiv_status_recoded <- as.factor(darkolivegreen_boxplot_IGF1$hiv_status_recoded)

t_test_p7 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkolivegreen_boxplot_IGF1)

t_test_p7 <- t_test_p7 %>% mutate(y.position = 2.35)


ggboxplot(darkolivegreen_boxplot_IGF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "IGFBP-3 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p7, label = "p.adj") +
  scale_fill_startrek() 

```


### PM19

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(245)

p19 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM19")

hub_prot_19 <- p19 %>% dplyr::filter(Target == "SOD")

ggpubr::ggscatter(p19, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.1, label.y = 0.185, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM19", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p19, PS.hiv_status_abs > 0.1 & kWithin > 0.95),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.19, nudge_y = 0.00895,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p19$Target == "SOD", "red", "black"), size = 1.8) +
  
       geom_point(data = hub_prot_19,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)
              
```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}

darkred <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkred") 

darkred_boxplot <- darkred %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkred_boxplot_SOD <- darkred_boxplot %>% dplyr::filter(Target ==  "SOD")


darkred_boxplot_SOD$hiv_status_recoded <- factor(darkred_boxplot_SOD$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkred_boxplot_SOD$hiv_status_recoded <- as.factor(darkred_boxplot_SOD$hiv_status_recoded)

t_test_p19 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkred_boxplot_SOD)

t_test_p19 <- t_test_p19 %>% mutate(y.position = 6.0)


ggboxplot(darkred_boxplot_SOD,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "SOD levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p19, label = "p.adj") +
  scale_fill_startrek() 

```

### PM26

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(246)

p26 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM26")


hub_prot_26 <- p26 %>% dplyr::filter(Target == "COL11A2")

ggpubr::ggscatter(p26, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.6, label.y = 0.22, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM26", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p26, PS.hiv_status_abs > 0.12 & kWithin > 0.95),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = 0.01495,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p26$Target == "COL11A2", "red", "black"), size = 1.8) +
  
       geom_point(data = hub_prot_26,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}

cyan <- protein_annot_file %>% 
  dplyr::filter(Modules == "cyan") 

cyan_boxplot <- cyan %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


cyan_boxplot_COL11A2 <- cyan_boxplot %>% dplyr::filter(Target ==  "COL11A2")

cyan_boxplot_COL11A2$hiv_status_recoded <- factor(cyan_boxplot_COL11A2$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

cyan_boxplot_COL11A2$hiv_status_recoded <- as.factor(cyan_boxplot_COL11A2$hiv_status_recoded)

t_test_p26 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = cyan_boxplot_COL11A2)

t_test_p26 <- t_test_p26 %>% mutate(y.position = 2.95)


ggboxplot(cyan_boxplot_COL11A2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "COL11A2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p26, label = "p.adj") +
  scale_fill_startrek() 

```

### PM23

```{r message=FALSE, fig.width=2.8, fig.height=3.4}

set.seed(247)

p23 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM23")


p23 <- p23 %>% dplyr::mutate(Target = dplyr::recode(Target, "IFN-lambda 2" = "IFNL2"))

hub_prot_23 <- p23 %>% dplyr::filter(Target == "Enterokinase")

ggpubr::ggscatter(p23, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.15, label.y = 0.185, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM23", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p23, Target == "Enterokinase"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.19, nudge_y = 0.01295,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p23$Target == "Enterokinase", "red", "black"), size = 1.8) +
  
       geom_point(data = hub_prot_23,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}

grey60 <- protein_annot_file %>% 
  dplyr::filter(Modules == "grey60") 

grey60_boxplot <- grey60 %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


grey60_boxplot_Enterokinase <- grey60_boxplot %>% dplyr::filter(Target ==  "Enterokinase")



grey60_boxplot_Enterokinase$hiv_status_recoded <- factor(grey60_boxplot_Enterokinase$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

grey60_boxplot_Enterokinase$hiv_status_recoded <- as.factor(grey60_boxplot_Enterokinase$hiv_status_recoded)

t_test_p23 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = grey60_boxplot_Enterokinase)

t_test_p23 <- t_test_p23 %>% mutate(y.position = 6.6)


ggboxplot(grey60_boxplot_Enterokinase,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "Enterokinase levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p23, label = "p.adj") +
  scale_fill_startrek() 

```

### PM22

```{r message=FALSE, fig.width=2.8, fig.height=3.4}

set.seed(248)

p22 <- scaled_intramodular_connectivity_combined %>% 
  dplyr::filter(module_recoded_PM == "PM22")


hub_prot_22 <- p22 %>% dplyr::filter(AptName == "seq.14143.8")

ggpubr::ggscatter(p22, x = "kWithin", "PS.hiv_status_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.15, label.y = 0.185, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM22", y = "| ß-coefficient | of proteins & HIV ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p22, AptName == "seq.14143.8"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.19, nudge_y = 0.00895,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p22$AptName == "seq.14143.8", "red", "black"), size = 1.8) +
  
         geom_point(data = hub_prot_22,
              aes(x = kWithin, y = PS.hiv_status_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}

lightgreen <- protein_annot_file %>% 
  dplyr::filter(Modules == "lightgreen") 

lightgreen_boxplot <- lightgreen %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


lightgreen_boxplo_H2B1K <- lightgreen_boxplot %>% dplyr::filter(Target ==  "H2B2E")


lightgreen_boxplo_H2B1K$hiv_status_recoded <- factor(lightgreen_boxplo_H2B1K$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

lightgreen_boxplo_H2B1K$hiv_status_recoded <- as.factor(lightgreen_boxplo_H2B1K$hiv_status_recoded)

t_test_p22 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = lightgreen_boxplo_H2B1K)

t_test_p22 <- t_test_p22 %>% mutate(y.position = 3.8)


ggboxplot(lightgreen_boxplo_H2B1K,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "H2B2E levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p22, label = "p.adj") +
  scale_fill_startrek() 

```

##########################################################################################
# Association between HIV-related protein modules and 90-day post-discharge anthropometry
##########################################################################################


#### MUAC

```{r hiv_associated_modules_growth_muac90}

set.seed(20)
# Create an empty list to store the linear models
model_summaries_Adj <- list()

module_cols <- hiv_ncc_sam_growth_analysis[,significantly_associated_hiv_with_baseline_weights_only$Modules]


hiv_ncc_sam_growth_analysis$hiv_group <- factor(hiv_ncc_sam_growth_analysis$hiv_group,
                                                levels = c("HIV (-)", "HIV (+)" ))

discharge_MUAC <-  hiv_ncc_sam_growth_analysis$muac_disch
discharge_age <- hiv_ncc_sam_growth_analysis$discharge_age
sex <- hiv_ncc_sam_growth_analysis$sex_adm
site <- hiv_ncc_sam_growth_analysis$site
hiv_status <- hiv_ncc_sam_growth_analysis$hiv_group

# Iterate over each column of growthreg_data
for (col in colnames(module_cols)) {
  # Create the formula with the variable name
  formula <- as.formula(paste("hiv_ncc_sam_growth_analysis$muac_fu90 ~", col, 
                              "+ discharge_MUAC + discharge_age + sex + hiv_status + (1|site)"))

  # Fit the linear model with the variable, adjusting for confounders
  lm_result1 <- lmerTest::lmer(formula, data = module_cols)

  # Store the linear model in the list, with variable name as the list element name
  model_summaries_Adj[[col]] <- lm_result1
}

 

# Create an empty list to store the coefficient information
coef_estimates <- list()

```


```{r calculate_SEs, warning=FALSE}

for (col in colnames(module_cols)) {
  # Extract the fixed effects coefficients
  fixed_effects <- fixef(model_summaries_Adj[[col]])

  # Compute the standard errors
  SE <- sqrt(diag(vcov(model_summaries_Adj[[col]])))

  # Combine the coefficient and standard error data into a data frame
  coef_df <- data.frame(
    Protein = names(fixed_effects),
    Estimate = fixed_effects,
    SE = SE
  ) 

  # Use broom.mixed to extract p-values
  coef_summary <- broom.mixed::tidy(model_summaries_Adj[[col]])
  p_values <- coef_summary$P.value

  # Add p-values to the coefficient data frame
  coef_df$p_value <- p_values

  # Append the coefficient data frame to the list
  coef_estimates[[col]] <- coef_df
}

 

# Combine the coefficient data frames into a single data frame
adj_coefs_df1 <- do.call(rbind, coef_estimates)

```


```{r exclude_values}

exclude_values <- c("discharge_age", "sexMale", "(Intercept)", "discharge_MUAC")
adj_data <- subset(adj_coefs_df1, !Protein %in% exclude_values)

```

 

###### Calculate the upper and lower error bars

```{r generate_pvalues_CI}

adj_filtered <- adj_data %>%
  group_by(Protein) %>%
  mutate(
    Lower = Estimate - qnorm(0.975) * SE,
    Upper = Estimate + qnorm(0.975) * SE
  ) %>%
  ungroup()
```


###### Adding a significance column

```{r add_sigf_col}

adj_filtered$significance <- ""
adj_filtered[which(adj_filtered$Lower<0 | adj_filtered$Lower>0 & adj_filtered$Upper<0 |adj_filtered$Upper>0),"significance"] <- "Yes"

 

adj_filtered[which(adj_filtered$Lower<0 & adj_filtered$Upper>0),
"significance"] <- "No" 


significantly_associated_muac_90 <- adj_filtered %>% dplyr::filter(significance == "Yes")

```



```{r}

modules_associated_with_growth90_ordered <- adj_filtered %>% 
  arrange(Estimate)



modules_associated_with_growth90_ordered <- modules_associated_with_growth90_ordered %>%  
  dplyr::mutate(Protein = substring(Protein, 3))

modules_associated_with_growth90_ordered <- modules_associated_with_growth90_ordered %>% 
  dplyr::inner_join(., proteinModuleSizes_recoded, by = c("Protein" = "module"))


# level_order <- c(paste(modules_associated_with_growth90_ordered$module_recoded_PM, sep = ","))
# 
# 
modules_associated_with_growth90_ordered$module_recoded_PM <- factor(modules_associated_with_growth90_ordered$module_recoded_PM,
                                                 levels = c("PM1", "PM2", "PM3", "PM5", "PM6", "PM7", "PM8", "PM10", "PM12", 
                                                            "PM14", "PM16", "PM17", "PM18", "PM19", "PM22", "PM23", "PM25", "PM26",
                                                            "PM27",  "PM30", "PM31", "PM32", "PM33", "PM34", "PM37", "PM38", "PM39"))

```



```{r hiv_associated_forest_plt_growth, fig.width=2.1, fig.height=6.5}


modules_associated_with_growth90_plt <- 
  modules_associated_with_growth90_ordered  %>% 
  ggplot(aes(y = fct_rev(module_recoded_PM), col = significance)) + 
  geom_point(aes(x=Estimate), size= 0.45) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.25, linewidth = 0.45) +
  #scale_x_continuous(breaks = seq(-60, 20, by = 10)) +
  #geom_linerange(aes(xmin=Lower_bound_C1_2, xmax=Upper_bound_C1_2))   +
  geom_vline(xintercept = 0, linetype="dashed",  linewidth = 0.35 ) +
  theme_light()+
  labs(x="Estimates (95% CI)", y="Modules", title = "", color = "BP Signf status") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.905),
        text = element_text(size = 9.5),
        axis.text = element_text(face = "bold", size = 9.5),
        legend.title  =  element_text(size=5),
        axis.ticks.length=unit(0.05,"cm"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                colour = "grey95")) +
  ggsci::scale_color_lancet() +
  theme(legend.position = "none")
  
  
modules_associated_with_growth90_plt

```

#### HAZ

```{r hiv_associated_modules_growth_haz_90}

set.seed(291)
# Create an empty list to store the linear models
model_summaries_Adj <- list()

#module_cols <- Proteins_discharge_modules[,2:40]

discharge_HAZ <-  hiv_ncc_sam_growth_analysis$haz_disch
discharge_age <- hiv_ncc_sam_growth_analysis$discharge_age
sex <- hiv_ncc_sam_growth_analysis$sex_adm
site <- hiv_ncc_sam_growth_analysis$site

# Iterate over each column of growthreg_data
for (col in colnames(module_cols)) {
  # Create the formula with the variable name
  formula <- as.formula(paste("hiv_ncc_sam_growth_analysis$haz_fu90 ~", col, 
                              "+ discharge_HAZ + discharge_age + sex + (1|site)"))

  # Fit the linear model with the variable, adjusting for confounders
  lm_result1 <- lmerTest::lmer(formula, data = module_cols)

  # Store the linear model in the list, with variable name as the list element name
  model_summaries_Adj[[col]] <- lm_result1
}

 

# Create an empty list to store the coefficient information
coef_estimates <- list()

```


```{r calculate_SEs, warning=FALSE}

for (col in colnames(module_cols)) {
  # Extract the fixed effects coefficients
  fixed_effects <- fixef(model_summaries_Adj[[col]])

  # Compute the standard errors
  SE <- sqrt(diag(vcov(model_summaries_Adj[[col]])))

  # Combine the coefficient and standard error data into a data frame
  coef_df <- data.frame(
    Protein = names(fixed_effects),
    Estimate = fixed_effects,
    SE = SE
  ) 

  # Use broom.mixed to extract p-values
  coef_summary <- broom.mixed::tidy(model_summaries_Adj[[col]])
  p_values <- coef_summary$P.value

  # Add p-values to the coefficient data frame
  coef_df$p_value <- p_values

  # Append the coefficient data frame to the list
  coef_estimates[[col]] <- coef_df
}

 

# Combine the coefficient data frames into a single data frame
adj_coefs_df1 <- do.call(rbind, coef_estimates)

```


```{r exclude_values}

exclude_values <- c("discharge_age", "sexMale", "(Intercept)", "discharge_HAZ")
adj_data <- subset(adj_coefs_df1, !Protein %in% exclude_values)

```

 
##### Calculate the upper and lower error bars

```{r generate_pvalues_CI}

adj_filtered <- adj_data %>%
  group_by(Protein) %>%
  mutate(
    Lower = Estimate - qnorm(0.975) * SE,
    Upper = Estimate + qnorm(0.975) * SE
  ) %>%
  ungroup()
```


##### Adding a significance column

```{r add_sigf_col}

adj_filtered$significance <- ""
adj_filtered[which(adj_filtered$Lower<0 | adj_filtered$Lower>0 & adj_filtered$Upper<0 |adj_filtered$Upper>0),"significance"] <- "Yes"

 

adj_filtered[which(adj_filtered$Lower<0 & adj_filtered$Upper>0),
"significance"] <- "No" 


significantly_associated_haz_90 <- adj_filtered %>% dplyr::filter(significance == "Yes")

```



```{r}

significantly_associated_haz_90_ordered <- adj_filtered %>% 
  arrange(Estimate)



significantly_associated_haz_90_ordered <- significantly_associated_haz_90_ordered %>%  
  dplyr::mutate(Protein = substring(Protein, 3))

significantly_associated_haz_90_ordered <- significantly_associated_haz_90_ordered %>% 
  dplyr::inner_join(., proteinModuleSizes_recoded, by = c("Protein" = "module"))


#level_order <- c(paste(significantly_associated_haz_90_ordered$module_recoded_PM, sep = ","))


# significantly_associated_haz_90_ordered$module_recoded_PM <- factor(significantly_associated_haz_90_ordered$module_recoded_PM,
#                                                  levels = level_order)

significantly_associated_haz_90_ordered$module_recoded_PM <- factor(significantly_associated_haz_90_ordered$module_recoded_PM,
                                                 levels = c("PM1", "PM2", "PM3", "PM5", "PM6", "PM7", "PM8", "PM10", "PM12", 
                                                            "PM14", "PM16", "PM17", "PM18", "PM19", "PM22", "PM23", "PM25", "PM26",
                                                            "PM27",  "PM30", "PM31", "PM32", "PM33", "PM34", "PM37", "PM38", "PM39"))


```



```{r hiv_associated_forest_plt_growth_HAZ, fig.width=1.9, fig.height=6.5}


modules_associated_with_growth90_plt <- 
  significantly_associated_haz_90_ordered  %>% 
  ggplot(aes(y = fct_rev(module_recoded_PM), col = significance)) + 
  geom_point(aes(x=Estimate), size= 0.45) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.25, linewidth = 0.45) +
  #scale_x_continuous(breaks = seq(-60, 20, by = 10)) +
  #geom_linerange(aes(xmin=Lower_bound_C1_2, xmax=Upper_bound_C1_2))   +
  geom_vline(xintercept = 0, linetype="dashed",  linewidth = 0.35 ) +
  theme_light()+
  labs(x="Estimates (95% CI)", y="", title = "", color = "BP Signf status") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.905),
        text = element_text(size = 9.5),
        axis.text = element_text(face = "bold", size = 9.5),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(.025, "cm"),
        legend.title  =  element_text(size=5),
        axis.ticks.length=unit(0.05,"cm"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                colour = "grey95")) +
  ggsci::scale_color_lancet() +
  theme(legend.position = "none")
  
  
modules_associated_with_growth90_plt

```

#### WAZ

```{r hiv_associated_modules_growth_WAZ90}

set.seed(211)
# Create an empty list to store the linear models
model_summaries_Adj <- list()

#module_cols <- Proteins_discharge_modules[,2:40]

discharge_WAZ <-  hiv_ncc_sam_growth_analysis$waz_disch
discharge_age <- hiv_ncc_sam_growth_analysis$discharge_age
sex <- hiv_ncc_sam_growth_analysis$sex_adm
site <- hiv_ncc_sam_growth_analysis$site

# Iterate over each column of growthreg_data
for (col in colnames(module_cols)) {
  # Create the formula with the variable name
  formula <- as.formula(paste("hiv_ncc_sam_growth_analysis$waz_fu90 ~", col, 
                              "+ discharge_WAZ + discharge_age + sex + (1|site)"))

  # Fit the linear model with the variable, adjusting for confounders
  lm_result1 <- lmerTest::lmer(formula, data = module_cols)

  # Store the linear model in the list, with variable name as the list element name
  model_summaries_Adj[[col]] <- lm_result1
}

 

# Create an empty list to store the coefficient information
coef_estimates <- list()

```



```{r calculate_SEs, warning=FALSE}

for (col in colnames(module_cols)) {
  # Extract the fixed effects coefficients
  fixed_effects <- fixef(model_summaries_Adj[[col]])

  # Compute the standard errors
  SE <- sqrt(diag(vcov(model_summaries_Adj[[col]])))

  # Combine the coefficient and standard error data into a data frame
  coef_df <- data.frame(
    Protein = names(fixed_effects),
    Estimate = fixed_effects,
    SE = SE
  ) 

  # Use broom.mixed to extract p-values
  coef_summary <- broom.mixed::tidy(model_summaries_Adj[[col]])
  p_values <- coef_summary$P.value

  # Add p-values to the coefficient data frame
  coef_df$p_value <- p_values

  # Append the coefficient data frame to the list
  coef_estimates[[col]] <- coef_df
}

 

# Combine the coefficient data frames into a single data frame
adj_coefs_df1 <- do.call(rbind, coef_estimates)

```


```{r exclude_values}

exclude_values <- c("discharge_age", "sexMale", "(Intercept)", "discharge_WAZ")
adj_data <- subset(adj_coefs_df1, !Protein %in% exclude_values)

```


##### Calculate the upper and lower error bars

```{r generate_pvalues_CI}

adj_filtered <- adj_data %>%
  group_by(Protein) %>%
  mutate(
    Lower = Estimate - qnorm(0.975) * SE,
    Upper = Estimate + qnorm(0.975) * SE
  ) %>%
  ungroup()
```


##### Adding a significance column

```{r add_sigf_col}

adj_filtered$significance <- ""
adj_filtered[which(adj_filtered$Lower<0 | adj_filtered$Lower>0 & adj_filtered$Upper<0 |adj_filtered$Upper>0),"significance"] <- "Yes"

 

adj_filtered[which(adj_filtered$Lower<0 & adj_filtered$Upper>0),
"significance"] <- "No" 


significantly_associated_waz_90 <- adj_filtered %>% dplyr::filter(significance == "Yes")

```



```{r}

significantly_associated_waz_90_ordered <- adj_filtered %>% 
  arrange(Estimate)



significantly_associated_waz_90_ordered <- significantly_associated_waz_90_ordered %>%  
  dplyr::mutate(Protein = substring(Protein, 3))

significantly_associated_waz_90_ordered <- significantly_associated_waz_90_ordered %>% 
  dplyr::inner_join(., proteinModuleSizes_recoded, by = c("Protein" = "module"))


# level_order <- c(paste(significantly_associated_waz_90_ordered$module_recoded_PM, sep = ","))
# 
# 
# significantly_associated_waz_90_ordered$module_recoded_PM <- factor(significantly_associated_waz_90_ordered$module_recoded_PM,
#                                                  levels = level_order)

significantly_associated_waz_90_ordered$module_recoded_PM <- factor(significantly_associated_waz_90_ordered$module_recoded_PM,
                                                 levels = c("PM1", "PM2", "PM3", "PM5", "PM6", "PM7", "PM8", "PM10", "PM12", 
                                                            "PM14", "PM16", "PM17", "PM18", "PM19", "PM22", "PM23", "PM25", "PM26",
                                                            "PM27",  "PM30", "PM31", "PM32", "PM33", "PM34", "PM37", "PM38", "PM39"))

```



```{r hiv_associated_forest_plt_growth_HAZ, fig.width=1.8, fig.height=6.5}


modules_associated_with_growth90_plt <- 
  significantly_associated_waz_90_ordered  %>% 
  ggplot(aes(y = fct_rev(module_recoded_PM), col = significance)) + 
  geom_point(aes(x=Estimate), size= 0.45) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.25, linewidth = 0.45) +
  #scale_x_continuous(breaks = seq(-60, 20, by = 10)) +
  #geom_linerange(aes(xmin=Lower_bound_C1_2, xmax=Upper_bound_C1_2))   +
  geom_vline(xintercept = 0, linetype="dashed",  linewidth = 0.35 ) +
  theme_light()+
  labs(x="Estimates (95% CI)", y="", title = "", color = "BP Signf status") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.905),
        text = element_text(size = 9.5),
        axis.text = element_text(face = "bold", size = 9.5),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(.025, "cm"),
        legend.title  =  element_text(size=5),
        axis.ticks.length=unit(0.05,"cm"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                colour = "grey95")) +
  ggsci::scale_color_lancet() +
  theme(legend.position = "none")
  
  
modules_associated_with_growth90_plt

```

#### WHZ

```{r hiv_associated_modules_growth}

set.seed(213)
# Create an empty list to store the linear models
model_summaries_Adj <- list()

#module_cols <- Proteins_discharge_modules[,2:40]

discharge_WHZ <-  hiv_ncc_sam_growth_analysis$whz_disch
discharge_age <- hiv_ncc_sam_growth_analysis$discharge_age
sex <- hiv_ncc_sam_growth_analysis$sex_adm
site <- hiv_ncc_sam_growth_analysis$site

# Iterate over each column of growthreg_data
for (col in colnames(module_cols)) {
  # Create the formula with the variable name
  formula <- as.formula(paste("hiv_ncc_sam_growth_analysis$whz_fu90 ~", col, 
                              "+ discharge_WHZ + discharge_age + sex + (1|site)"))

  # Fit the linear model with the variable, adjusting for confounders
  lm_result1 <- lmerTest::lmer(formula, data = module_cols)

  # Store the linear model in the list, with variable name as the list element name
  model_summaries_Adj[[col]] <- lm_result1
}

 

# Create an empty list to store the coefficient information
coef_estimates <- list()

```

```{r calculate_SEs, warning=FALSE}

for (col in colnames(module_cols)) {
  # Extract the fixed effects coefficients
  fixed_effects <- fixef(model_summaries_Adj[[col]])

  # Compute the standard errors
  SE <- sqrt(diag(vcov(model_summaries_Adj[[col]])))

  # Combine the coefficient and standard error data into a data frame
  coef_df <- data.frame(
    Protein = names(fixed_effects),
    Estimate = fixed_effects,
    SE = SE
  ) 

  # Use broom.mixed to extract p-values
  coef_summary <- broom.mixed::tidy(model_summaries_Adj[[col]])
  p_values <- coef_summary$P.value

  # Add p-values to the coefficient data frame
  coef_df$p_value <- p_values

  # Append the coefficient data frame to the list
  coef_estimates[[col]] <- coef_df
}

 

# Combine the coefficient data frames into a single data frame
adj_coefs_df1 <- do.call(rbind, coef_estimates)

```


```{r exclude_values}

exclude_values <- c("discharge_age", "sexMale", "(Intercept)", "discharge_WHZ")
adj_data <- subset(adj_coefs_df1, !Protein %in% exclude_values)

```

##### Calculate the upper and lower error bars

```{r generate_pvalues_CI}

adj_filtered <- adj_data %>%
  group_by(Protein) %>%
  mutate(
    Lower = Estimate - qnorm(0.975) * SE,
    Upper = Estimate + qnorm(0.975) * SE
  ) %>%
  ungroup()
```


##### Adding a significance column

```{r add_sigf_col}

adj_filtered$significance <- ""
adj_filtered[which(adj_filtered$Lower<0 | adj_filtered$Lower>0 & adj_filtered$Upper<0 |adj_filtered$Upper>0),"significance"] <- "Yes"

 

adj_filtered[which(adj_filtered$Lower<0 & adj_filtered$Upper>0),
"significance"] <- "No" 


significantly_associated_whz_90 <- adj_filtered %>% dplyr::filter(significance == "Yes")

```


```{r}

significantly_associated_whz_90_ordered <- adj_filtered %>% 
  arrange(Estimate)



significantly_associated_whz_90_ordered <- significantly_associated_whz_90_ordered %>%  
  dplyr::mutate(Protein = substring(Protein, 3))

significantly_associated_whz_90_ordered <- significantly_associated_whz_90_ordered %>% 
  dplyr::inner_join(., proteinModuleSizes_recoded, by = c("Protein" = "module"))


# level_order <- c(paste(significantly_associated_whz_90_ordered$module_recoded_PM, sep = ","))
# 
# 
# significantly_associated_whz_90_ordered$module_recoded_PM <- factor(significantly_associated_whz_90_ordered$module_recoded_PM,
#                                                  levels = level_order)


significantly_associated_whz_90_ordered$module_recoded_PM <- factor(significantly_associated_whz_90_ordered$module_recoded_PM,
                                                 levels = c("PM1", "PM2", "PM3", "PM5", "PM6", "PM7", "PM8", "PM10", "PM12", 
                                                            "PM14", "PM16", "PM17", "PM18", "PM19", "PM22", "PM23", "PM25", "PM26",
                                                            "PM27",  "PM30", "PM31", "PM32", "PM33", "PM34", "PM37", "PM38", "PM39"))

```



```{r hiv_associated_forest_plt_growth_WHZ, fig.width=1.9, fig.height=6.5}


modules_associated_with_growth90_plt <- 
  significantly_associated_whz_90_ordered  %>% 
  ggplot(aes(y = fct_rev(module_recoded_PM), col = significance)) + 
  geom_point(aes(x=Estimate), size= 0.45) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), 
                 height = 0.25, linewidth = 0.45) +
  #scale_x_continuous(breaks = seq(-60, 20, by = 10)) +
  #geom_linerange(aes(xmin=Lower_bound_C1_2, xmax=Upper_bound_C1_2))   +
  geom_vline(xintercept = 0, linetype="dashed",  linewidth = 0.35 ) +
  theme_light()+
  labs(x="Estimates (95% CI)", y="", title = "", color = "BP Signf status") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.905),
        text = element_text(size = 9.5),
        axis.text = element_text(face = "bold", size = 9.5),
        legend.title  =  element_text(size=5),
        axis.text.y = element_blank(),
        axis.ticks.length.y = unit(.025, "cm"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "white"),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                colour = "grey95")) +
  ggsci::scale_color_lancet() +
  theme(legend.position = "none")
  
  
modules_associated_with_growth90_plt

```

##########################################################################################
	# Characteristics of the participants with plasma proteome data at discharge, who were included in the HIV mechanistic analysis
##########################################################################################


```{r}

hiv_ncc_sam_growth_analysis_ids <- hiv_ncc_sam_growth_analysis %>% dplyr::select(record_id)

hiv_ncc_sam_growth_analysis_ids_table <- hiv_ncc_sam_growth_analysis_ids %>% 
  dplyr::inner_join(., HIV_ncc_growth_Table1_recodIDG, by = "record_id")


hiv_ncc_sam_growth_analysis_ids_tableIDs <- hiv_ncc_sam_growth_analysis_ids_table %>% 
  dplyr::select(-c(record_id)) #resp_diag_pulmonarytb_disch

hiv_ncc_sam_growth_analysis_ids_tableIDs %>% 
  tbl_summary(by = hiv_status,
              type = list(#age_adm ~ 'continuous2',
                          discharge_age ~ 'continuous2',
                          muac_disch ~ 'continuous2',
                          haz_disch ~ 'continuous2',
                          waz_disch ~ 'continuous2',
                          whz_disch ~ 'continuous2',
                          head_circ_disch ~ 'continuous2'),
              statistic = list(all_continuous() ~ c("{median} ({p25} - {p75})")),
              missing =  "no",
              label = list(sex_adm ~ "Males (%)",
                           #age_adm ~ "Admission age, months",
                           discharge_age ~ " Discharge age, months",
                           site ~ "Site",
                           nutritional_category ~ "Nutritional status at discharge",
                           diarrhoea ~ "Diarrhoea",
                           muac_disch ~ "MUAC (cm)",
                           haz_disch ~ "HAZ score",
                           waz_disch ~ "WAZ score",
                           whz_disch ~ "WHZ score",
                           head_circ_disch ~ "Head circumference",
                           oedema_disch ~ "Oedema - Yes",
                           infec_diag_conrmedmalaria_disch_processed ~  "Malaria Positive (RDT)",
                           infec_diag_measles_disch_processed ~ "Measles", 
                           infec_diag_sepsis_disch_processed ~ "Sepsis",
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
  writexl::write_xlsx(., "../Proteomics/Somalogic-2/Significant_modules/NCC_HIV_substudy/Tables/Growth_study_participants_xtics.xlsx")



```


######################################################################################################
	# Intramodular connectivity of HIV-association modules with anthropometry at day 90 post-discharge
######################################################################################################


## MUAC

### PM12

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(240)

p12 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM12")


hub_prot_12_muac <- p12 %>% dplyr::filter(Target == "ITM2A")

ggpubr::ggscatter(p12, x = "kWithin", "P.S.muac_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.35, label.y = 0.115, size = 5.5,
           label.sep='\n') +
  labs(x = "Connectivity within PM12", y = "| ß-coefficient | of proteins & MUAC ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p12, Target == "ITM2A"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = -0.00225,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p12$Target == "ITM2A", "red", "black"), size = 1.8) +
    
         geom_point(data = hub_prot_12_muac,
              aes(x = kWithin, y = P.S.muac_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.65, fig.height=3.6, warning=FALSE}


skyblue <- protein_annot_file %>% 
  dplyr::filter(Modules == "skyblue") 

skyblue_boxplot <- skyblue %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


skyblue_boxplot_MOB3B <- skyblue_boxplot %>% dplyr::filter(Target ==  "ITM2A")


skyblue_boxplot_MOB3B$hiv_status_recoded <- factor(skyblue_boxplot_MOB3B$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

skyblue_boxplot_MOB3B$hiv_status_recoded <- as.factor(skyblue_boxplot_MOB3B$hiv_status_recoded)

t_test_p12 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = skyblue_boxplot_MOB3B)

t_test_p12 <- t_test_p12 %>% mutate(y.position = 6.8)


ggboxplot(skyblue_boxplot_MOB3B,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "ITM2A levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p12, label = "p.adj") +
  scale_fill_startrek() 

```

### PM31

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(281)

p31 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM31")

p31 <- p31 %>% dplyr::mutate(Target = dplyr::recode(Target, "cAMP-regulated phosphoprotein 21" = "ARPP21"))

hub_prot_31_muac <- p31 %>% dplyr::filter(Target == "ARPP21")

ggpubr::ggscatter(p31, x = "kWithin", "P.S.muac_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.45, label.y = 0.13, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM31", y = "| ß-coefficient | of proteins & MUAC ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p31, Target == "ARPP21"),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = -0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 # geom_point(color = ifelse(p31$Target == "ARPP21", "red", "black"), size = 1.5) +
         geom_point(data = hub_prot_31_muac,
              aes(x = kWithin, y = P.S.muac_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}


magenta <- protein_annot_file %>% 
  dplyr::filter(Modules == "magenta") 

magenta_boxplot <- magenta %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


magenta_boxplot_LRRT4 <- magenta_boxplot %>% dplyr::filter(Target == "cAMP-regulated phosphoprotein 21")

magenta_boxplot_LRRT4$hiv_status_recoded <- factor(magenta_boxplot_LRRT4$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

magenta_boxplot_LRRT4$hiv_status_recoded <- as.factor(magenta_boxplot_LRRT4$hiv_status_recoded)

t_test_p31 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = magenta_boxplot_LRRT4)

t_test_p31 <- t_test_p31 %>% mutate(y.position = 11.4)


ggboxplot(magenta_boxplot_LRRT4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "ARPP21 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p31, label = "p.adj") +
  scale_fill_startrek() 

```

### PM37

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(230)

brown_ki_cmbd_protein <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM37")

hub_prot_37_muac <- brown_ki_cmbd_protein %>% dplyr::filter(Target == "NGRN")

ggpubr::ggscatter(brown_ki_cmbd_protein, x = "kWithin", "P.S.muac_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.45, label.y = 0.165, size = 5,
           label.sep='\n') +
  labs(x = "Connectivity within PM37", y = "| ß-coefficient | of proteins & MUAC ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(brown_ki_cmbd_protein, P.S.muac_fu90_abs > 0.1 & kWithin > 0.95),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
  #geom_point(color = ifelse(brown_ki_cmbd_protein$Target == "NGRN", "red", "black"), size = 1.5) +
  
    geom_point(data = hub_prot_37_muac,
              aes(x = kWithin, y = P.S.muac_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


brown <- protein_annot_file %>% 
  dplyr::filter(Modules == "brown") 

brown_boxplot <- brown %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


brown_boxplot_SPEF1 <- brown_boxplot %>% dplyr::filter(Target == "NGRN")


brown_boxplot_SPEF1$hiv_status_recoded <- factor(brown_boxplot_SPEF1$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

brown_boxplot_SPEF1$hiv_status_recoded <- as.factor(brown_boxplot_SPEF1$hiv_status_recoded)

t_test_p37 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = brown_boxplot_SPEF1)

t_test_p37 <- t_test_p37 %>% mutate(y.position = 10.1)


ggboxplot(brown_boxplot_SPEF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "NGRN levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p37, label = "p.adj") +
  scale_fill_startrek() 

```

### PM26

```{r message=FALSE, fig.width=2.6, fig.height=3.5}

set.seed(246)

p26 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM26")


hub_prot_26_muac <- p26 %>% dplyr::filter(Target == "CO9A1")

ggpubr::ggscatter(p26, x = "kWithin", "P.S.muac_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.05, label.y = 0.175, size = 5,
           label.sep=',') +
  #ylim(0,0.16) +
  labs(x = "Connectivity within PM26", y = "| ß-coefficient | of proteins & MUAC ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p26, P.S.muac_fu90_abs < 0.012 & kWithin > 0.95),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = -0.0125,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p26$Target == "CO9A1", "red", "black"), size = 1.8) +
  
     geom_point(data = hub_prot_26_muac,
              aes(x = kWithin, y = P.S.muac_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.65, fig.height=3.6, warning=FALSE}

cyan <- protein_annot_file %>% 
  dplyr::filter(Modules == "cyan") 

cyan_boxplot <- cyan %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


cyan_boxplot_COL11A2 <- cyan_boxplot %>% dplyr::filter(Target ==  "CO9A1")


cyan_boxplot_COL11A2$hiv_status_recoded <- factor(cyan_boxplot_COL11A2$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

cyan_boxplot_COL11A2$hiv_status_recoded <- as.factor(cyan_boxplot_COL11A2$hiv_status_recoded)

t_test_p26 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = cyan_boxplot_COL11A2)

t_test_p26 <- t_test_p26 %>% mutate(y.position = 2.95)


ggboxplot(cyan_boxplot_COL11A2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "CO9A1 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p26, label = "p.adj") +
  scale_fill_startrek() 

```

## HAZ


### PM3

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(243)

p3 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM3")


p3 <- p3 %>% dplyr::mutate(Target = dplyr::recode(Target, "HSP 70" = "HSP70"))

hub_prot_3_haz <- p3 %>% dplyr::filter(Target == "TPGS2")

ggpubr::ggscatter(p3, x = "kWithin", "P.S.haz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.04, label.y = 0.195, size = 5,
           label.sep=',') + 
  labs(x = "Connectivity within PM3", y = "| ß-coefficient | of proteins & HAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p3, P.S.haz_fu90_abs > 0.2 & kWithin > 0.9),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = 0.00495,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p3$Target == "TPGS2", "red", "black"), size = 1.8) =
    
     geom_point(data = hub_prot_3_haz,
              aes(x = kWithin, y = P.S.haz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


skyblue3 <- protein_annot_file %>% 
  dplyr::filter(Modules == "skyblue3") 

skyblue3_boxplot <- skyblue3 %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


skyblue3_boxplot_KNG2 <- skyblue3_boxplot %>% dplyr::filter(Target ==  "TPGS2")


skyblue3_boxplot_KNG2$hiv_status_recoded <- factor(skyblue3_boxplot_KNG2$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

skyblue3_boxplot_KNG2$hiv_status_recoded <- as.factor(skyblue3_boxplot_KNG2$hiv_status_recoded)

t_test_p3 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = skyblue3_boxplot_KNG2)

t_test_p3<- t_test_p3 %>% mutate(y.position = 3.4)


ggboxplot(skyblue3_boxplot_KNG2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "TPGS2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p3, label = "p.adj") +
  scale_fill_startrek() 

```

### PM12

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(240)

p12 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM12")

hub_prot_12_haz <- p12 %>% dplyr::filter(Target == "MOB3B")

ggpubr::ggscatter(p12, x = "kWithin", "P.S.haz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.2, label.y = 0.155, size = 4,
           label.sep=',') +
  labs(x = "Connectivity within PM12", y = "| ß-coefficient | of proteins & HAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p12, P.S.haz_fu90_abs > 0.13 & kWithin > 0.9),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = -0.035, nudge_y = -0.0105,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p12$Target == "MOB3B", "red", "black"), size = 1.8) +
  
       geom_point(data = hub_prot_12_haz,
              aes(x = kWithin, y = P.S.haz_fu90_abs),
              color = "red", size = 1.25)
              
```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}


skyblue <- protein_annot_file %>% 
  dplyr::filter(Modules == "skyblue") 

skyblue_boxplot <- skyblue %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


skyblue_boxplot_MOB3B <- skyblue_boxplot %>% dplyr::filter(Target ==  "MOB3B")


skyblue_boxplot_MOB3B$hiv_status_recoded <- factor(skyblue_boxplot_MOB3B$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

skyblue_boxplot_MOB3B$hiv_status_recoded <- as.factor(skyblue_boxplot_MOB3B$hiv_status_recoded)

t_test_p12 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = skyblue_boxplot_MOB3B)

t_test_p12 <- t_test_p12 %>% mutate(y.position = 4.2)


ggboxplot(skyblue_boxplot_MOB3B,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "MOB3B levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p12, label = "p.adj") +
  scale_fill_startrek() 

```

### PM18

```{r message=FALSE, fig.width=2.6, fig.height=3.2}

set.seed(241)

p18 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM18")

p18 <- p18 %>% dplyr::filter(Target != "No Protein")


p18 <- p18 %>% dplyr::mutate(Target = dplyr::recode(Target, "Carbonic anhydrase VII" = "CA7"))

hub_prot_18_haz <- p18 %>% dplyr::filter(Target == "LRRT3")

ggpubr::ggscatter(p18, x = "kWithin", "P.S.haz_fu90_abs",
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.2, label.y = 0.11, size = 5,
           label.sep=',') +
  labs(x = "Connectivity within PM18", y = "| ß-coefficient | of proteins & HAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p18, P.S.haz_fu90_abs > 0.06 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = 0.005,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p18$Target == "LRRT3", "red", "black"), size = 1.8) +
  
      geom_point(data = hub_prot_18_haz,
              aes(x = kWithin, y = P.S.haz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}


darkgreen <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkgreen") 

darkgreen_boxplot <- darkgreen %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkgreen_boxplot_CA7 <- darkgreen_boxplot %>% dplyr::filter(Target ==  "LRRT3")


darkgreen_boxplot_CA7$hiv_status_recoded <- factor(darkgreen_boxplot_CA7$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkgreen_boxplot_CA7$hiv_status_recoded <- as.factor(darkgreen_boxplot_CA7$hiv_status_recoded)

t_test_p18 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkgreen_boxplot_CA7)

t_test_p18 <- t_test_p18 %>% mutate(y.position = 9.2)


ggboxplot(darkgreen_boxplot_CA7,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "LRRT3 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p18, label = "p.adj") +
  scale_fill_startrek() 

```

## WAZ


### PM31

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(281)

p31 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM31")

p31 <- p31 %>% dplyr::mutate(Target = dplyr::recode(Target, "cAMP-regulated phosphoprotein 21" = "ARPP21"))

hub_prot_31_waz <- p31 %>% dplyr::filter(Target == "SMDC1")

ggpubr::ggscatter(p31, x = "kWithin", "P.S.waz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.2, label.y = 0.12, size = 4.5,
           label.sep=',') +
  labs(x = "Connectivity within PM31", y = "| ß-coefficient | of proteins & WAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p31, P.S.waz_fu90_abs < 0.01 & kWithin > 0.7),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = -0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  #geom_point(color = ifelse(p31$Target == "SMDC1", "red", "black"), size = 1.5) +
  
      geom_point(data = hub_prot_31_waz,
              aes(x = kWithin, y = P.S.waz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


magenta <- protein_annot_file %>% 
  dplyr::filter(Modules == "magenta") 

magenta_boxplot <- magenta %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


magenta_boxplot_LRRT4 <- magenta_boxplot %>% dplyr::filter(Target == "SMDC1")


magenta_boxplot_LRRT4$hiv_status_recoded <- factor(magenta_boxplot_LRRT4$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

magenta_boxplot_LRRT4$hiv_status_recoded <- as.factor(magenta_boxplot_LRRT4$hiv_status_recoded)

t_test_p31 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = magenta_boxplot_LRRT4)

t_test_p31 <- t_test_p31 %>% mutate(y.position = 11.4)


ggboxplot(magenta_boxplot_LRRT4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "SMDC1 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p31, label = "p.adj") +
  scale_fill_startrek() 

```


### PM37

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(230)

brown_ki_cmbd_protein <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM37")

hub_prot_37_waz <- brown_ki_cmbd_protein %>% dplyr::filter(Target == "ISL1")

ggpubr::ggscatter(brown_ki_cmbd_protein, x = "kWithin", "P.S.waz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.2, label.y = 0.155, size = 5,
           label.sep=',') +
  labs(x = "Connectivity within PM37", y = "| ß-coefficient | of proteins & WAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(brown_ki_cmbd_protein, P.S.waz_fu90_abs > 0.1 & kWithin > 0.85),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
  #geom_point(color = ifelse(brown_ki_cmbd_protein$Target == "ISL1", "red", "black"), size = 1.5) +
  
        geom_point(data = hub_prot_37_waz,
              aes(x = kWithin, y = P.S.waz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


brown <- protein_annot_file %>% 
  dplyr::filter(Modules == "brown") 

brown_boxplot <- brown %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


brown_boxplot_SPEF1 <- brown_boxplot %>% dplyr::filter(Target == "ISL1")


brown_boxplot_SPEF1$hiv_status_recoded <- factor(brown_boxplot_SPEF1$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

brown_boxplot_SPEF1$hiv_status_recoded <- as.factor(brown_boxplot_SPEF1$hiv_status_recoded)

t_test_p37 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = brown_boxplot_SPEF1)

t_test_p37 <- t_test_p37 %>% mutate(y.position = 8.6)


ggboxplot(brown_boxplot_SPEF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "ISL1 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p37, label = "p.adj") +
  scale_fill_startrek() 

```

### PM26

```{r message=FALSE, fig.width=2.6, fig.height=3.5}

set.seed(246)

p26 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM26")


hub_prot_26_waz <- p26 %>% dplyr::filter(Target == "CILP2")

ggpubr::ggscatter(p26, x = "kWithin", "P.S.waz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +

  stat_cor(method  = "pearson", label.x = 0.05, label.y = 0.23, size = 5,
           label.sep=',') +
  #ylim(0,0.16) +
  labs(x = "Connectivity within PM26", y = "| ß-coefficient | of proteins & WAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p26, P.S.waz_fu90_abs > 0.15 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.1, nudge_y = -0.0125,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p26$Target == "CILP2", "red", "black"), size = 1.8) +
  
          geom_point(data = hub_prot_26_waz,
              aes(x = kWithin, y = P.S.waz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}

cyan <- protein_annot_file %>% 
  dplyr::filter(Modules == "cyan") 

cyan_boxplot <- cyan %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


cyan_boxplot_COL11A2 <- cyan_boxplot %>% dplyr::filter(Target ==  "CILP2")

cyan_boxplot_COL11A2$hiv_status_recoded <- factor(cyan_boxplot_COL11A2$hiv_status_recoded,
                                                   levels = c("HIV+","HIV-"))


my_comparisons <- list( c("HIV-", "HIV+"))

cyan_boxplot_COL11A2$hiv_status_recoded <- as.factor(cyan_boxplot_COL11A2$hiv_status_recoded)

t_test_p26 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = cyan_boxplot_COL11A2)

t_test_p26 <- t_test_p26 %>% mutate(y.position = 2.5)


ggboxplot(cyan_boxplot_COL11A2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "CILP2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p26, label = "p.adj") +
  scale_fill_startrek() 

```

### PM7

```{r message=FALSE, fig.width=2.5, fig.height=3.2}

set.seed(244)

p7 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM7")


hub_prot_7_waz <- p7 %>% dplyr::filter(Target == "IGF-I")

ggpubr::ggscatter(p7, x = "kWithin", "P.S.waz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.1, label.y = 0.31, size = 5,
           label.sep=',') +
  labs(x = "Connectivity within PM7", y = "| ß-coefficient | of proteins & WAZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p7, P.S.waz_fu90_abs > 0.25 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.19, nudge_y = 0.0195,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p7$Target == "IGF-I", "red", "black"), size = 1.8) +
    
          geom_point(data = hub_prot_7_waz,
              aes(x = kWithin, y = P.S.waz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}

darkolivegreen <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkolivegreen") 

darkolivegreen_boxplot <- darkolivegreen %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkolivegreen_boxplot_IGF1 <- darkolivegreen_boxplot %>% dplyr::filter(Target ==  "IGF-I")


darkolivegreen_boxplot_IGF1$hiv_status_recoded <- factor(darkolivegreen_boxplot_IGF1$hiv_status_recoded,
                                                         levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkolivegreen_boxplot_IGF1$hiv_status_recoded <- as.factor(darkolivegreen_boxplot_IGF1$hiv_status_recoded)

t_test_p7 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkolivegreen_boxplot_IGF1)

t_test_p7 <- t_test_p7 %>% mutate(y.position = 2.5)


ggboxplot(darkolivegreen_boxplot_IGF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "IGF-I levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p7, label = "p.adj") +
  scale_fill_startrek() 

```

### PM31

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(281)

p31 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM31")

p31 <- p31 %>% dplyr::mutate(Target = dplyr::recode(Target, "cAMP-regulated phosphoprotein 21" = "ARPP21"))

hub_prot_31_whz <- p31 %>% dplyr::filter(Target == "BRD2")

ggpubr::ggscatter(p31, x = "kWithin", "P.S.whz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.2, label.y = 0.13, size = 4.5,
           label.sep=',') +
  labs(x = "Connectivity within PM31", y = "| ß-coefficient | of proteins & WHZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p31, P.S.whz_fu90_abs > 0.05 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = -0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 # geom_point(color = ifelse(p31$Target == "BRD2", "red", "black"), size = 1.5) +
  
          geom_point(data = hub_prot_31_whz,
              aes(x = kWithin, y = P.S.whz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}


magenta <- protein_annot_file %>% 
  dplyr::filter(Modules == "magenta") 

magenta_boxplot <- magenta %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


magenta_boxplot_LRRT4 <- magenta_boxplot %>% dplyr::filter(Target == "BRD2")


magenta_boxplot_LRRT4$hiv_status_recoded <- factor(magenta_boxplot_LRRT4$hiv_status_recoded,
                                                         levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

magenta_boxplot_LRRT4$hiv_status_recoded <- as.factor(magenta_boxplot_LRRT4$hiv_status_recoded)

t_test_p31 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = magenta_boxplot_LRRT4)

t_test_p31 <- t_test_p31 %>% mutate(y.position = 9.9)


ggboxplot(magenta_boxplot_LRRT4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "BRD2 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p31, label = "p.adj") +
  scale_fill_startrek() 

```

### PM37

```{r message=FALSE, fig.width=2.6, fig.height=3.4}

set.seed(230)

brown_ki_cmbd_protein <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM37")

hub_prot_37_whz <- brown_ki_cmbd_protein %>% dplyr::filter(Target == "ISL1")

ggpubr::ggscatter(brown_ki_cmbd_protein, x = "kWithin", "P.S.whz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.2, label.y = 0.155, size = 5,
           label.sep=',') +
  labs(x = "Connectivity within PM37", y = "| ß-coefficient | of proteins & WHZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(brown_ki_cmbd_protein, P.S.whz_fu90_abs > 0.1 & kWithin > 0.85),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.185, nudge_y = 0.0075,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20))+
 # geom_point(color = ifelse(brown_ki_cmbd_protein$Target == "ISL1", "red", "black"), size = 1.5) +
    
          geom_point(data = hub_prot_37_whz,
              aes(x = kWithin, y = P.S.whz_fu90_abs),
              color = "red", size = 1.25)
              
```

#### Box plot

```{r fig.width=1.6, fig.height=3.6, warning=FALSE}


brown <- protein_annot_file %>% 
  dplyr::filter(Modules == "brown") 

brown_boxplot <- brown %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


brown_boxplot_SPEF1 <- brown_boxplot %>% dplyr::filter(Target == "ISL1")


brown_boxplot_SPEF1$hiv_status_recoded <- factor(brown_boxplot_SPEF1$hiv_status_recoded,
                                                         levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

brown_boxplot_SPEF1$hiv_status_recoded <- as.factor(brown_boxplot_SPEF1$hiv_status_recoded)

t_test_p37 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = brown_boxplot_SPEF1)

t_test_p37 <- t_test_p37 %>% mutate(y.position = 8.6)


ggboxplot(brown_boxplot_SPEF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "ISL1 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p37, label = "p.adj") +
  scale_fill_startrek() 

```

### PM26

```{r message=FALSE, fig.width=2.6, fig.height=3.5}

set.seed(246)

p26 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM26")


hub_prot_26_whz <- p26 %>% dplyr::filter(Target == "TSP3")

ggpubr::ggscatter(p26, x = "kWithin", "P.S.whz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  #geom_smooth(se = FALSE, cor.coef = FALSE) +
  stat_cor(method  = "pearson", label.x = 0.05, label.y = 0.2, size = 5,
           label.sep=',') +
  labs(x = "Connectivity within PM26", y = "| ß-coefficient | of proteins & WHZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p26, P.S.whz_fu90_abs > 0.1 & kWithin > 0.85),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.165, nudge_y = -0.0125,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p26$Target == "TSP3", "red", "black"), size = 1.8) +
  
          geom_point(data = hub_prot_26_whz,
              aes(x = kWithin, y = P.S.whz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}

cyan <- protein_annot_file %>% 
  dplyr::filter(Modules == "cyan") 

cyan_boxplot <- cyan %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


cyan_boxplot_COL11A2 <- cyan_boxplot %>% dplyr::filter(Target ==  "TSP3")

cyan_boxplot_COL11A2$hiv_status_recoded <- factor(cyan_boxplot_COL11A2$hiv_status_recoded,
                                                         levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

cyan_boxplot_COL11A2$hiv_status_recoded <- as.factor(cyan_boxplot_COL11A2$hiv_status_recoded)

t_test_p26 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = cyan_boxplot_COL11A2)

t_test_p26 <- t_test_p26 %>% mutate(y.position = 3.2)


ggboxplot(cyan_boxplot_COL11A2,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "TSP3 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p26, label = "p.adj") +
  scale_fill_startrek() 

```

### PM7

```{r message=FALSE, fig.width=2.5, fig.height=3.2}

set.seed(244)

p7 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM7")


hub_prot_7_whz <- p7 %>% dplyr::filter(Target == "IGF-I")

ggpubr::ggscatter(p7, x = "kWithin", "P.S.whz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.1, label.y = 0.21, size = 5,
           label.sep=',') +
  labs(x = "Connectivity within PM7", y = "| ß-coefficient | of proteins & WHZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p7, P.S.whz_fu90_abs > 0.12 & kWithin > 0.8),
                           aes(label = Target),
                           size = 4,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.19, nudge_y = 0.0155,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
 #geom_point(color = ifelse(p7$Target == "IGF-I", "red", "black"), size = 1.8) +
  
          geom_point(data = hub_prot_7_whz,
              aes(x = kWithin, y = P.S.whz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.5, fig.height=3.6, warning=FALSE}

darkolivegreen <- protein_annot_file %>% 
  dplyr::filter(Modules == "darkolivegreen") 

darkolivegreen_boxplot <- darkolivegreen %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


darkolivegreen_boxplot_IGF1 <- darkolivegreen_boxplot %>% dplyr::filter(Target ==  "IGF-I")

darkolivegreen_boxplot_IGF1$hiv_status_recoded <- factor(darkolivegreen_boxplot_IGF1$hiv_status_recoded,
                                                         levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

darkolivegreen_boxplot_IGF1$hiv_status_recoded <- as.factor(darkolivegreen_boxplot_IGF1$hiv_status_recoded)

t_test_p7 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = darkolivegreen_boxplot_IGF1)

t_test_p7 <- t_test_p7 %>% mutate(y.position = 2.5)


ggboxplot(darkolivegreen_boxplot_IGF1,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "IGF-I levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p7, label = "p.adj") +
  scale_fill_startrek() 

```

### PM33

```{r message=FALSE, fig.width=2.8, fig.height=3.5}

set.seed(244)

p33 <- scaled_intramodular_connectivity_combined_growth %>% 
  dplyr::filter(module_recoded_PM == "PM33")


hub_prot_33_whz <- p33 %>% dplyr::filter(Target == "STMN4")

ggpubr::ggscatter(p33, x = "kWithin", "P.S.whz_fu90_abs", 
                  conf.int = TRUE, add = "reg.line", size = 0.95,
                  #label = "Target", 
                  repel = FALSE,
                  add.params = list(color = "blue", fill = "lightgray", size = 0.65)) +
  stat_cor(method  = "pearson", label.x = 0.1, label.y = 0.15, size = 4,
           label.sep=',') +
  labs(x = "Connectivity within PM33", y = "| ß-coefficient | of proteins & WHZ ") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(face = "plain")) +
    ggrepel::geom_text_repel(data=subset(p33, P.S.whz_fu90_abs > 0.09 & kWithin > 0.9),
                           aes(label = Target),
                           size = 3,
                           #arrow = arrow(length = unit(0.025, "inches"), type = "open", angle = 30, ends="last"),
                           nudge_x = 0.05, nudge_y = 0.0115,
                           box.padding = 0.05, seed = 1934, parse = TRUE,
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
# geom_point(color = ifelse(p33$Target == "STMN4", "red", "black"), size = 1.8) +
    
          geom_point(data = hub_prot_33_whz,
              aes(x = kWithin, y = P.S.whz_fu90_abs),
              color = "red", size = 1.25)

```

#### Box plot

```{r fig.width=1.7, fig.height=3.6, warning=FALSE}

black <- protein_annot_file %>% 
  dplyr::filter(Modules == "black") 

black_boxplot <- black %>% dplyr::inner_join(., log_transformed_autoscaled_raw_protein_D0_long, 
                                                             by = "AptName")


black_boxplot_STMN4 <- black_boxplot %>% dplyr::filter(Target ==  "STMN4")

black_boxplot_STMN4$hiv_status_recoded <- factor(black_boxplot_STMN4$hiv_status_recoded,
                                                         levels = c("HIV+","HIV-"))

my_comparisons <- list( c("HIV-", "HIV+"))

black_boxplot_STMN4$hiv_status_recoded <- as.factor(black_boxplot_STMN4$hiv_status_recoded)

t_test_p33 <- compare_means(Protein_levels ~ hiv_status_recoded, 
                            comparisons = my_comparisons, 
                            p.adjust.method = "bonferroni", 
                            method='t.test', 
                            data = black_boxplot_STMN4)

t_test_p33 <- t_test_p33 %>% mutate(y.position = 4.0)


ggboxplot(black_boxplot_STMN4,
          x = "hiv_status_recoded", y = "Protein_levels", 
          fill = "hiv_status_recoded", bxp.errorbar.width = 0.2,
          bxp.errorbar = TRUE, size = 0.25) +
  labs(x = NULL, y =  "STMN4 levels (log)") +
  theme(legend.position = "na",
        text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "plain")) +
  stat_pvalue_manual(t_test_p33, label = "p.adj") +
  scale_fill_startrek() 

```

##########################################################################################
	# Structural equation modelling
##########################################################################################


## MUAC SEM MODELS

```{r all_muac_sem, fig.width=4, fig.height=8, warning=TRUE}

set.seed(147)

all_muac_sem_1 <- '
# regression 

MEskyblue ~  hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEcyan ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEmagenta ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEbrown ~ hiv_status_recored + discharge_age + sex_recoded + site_recoded 


muac_disch ~ MEskyblue + MEcyan + MEmagenta  + MEbrown + sex_recoded + hiv_status_recored + discharge_age #+ site_recoded

muac_fu90 ~ MEskyblue + MEcyan + MEmagenta  + MEbrown + muac_disch + hiv_status_recored + discharge_age + sex_recoded #+site_recoded


# covariance

MEskyblue ~~ MEcyan

MEcyan ~~ MEmagenta

MEcyan ~~ MEbrown

MEskyblue ~~ MEmagenta

MEskyblue ~~ MEbrown

MEmagenta ~~ MEbrown

'
all_muac_sem_fit <- sem(all_muac_sem_1, 
                        missing = "ML",
                        likelihood = "wishart",
                        sample.cov.rescale = TRUE,
                        data = hiv_ncc_sam_growth_analysis_SEM)

summary(all_muac_sem_fit, fit.measures = TRUE, standardized = TRUE)


node_labeles = list(muac_disch = "Baseline MUAC", 
                    muac_fu90 = "Post-discharge MUAC at Day 90", 
                    MEskyblue = "PM12",
                    MEcyan = "PM26", 
                    MEmagenta = "PM31",
                    MEbrown = "PM37",
                    site_recoded = "Site",
                    hiv_status_recored = "HIV", 
                    discharge_age = "Discharge age", 
                    sex_recoded = "Sex"
                    )


all_muac_sem_plt <- lavaanPlot(model = all_muac_sem_fit, coefs = TRUE,
                               stand = TRUE, #sig = 0.05, 
                               labels = node_labeles,
                               covs = TRUE,
                               edge_options = list(color = "grey50"),
                               graph_options = list(rankdir = "TB"),
                               stars = "regress")

all_muac_sem_plt


varTable(all_muac_sem_fit)

modindices(all_muac_sem_fit, sort = TRUE)

#lavTestLRT(all_muac_sem_fit, all_muac_sem_fit1)

```

## WAZ SEM MODEL


```{r all_waz_sem, fig.width=4, fig.height=8, warning=FALSE}

set.seed(151)

all_waz_sem <- '
# regression 

MEdarkolivegreen ~  hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEcyan ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEmagenta ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEbrown ~ hiv_status_recored + discharge_age + sex_recoded + site_recoded 

waz_disch ~ MEdarkolivegreen + MEcyan + MEmagenta  + MEbrown + sex_recoded + hiv_status_recored + discharge_age 

waz_fu90 ~ MEdarkolivegreen + MEcyan + MEmagenta  + MEbrown + waz_disch + hiv_status_recored + discharge_age + sex_recoded


# covariance

MEdarkolivegreen ~~ MEcyan

MEcyan ~~ MEmagenta

MEcyan ~~ MEbrown

MEdarkolivegreen ~~ MEmagenta

MEdarkolivegreen ~~ MEbrown

MEmagenta ~~ MEbrown

'

all_waz_sem_fit <- sem(all_waz_sem,
                        missing = "ML",
                        likelihood = "wishart",
                        sample.cov.rescale = TRUE,
                        data = hiv_ncc_sam_growth_analysis_normalized_WAZ)


summary(all_waz_sem_fit, fit.measures = TRUE, standardized = TRUE)


node_labeles = list(waz_disch = "Baseline WAZ", 
                    waz_fu90 = "Post-discharge WAZ at Day 90", 
                    MEdarkolivegreen = "PM7",
                    MEcyan = "PM26", 
                    MEmagenta = "PM31",
                    MEbrown = "PM37",
                    #MEyellow = "PM36",
                    site_recoded = "Site",
                    #malaria = "Malaria", diarrhoea = "Diarrhoea", severe_pneumonia = "Pneumonia",
                    hiv_status_recored = "HIV", 
                    discharge_age = "Discharge age", 
                    sex_recoded = "Sex"
                    )


all_waz_sem_plt <- lavaanPlot(model = all_waz_sem_fit, coefs = TRUE,
                               stand = TRUE, #sig = 0.05, 
                               labels = node_labeles,
                               covs = TRUE,
                               edge_options = list(color = "grey50"),
                               graph_options = list(rankdir = "TB"),
                               stars = "regress")


all_waz_sem_plt


varTable(all_waz_sem_fit)

modindices(all_waz_sem_fit, sort = TRUE)


```

## WHZ SEM MODELS


```{r all_whz_sem, fig.width=4, fig.height=8, warning=FALSE}

set.seed(153)

all_whz_sem <- '
# regression 

MEdarkolivegreen ~  hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEcyan ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEmagenta ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEblack ~  hiv_status_recored + discharge_age + sex_recoded + site_recoded 

MEbrown ~ hiv_status_recored + discharge_age + sex_recoded + site_recoded 


# covariance

MEdarkolivegreen ~~ MEcyan

MEcyan ~~ MEmagenta

MEcyan ~~ MEblack

MEcyan ~~ MEbrown

MEdarkolivegreen ~~ MEmagenta

MEdarkolivegreen ~~ MEblack

MEdarkolivegreen ~~ MEbrown

MEmagenta ~~ MEblack

MEmagenta ~~ MEbrown

MEblack ~~ MEbrown

'

all_whz_sem_fit <- sem(all_whz_sem, 
                        missing = "ML",
                        likelihood = "wishart",
                        sample.cov.rescale = TRUE,
                        data = hiv_ncc_sam_growth_analysis_normalized_WAZ)


summary(all_whz_sem_fit, fit.measures = TRUE, standardized = TRUE)


node_labeles = list(whz_disch = "Baseline WHZ", 
                    whz_fu90 = "Post-discharge WHZ at Day 90", 
                    MEdarkolivegreen = "PM7",
                    MEcyan = "PM26", 
                    MEmagenta = "PM31",
                    MEblack = "PM33",
                    MEbrown = "PM37",
                    site_recoded = "Site",
                    #malaria = "Malaria", diarrhoea = "Diarrhoea", severe_pneumonia = "Pneumonia",
                    hiv_status_recored = "HIV", 
                    discharge_age = "Discharge age", 
                    sex_recoded = "Sex"
                    )


all_whz_sem_plt <- lavaanPlot(model = all_whz_sem_fit, coefs = TRUE,
                               stand = TRUE, #sig = 0.05, 
                               labels = node_labeles,
                               covs = TRUE,
                               edge_options = list(color = "grey50"),
                               graph_options = list(rankdir = "TB"),
                               stars = "regress")

all_whz_sem_plt


varTable(all_whz_sem_fit)

modindices(all_whz_sem_fit, sort = TRUE)


```

## HAZ SEM MODELS


```{r all_haz_sem, fig.width=4, fig.height=8, warning=FALSE}

set.seed(155)

all_haz_sem <- '
# regression 

MEskyblue3 ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEskyblue ~ hiv_status_recored + discharge_age + sex_recoded  + site_recoded 

MEdarkgreen ~  hiv_status_recored + discharge_age + sex_recoded + site_recoded 

haz_disch ~ MEskyblue3 + MEskyblue + MEdarkgreen + sex_recoded + hiv_status_recored + discharge_age 

haz_fu90 ~ MEskyblue3 + MEskyblue + MEdarkgreen + haz_disch + hiv_status_recored + discharge_age + sex_recoded


# covariance

MEskyblue3 ~~ MEskyblue

MEskyblue3 ~~ MEdarkgreen

MEskyblue ~~ MEdarkgreen

'
all_haz_sem_fit <- sem(all_haz_sem, 
                        missing = "ML",
                        likelihood = "wishart",
                        sample.cov.rescale = TRUE,
                        data = hiv_ncc_sam_growth_analysis_normalized_WAZ)


summary(all_haz_sem_fit, fit.measures = TRUE, standardized = TRUE)


node_labeles = list(haz_disch = "Baseline HAZ", 
                    haz_fu90 = "Post-discharge HAZ at Day 90", 
                    MEdarkgreen = "PM18",
                    MEskyblue = "PM12", 
                    MEskyblue3 = "PM3",
                    site_recoded = "Site",
                    hiv_status_recored = "HIV", 
                    discharge_age = "Discharge age", 
                    sex_recoded = "Sex"
                    )

all_haz_sem_plt <- lavaanPlot(model = all_haz_sem_fit, coefs = TRUE,
                               stand = TRUE, #sig = 0.05, 
                               labels = node_labeles,
                               covs = TRUE,
                               edge_options = list(color = "grey50"),
                               graph_options = list(rankdir = "TB"),
                               stars = "regress")

all_haz_sem_plt


varTable(all_haz_sem_fit)

modindices(all_haz_sem_fit, sort = TRUE)


```










