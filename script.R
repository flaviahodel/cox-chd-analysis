# SESSION INFORMATION ==========================================================

# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6

# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

# locale:
# [1] fr_CH.UTF-8/fr_CH.UTF-8/fr_CH.UTF-8/C/fr_CH.UTF-8/fr_CH.UTF-8
 

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] ggcorrplot_0.1.4  psych_2.2.9       corrplot_0.92     DescTools_0.99.47 survminer_0.4.9   ggpubr_0.5.0     
# [7] survival_3.3-1    ggmosaic_0.3.3    forcats_0.5.2     stringr_1.4.1     purrr_0.3.5       readr_2.1.3      
# [13] tidyr_1.2.1       tibble_3.1.8      ggplot2_3.4.0     tidyverse_1.3.2   dplyr_1.0.10     




# 1. LOAD LIBS AND DATA, AND COMBINE ===========================================

# Libs
library(dplyr)
library(tidyverse)
library(ggmosaic)
library(tidyr)
library(ggplot2)
library(survival)
library(ggpubr)
library(survminer)
library(DescTools)
library(corrplot)
library(psych)
library(ggcorrplot)


# Theme
ggthemr::ggthemr("dust")
#   Grey:    #9A8A76
#   Red:     #DB735C
#   Orange:  #EFA86E
#   Yellow:  #F3C57B


# Survival data 
# Update 'stop_days_fromBaseline' if there is a more recent event in CV dataset
setwd("data")
survival.data <- read.table("./data_cox_complete.txt", header = TRUE) %>%
  dplyr::select(-sex)
survival.data.updated <- survival.data %>% 
  transform(stop_days_fromBaseline = pmax(coro_days_fromBaseline, aomi_days_fromBaseline, avc_days_fromBaseline, stop_days_fromBaseline, na.rm = TRUE)) 

# Add an event description
# if there is no event, the last follow-up is the censor date
survival.data.updated$event_coro <- ifelse(is.na(survival.data.updated$coro_days_fromBaseline), 0, 1)
survival.data.updated$coro_days_fromBaseline[which(is.na(survival.data.updated$coro_days_fromBaseline))] <- survival.data.updated$stop_days_fromBaseline[which(is.na(survival.data.updated$coro_days_fromBaseline))]

# Serologies, pcs and inflammation biomarkers 
sero.pc.inflammation.data <- read.csv("./colaus.sero.all.csv") %>%
  dplyr::select(-sex, -age, -BMI, -BMI, -phyact, -smoking, -PRS, -hs.CRP, -IL.1b, -IL.6, -TNF.a) %>%
  rename(hs.CRP.log10 = "log10.hs.CRP.", 
         IL.1b.log10 = "log10.IL.1b.",
         IL.6.log10 = "log10.IL.6.", 
         TNF.a.log10 = "log10.TNF.a.")

# Polygenic risk scores 
prs_coro.data <- read.table("./CoLaus.PRS.v2.best", header = TRUE) %>%
  dplyr::select(-FID, -In_Regression) %>%
  dplyr::rename(pt = "IID", PRS_coro = "PRS")

# Covariates
cov.data <- read.table("./CoLaus_var.txt", header = TRUE) %>%
  dplyr::select(-hsCRP, -IL1b, -IL6, -TNFa, -Death, -HDL_cholesterol, -Cholesterol, -LDL_cholesterol) %>%
  dplyr::rename(pt = "ID")

# Hdl cholesterol
hdl <- read.table("./hdl_cholesterol.txt", header = TRUE) %>%
  dplyr::rename(HDL_cholesterol = "hdlch", LDL_cholesterol = "ldlch", Total_cholesterol = "chol")

# Score2
score2 <- read.table("./SCORE2.values.txt", header = TRUE)
score2$score2_z <- score2$score2
score2$score2_z <- as.numeric(scale(score2$score2_z))

# Pathogen burden
tot.burden <- read.table("./colaus.sero.all_burden.txt", header = TRUE)

# Social covariates
social <- read.table("./CoLaus.social_covariates.txt", header = TRUE) %>%
  dplyr::rename(pt = "F1pt")


# Merge
MyMerge <- function(x, y){
  df <- merge(x, y, by = "pt", all.x = FALSE, all.y = FALSE)
  return(df)}
data <- Reduce(MyMerge, list(survival.data.updated, sero.pc.inflammation.data, tot.burden, prs_coro.data, cov.data, hdl, score2, social))
#rm(survival.data, survival.data.updated, sero.pc.inflammation.data, tot.burden, prs_coro.data, cov.data, hdl, score2, social)

# Columns of interest
data <- data %>%
  dplyr::select(pt, event_coro, coro_days_fromBaseline, Age, score2_z, PRS_coro, 
                hs.CRP.log10, IL.1b.log10, IL.6.log10, TNF.a.log10, 
                Statin, Monthlyhouseholdgrossincome, AVS, AI, Professional_activity,
                BKPyV, JCPyV, WUPyV, HPyV6, HSV_1, HSV_2, VZV, HHV_7, 
                C_trachomatis, S_gallolyticus, VP1_unique, EBV, KSHV,
                HCMV,  HHV_6A, HHV_6B, T_gondii, H_pylori, F_nucleatum,
                Total_burden)

# Convert factors to factors 
factors <- c("Statin", "AVS", "AI", "Professional_activity", 
             "BKPyV", "JCPyV", "WUPyV", "HPyV6", "HSV_1", "HSV_2", 
             "VZV", "HHV_7", "C_trachomatis", "S_gallolyticus", 
             "VP1_unique", "EBV", "KSHV", "HCMV",  "HHV_6A", "HHV_6B", 
             "T_gondii", "H_pylori", "F_nucleatum")
data[factors] <- lapply(data[factors], factor); rm(factors)




### 2. FILTERING DATA ----------------------------------------------------------

# Remove individuals with event BEFORE start of study (baseline)
data.coro <- data %>%                                                   # before, N=3681
  filter(is.na(coro_days_fromBaseline) | coro_days_fromBaseline > 0 )   # after,  N=3585

# Remove individuals with event AFTER end of study (day 4500)
data.coro <- data.coro %>%                 # before, N = 3585
  filter(coro_days_fromBaseline < 4500 )   # after,  N = 3459

# Nb of events for individuals with CHD event in 12 years
ind_in_study <- data.coro$pt
survival.data <- read.table("./colaus_cv_dates.txt", header = TRUE)
survival.data <- survival.data[survival.data$pt %in% ind_in_study,] %>%
  dplyr::select(pt, datecoro)  %>%
  na.omit()

# How many events per individual?
occurrence_tot <- as.data.frame(t(table(table(survival.data['pt']))))[2:3]
names(occurrence_tot) <- c("Occurrence", "Count")
rm(survival.data, survival.data, ind_in_study, occurrence_id)
occurrence_tot

# Histogram
ggplot(data=occurrence_tot, aes(x=Occurrence, y=Count)) +
  geom_bar(stat="identity") +
  xlab("Number of CHD events") +
  ylab("Individuals") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1))


# Nb of individuals by monthly income

# Total
round(table(data.coro$Monthlyhouseholdgrossincome)/nrow(data.coro)*100,2)
round((nrow(data.coro)-(178+452+552+504+338+344))/nrow(data.coro)*100,2)

# Controls
table(is.na(data.coro$Monthlyhouseholdgrossincome[data.coro$event_coro=="0"]))
table(data.coro[data.coro$event_coro=="0", ]$Monthlyhouseholdgrossincome, exclude = NULL)
round(table(data.coro[data.coro$event_coro=="0", ]$Monthlyhouseholdgrossincome, exclude = NULL)/nrow(data.coro[data.coro$event_coro=="0", ])*100, 2)

# Cases 
table(is.na(data.coro$Monthlyhouseholdgrossincome[data.coro$event_coro=="1"]))
table(data.coro[data.coro$event_coro=="1", ]$Monthlyhouseholdgrossincome, exclude = NULL)
round(table(data.coro[data.coro$event_coro=="1", ]$Monthlyhouseholdgrossincome, exclude = NULL)/nrow(data.coro[data.coro$event_coro=="1", ])*100, 2)

# Table 1
tableone::CreateTableOne(vars = names(data.coro), strata = "event_coro" , 
                         factorVars = c("event_coro", "Statin", "Monthlyhouseholdgrossincome"),
                         includeNA = TRUE,
                         addOverall = TRUE,
                         data = data.coro)



### 3. UNIVARIABLE COX REGRESSIONS ----------------------------------

# For all non-missing columns
lapply(data.coro[,c(4:7, 11, 16:35)], function(x){
  summary(coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ x, 
                data = data.coro))$coefficients[1, c(1,5)]})


# To get HR etc
forestmodel::forest_model(coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ Total_burden, 
                                data = data.coro), factor_separate_line = TRUE)

summary(coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ Total_burden,
              data = data.coro))$coefficient

# Variables with NAs:
# a) Biomarkers of inflammation
# IL.1b.log10, IL.6.log10, and TNF.a.log10

data.IL.1b <- data.coro[!is.na(data.coro$IL.1b.log10), ]
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ IL.1b.log10, data = data.IL.1b) %>%
  gtsummary::tbl_regression(); rm(data.IL.1b)

data.il6 <- data.coro[!is.na(data.coro$IL.6.log10), ]
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ IL.6.log10, data = data.il6) %>%
  gtsummary::tbl_regression(); rm(data.il6)

data.tnf.a <- data.coro[!is.na(data.coro$TNF.a.log10), ]
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ TNF.a.log10, data = data.tnf.a) %>%
  gtsummary::tbl_regression(); rm(data.tnf.a)


# b) Socio-economic factors
# Monthlyhouseholdgrossincome, AVS, AI, Professional_activity
data.MHGI <- data.coro[!is.na(data.coro$Monthlyhouseholdgrossincome), ] # 1091 NAs
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ Monthlyhouseholdgrossincome, data = data.MHGI[data.MHGI$Statin == "0", ]) %>%
  gtsummary::tbl_regression(); rm(data.MHGI)

data.AVS <- data.coro[!is.na(data.coro$AVS), ] # 612 NAs
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ AVS, data = data.AVS[data.AVS$Statin == "0", ]) %>%
  gtsummary::tbl_regression(); rm(data.AVS)

data.AI <- data.coro[!is.na(data.coro$AI), ] # 672 NAs
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ AI, data = data.AI) %>%
  gtsummary::tbl_regression(); rm(data.AI)

data.PA <- data.coro[!is.na(data.coro$Professional_activity), ] # 587 NAs
coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ Professional_activity, data = data.PA) %>%
  gtsummary::tbl_regression(); rm(data.PA)




### 4. MULTIVARIABLE COX REGRESSION --------------------------------------------
# Exclude individuals with NAs in TNF-a and MHGI
data.coro.na <- data.coro %>% 
  tidyr::drop_na("Monthlyhouseholdgrossincome", "TNF.a.log10") 

# Rename variables for esthetics  
cox.data <- data.coro.na
colnames(x = cox.data)[colnames(x = cox.data)== "score2_z"] <- 'SCORE2'
colnames(x = cox.data)[colnames(x = cox.data)== "Monthlyhouseholdgrossincome"] <- 'Income'
colnames(x = cox.data)[colnames(x = cox.data)== "PRS_coro"] <- 'CHD-PRS'
colnames(x = cox.data)[colnames(x = cox.data) == "hs.CRP.log10"] <- 'hs-CRP'
colnames(x = cox.data)[colnames(x = cox.data) == "TNF.a.log10"] <- 'TNF-a'
colnames(x = cox.data)[colnames(x = cox.data) == "HSV_1"] <- 'HSV-1'
colnames(x = cox.data)[colnames(x = cox.data) == "C_trachomatis"] <- 'C. trachomatis'
colnames(x = cox.data)[colnames(x = cox.data) == "F_nucleatum"] <- 'F. nucleatum'
colnames(x = cox.data)[colnames(x = cox.data) == "HHV_6A"] <- 'HHV-6A'

# Multivariable cox regression 
cox.model.2324 <- coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ SCORE2 +
                          Statin + Income +
                          `CHD-PRS` + `hs-CRP` + `TNF-a` + 
                          `F. nucleatum`  + `C. trachomatis` + HPyV6 + `HSV-1` +  `HHV-6A` +   VZV,
                        data = cox.data) 
gtsummary::tbl_regression(cox.model.2324)

# A Stepwise Reduction of our model
cox.model.2324.step <- stats::step(cox.model.2324, direction = "back")
gtsummary::tbl_regression(cox.model.2324.step)
summary(cox.model.2324.step)
# --> Significant: SCORE2, PRS, F. nucleatum 


# Forest model
fig2 <- forestmodel::forest_model(cox.model.2324.step, 
                                  factor_separate_line = TRUE) +
  theme(axis.text.x = element_text(size=17))

build_fig2 <- ggplot_build(fig2)



### 5. CHECK PH ASSUMPTIONS ----------------------------------------------------

# Testing the Proportional Hazards Assumption
cox.zph(cox.model.2324.step, transform = "km", global = TRUE)
fig_s3 <- ggcoxzph(cox.zph(cox.model.2324.step), ggtheme = theme_survminer(), 
                   point.col = "#DB735C")

# Get residuals and time
zph <- cox.zph(cox.model.2324.step)
zph$time
zph$y
plot(zph$time, zph$y[1:nrow(zph$y)], ylim = c(-10,10))


# Assessing co-linearity
# Perhaps we have some co-linearity here, which might imply that we could sensibly fit a smaller model, 
# We should probably be sticking to a model with no more than 2 or perhaps as many as 3 coefficients to be estimated.
rms::vif(cox.model.2324.step)





### 6. COX FOR PATHOGEN BURDEN -------------------------------------------------

# Plot
fig_S4 <- ggplot(data.coro.na, aes(Total_burden)) +
  geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1, size = 6) + 
  xlab("Pathogen burden") +
  ylab("Count") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1))

ggplot_build(fig_S4)

cox.model.burden <- coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ score2_z +
                            PRS_coro + hs.CRP.log10 + TNF.a.log10 + Statin +  
                            Monthlyhouseholdgrossincome + Total_burden,
                          data = data.coro.na) 
gtsummary::tbl_regression(cox.model.burden)

mycolors = c(RColorBrewer::brewer.pal(name="Dark2", n = 9), RColorBrewer::brewer.pal(name="Paired", n = 9))
fit.burden <- survfit(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ 
                        Total_burden,
                      data = data.coro.na) 
#survminer::ggsurvplot(fit.burden, data = data.coro.na, color = mycolors)

## Boxplots
boxplot(data.coro$Age ~ data.coro$AVS)
boxplot(data.coro$Age ~ data.coro$AVS)





### 7. ADDITIONAL FIGURES ----------------------------------------------------------

### Suppl. figure 1 = Correlation: quantitative variables 
dat_pred_quant <- data.coro.na %>% dplyr::select( "score2_z", "PRS_coro", "hs.CRP.log10", "TNF.a.log10")

names(dat_pred_quant) <- recode(names(dat_pred_quant), 
                                PRS_coro = "CHD-PRS",
                                hs.CRP.log10 = "hs-CRP",
                                TNF.a.log10 = "TNF-a",
                                score2_z = "SCORE2")

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)}
pairs(dat_pred_quant, diag.panel=panel.hist, pch=16, col="steelblue")


fig_S1 <- pairs.panels(dat_pred_quant, 
                       method = "pearson", 
                       hist.col = "#DB735C",
                       density = TRUE,  
                       ellipses = FALSE, 
                       lm = TRUE, 
                       rug = FALSE, 
                       breaks = 12, 
                       cex.cor = 0.4, 
                       ci = TRUE, 
                       col = "#DB735C", 
                       pch = 20, 
                       bg = "lightgray", 
                       cex.axis = 1.5, 
                       digits = 2)

ggplot_build(fig_S1)


### Suppl. figure 2 = Correlation: categorical variables 
dat_pred_cat <- data.coro %>% dplyr::select("Statin", "Monthlyhouseholdgrossincome", "HSV_1", "HHV_6A", "HPyV6", "F_nucleatum", "C_trachomatis", "VZV")  %>% mutate_all(as.character)

corr_categorical <- DescTools::PairApply(dat_pred_cat, DescTools::CramerV)
corrplot::corrplot(DescTools::PairApply(dat_pred_cat, DescTools::CramerV), diag = FALSE, type = "upper")

fig_S2 <- corr_categorical %>% 
  ggcorrplot(show.diag = F, type="lower", lab=TRUE, lab_size=7,
             #colors = c("white", "#DB735C"),
             method = "square",
             theme_replace()) + 
  scale_fill_gradient2(limit = c(0,1), low = "white", high =  "#DB735C") + 
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, 'cm')) +
  guides(fill=guide_legend(title="Cramer's V")) +
  scale_x_discrete(labels=c("Monthlyhouseholdgrossincome" = "Income", "HSV_1" = "HSV-1", "HHV_6A" = "HHV-6A", "HPyV6" = "HPyV6", "F_nucleatum" = "F. nucleatum", "C_trachomatis" = "C. trachomatis", "VZV" = "VZV"), 
                   guide=guide_axis(angle=20)) +
  scale_y_discrete(labels=c("Statin" = "Statin", "Monthlyhouseholdgrossincome" = "Income", "HSV_1" = "HSV-1", "HHV_6A" = "HHV-6A", "HPyV6" = "HPyV6", "F_nucleatum" = "F. nucleatum", "C_trachomatis" = "C. trachomatis"))
ggplot_build(fig_S2)


### CRP plot
# Histogram
ggplot(data.coro, aes(hs.CRP.log10)) +
  geom_bar(bins = 10, ) +
  scale_x_binned() +
  xlab("hsCRP levels") +
  ylab("Individuals") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1))


### 8. ADDITIONAL WORK ----------------------------------------------------------

## hs-CRP
data <- data.coro[, c(34, 7)]
data <- na.omit(data)
ggboxplot(data, x = "F_nucleatum", y = "hs.CRP.log10",
          color = "F_nucleatum", palette = "jco",
          add = "jitter") + 
  stat_compare_means(method = "t.test")

## TNF-a
data <- data.coro[, c(34, 10)]
data <- na.omit(data)
ggboxplot(data, x = "F_nucleatum", y = "TNF.a.log10",
          color = "F_nucleatum", palette = "jco",
          add = "jitter") + 
  stat_compare_means(method = "t.test")


## Is inflammation is predictive of CHD in multivariable models excluding F. nucleatum?
summary(lm(coro_days_fromBaseline ~ hs.CRP.log10, data = data.coro))
summary(lm(coro_days_fromBaseline ~ hs.CRP.log10 * F_nucleatum, data = data.coro))


# Multivariable cox regression without Fn
cox.model.2324 <- coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ `hs-CRP` + `F. nucleatum`,
                        data = cox.data) 
gtsummary::tbl_regression(cox.model.2324)


## Testing interaction 
# Multivariable cox regression 
cox.model.2324 <- coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ SCORE2 *
                          Statin * Income *
                          `CHD-PRS` * `hs-CRP` * `TNF-a` * 
                          `F. nucleatum`  * `C. trachomatis` * HPyV6 * `HSV-1` *  `HHV-6A` *   VZV,
                        data = cox.data) 
gtsummary::tbl_regression(cox.model.2324)






### 8. RESET ===================================================================

# Reset theme
ggthemr::ggthemr_reset()

