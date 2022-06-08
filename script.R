# SESSION INFORMATION ==========================================================

# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 20.04.3 LTS

# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/liblapack.so.3

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
# [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
# [10] LC_TELEPHONE=C            LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

# attached base packages:
#  [1] stats graphics grDevices utils datasets methods base     

# loaded via a namespace (and not attached):
#  [1] compiler_4.1.1 tools_4.1.1  





# 1. LOAD LIBS AND DATA, AND COMBINE ===========================================

# Libs
library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(ggpubr)

# Theme
ggthemr::ggthemr("dust")
# Grey:    #9A8A76
# Red:     #DB735C
# Orange:  #EFA86E
# Yellow:  #F3C57B

# Merge datasets containing variables
MyMerge <- function(x, y){
  df <- merge(x, y, by = "pt", all.x = FALSE, all.y = FALSE)
  return(df)}
data <- Reduce(MyMerge, list(data1, data2, data3)); rm(c(data1, data2, data3))

# Convert factors to factors 
factors <- names(data)[c(11, 13:34)]
data[factors] <- lapply(data[factors], factor); rm(factors)





### 2. DATA FILTERING AND EXPLORATION ==========================================

# Remove individuals with CHD event before start of study (baseline)
# Remove individuals with CHD event after end of study (day 4500)
data.coro <- data %>%                                                   
  filter(is.na(coro_days_fromBaseline) | coro_days_fromBaseline > 0 ) %>%                 
  filter(coro_days_fromBaseline < 4500 ) 

# Nb of events for individuals with CHD event in 12 years
ind_in_study <- data.coro$pt
survival.data <- read.table("cv_dates.txt", header = TRUE)
survival.data <- survival.data[survival.data$pt %in% ind_in_study,] %>%
  select(pt, datecoro)  %>%
  na.omit()

# Count nb of events per individual
occurrence_tot <- as.data.frame(t(table(table(survival.data['pt']))))[2:3]
names(occurrence_tot) <- c("Occurrence", "Count")
occurrence_tot; rm(survival.data, survival.data, ind_in_study, occurrence_id)

# Plot histogram
ggplot(data = occurrence_tot, aes(x = Occurrence, y = Count)) +
  geom_bar(stat = "identity") +
  xlab("Number of CHD events") +
  ylab("Individuals") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1))

# Table 1 - Baseline characteristics
tableone::CreateTableOne(vars = names(data.coro), strata = "event_coro" , 
                         factorVars = c("event_coro", "Statin", 
                                        "Monthlyhouseholdgrossincome"),
                         includeNA = TRUE,
                         addOverall = TRUE,
                         data = data.coro)





### 3. UNIVARIABLE COX REGRESSIONS =============================================
# Variable selection using univariable Cox proportional HR analyses

# Univariable cox regression for all included variables, print P-values
results_univ <- lapply(data.coro[,c(4:35)], function(x){
  summary(coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ x, 
                data = data.coro))$coefficients[1, c(1,5)]})

# Print variables with P < 0.05
sign_univ <- names(results_univ[sapply(results_univ, function(x) x[2] < 0.05)])
select(data.coro, c("coro_days_fromBaseline", "event_coro", sign_univ))





### 4. MULTIVARIABLE COX REGRESSION ANALYSIS ===================================
# Multivariable Cox proportional HR analysis with stepwise selection

# Multivariable cox regression 
cox.model <- coxph(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ .,
                        data = select(data.coro, c("coro_days_fromBaseline", 
                                                   "event_coro", sign_univ))) 

# Print summary
gtsummary::tbl_regression(cox.model)

# Step-wise reduction of the model
cox.model.step <- stats::step(cox.model, direction = "back")
gtsummary::tbl_regression(cox.model.step)
summary(cox.model.step)

# Forest plot (figure 2)
forestmodel::forest_model(cox.model.step,
                          factor_separate_line = TRUE) +
  theme(axis.text.x = element_text(size = 17))





### 5. CHECK PROPORTIONAL HAZARD ASSUMPTIONS ===================================

# Test the PH assumption, and plot 
cox.zph(cox.model.step, transform = "km", global = TRUE)
ggcoxzph(cox.zph(cox.model.step), ggtheme = theme_survminer(),
         point.col = "#DB735C")
residuals(zph)

# Extract residuals and time
zph <- cox.zph(cox.model.step)
plot(zph$time, zph$y[1:nrow(zph$y)], ylim = c(-10,10))

# Assess co-linearity
rms::vif(cox.model.step)





### 6. COX REGRESSION ANALYSIS FOR PATHOGEN BURDEN =============================

# Barplot: pathogen burden count
ggplot(data.coro, aes(Total_burden)) +
  geom_bar() +
  geom_text(stat='count', aes(label = ..count..), vjust = -1, size = 6) + 
  xlab("Pathogen burden") +
  ylab("Count") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
        panel.spacing = unit(0.05, "lines"),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.justification = c(1, 1), 
        legend.position = c(1, 1))

# Cox regression 
cox.model.burden <- coxph(Surv(coro_days_fromBaseline, 
                               as.numeric(event_coro)) ~ .,
                          data = select(data.coro, 
                                        c("coro_days_fromBaseline",
                                          "event_coro", 
                                          sign_univ[1:6], 
                                          "Total_burden"))) 
gtsummary::tbl_regression(cox.model.burden)

# Kaplan-Meier curve for pathogen burden 
fit.burden <- survfit(Surv(coro_days_fromBaseline, as.numeric(event_coro)) ~ 
                        Total_burden,
                      data = data.coro) 
survminer::ggsurvplot(fit.burden, data = data.coro)





### 7. RESET ===================================================================

# Reset theme
ggthemr::ggthemr_reset()

