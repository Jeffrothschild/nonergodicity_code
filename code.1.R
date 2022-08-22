library(tidyverse)
library(modeltime)
library(tidymodels)
library(ggridges)
library(ggpmisc)
library(Cairo)
library(patchwork)

my.formula <- y ~ poly(x, 1, raw = TRUE)

data_subset <- read_rds("synthetic_data.rds")


arima_ci_fn <- function(dv, iv1, iv2, date = date, data, ...){
  
  df <- data %>% dplyr::select({{dv}}, {{date}}, {{iv1}}, {{iv2}},  ... )
  
  frmla1 <- formula(paste(colnames(df)[1], " ~ ", colnames(df)[2], "+",  colnames(df)[3],  "+",
                          colnames(df)[4] , collapse = " "))
  
  if(ncol(df) == 5){
    frmla1 <- formula(paste(colnames(df)[1], " ~ ", colnames(df)[2], "+",  colnames(df)[3],  "+",
                            colnames(df)[4] , "+",  colnames(df)[5] , collapse = " "))
  }
  
  if(ncol(df) == 6){
    frmla1 <- formula(paste(colnames(df)[1], " ~ ", colnames(df)[2], "+",  colnames(df)[3],  "+",
                            colnames(df)[4] , "+",  colnames(df)[5] , "+",  colnames(df)[6] , collapse = " "))
  }
  
  if(ncol(df) == 7){
    frmla1 <- formula(paste(colnames(df)[1], " ~ ", colnames(df)[2], "+",  colnames(df)[3],  "+",
                            colnames(df)[4] , "+",  colnames(df)[5] , "+",  colnames(df)[6] , "+",  colnames(df)[7] , collapse = " "))
  }
  
  other_vars <- df %>%  select(4:last_col()) %>% colnames() %>% paste0(., collapse = ", ")
  
  model_fit_auto_arima <- arima_reg() %>%
    set_engine("auto_arima") %>%
    fit(frmla1,
        data = data)
  
  test_pred <- 
    predict(model_fit_auto_arima, df) %>% 
    bind_cols(df) %>% 
    select(.pred, AM_zq_score) 
  
  reg_metrics <- metric_set(rsq)
  
  individ_rsq <- test_pred %>% reg_metrics(AM_zq_score, .pred) %>% 
    pull(.estimate) %>% round(.,2)
  
  model_fit_auto_arima$fit$models$model_1 %>% tidy(conf.int = T) %>% 
    filter(term == colnames(df)[3]) %>% 
    select(-1) %>% 
    rename("arima_estimate" = 1,
           "arima_se" = std.error ,
           "arima_ci_low" = conf.low,
           "arima_ci_high" = conf.high) %>% 
    mutate(
      arima_dv = colnames(df)[1],
      arima_iv = colnames(df)[3],
      arima_other_vars = other_vars, 
      arima_r2 = individ_rsq
      
    )
}
spear_w_ci_fn <- function(x,y) {DescTools::SpearmanRho(x,y, use = c("pairwise.complete.obs"), 
                                                       conf.level = .95) %>% as.data.frame() %>%
    rownames_to_column() %>% as_tibble() %>% pivot_wider(names_from = 1, values_from = 2) %>% 
    rename("estimate" = rho,
           "conf.low" =  lwr.ci,
           "conf.high" = upr.ci)
}
nested_corr_function <- function(df, x, y){
  x_quo = enquo(x)
  y_quo = enquo(y)
  

  
  df %>% 
    group_by(subject_id) %>%
    nest() %>%
    mutate(
      shapiro1 = map_dbl(data, ~ shapiro.test(.x[[!!x_quo]])$p.value),
      shapiro2 = map_dbl(data, ~ shapiro.test(.x[[!!y_quo]])$p.value),
      type = ifelse(shapiro1 < .05 | shapiro2 < .05, "spearman", "pearson"),
      cor_pear = map(data, ~ cor.test(.x[[!!x_quo]], .x[[!!y_quo]], use = "pairwise.complete.obs")),
      cor_pear = map(cor_pear, broom::tidy),
      cor_spear =  map(data,  ~spear_w_ci_fn(.x[[!!x_quo]], .x[[!!y_quo]]))
      
    ) %>%
    unnest(cor_pear) %>% select(-c(statistic, parameter, alternative, method, p.value)) %>%
    rename(cor_val_pear = "estimate",
           ci_low_pear = "conf.low",
           ci_high_pear = "conf.high") %>% 
    unnest(cor_spear) %>%
    rename(cor_val_spear = "estimate",
           ci_low_spear = "conf.low",
           ci_high_spear = "conf.high") %>% 
    mutate(
      cor_val = ifelse(type == "spearman", cor_val_spear, cor_val_pear),
      conf.low = ifelse(type == "spearman", ci_low_spear, ci_low_pear),
      conf.high = ifelse(type == "spearman", ci_high_spear, ci_high_pear),
    ) %>% 
    mutate(
      sign = case_when(
        conf.low * conf.high < 0 ~ "NS",
        conf.low * conf.high >0 & conf.high >0 ~ "+",
        conf.low * conf.high >0 & conf.high < 0  ~ "-"
      )) %>%
    
    ggplot(aes(cor_val, fct_rev(fct_reorder(factor(subject_id), cor_val)), color = sign, shape = type))+
    geom_vline(xintercept=0, linetype = "dotted") +
    geom_point(aes(), size = 2.5) +
    geom_errorbarh(aes(xmin=conf.low, xmax=conf.high,height = 0.4)) +
    scale_color_manual(values = c("-" = "#801d1d",
                                  "NS" = "grey82",
                                  "+"= "#1d801d")) +
    theme(
      legend.position = "none"
    )
  
}


# correlations ------------------------------------------------------------

# *carb and PRS scatterplot ---------------------------------------------------------

rm_cor <- rmcorr::rmcorr(factor(subject_id), AM_zq_score, lag1_diet_carb_g_kg, data_subset)
rm_r <- round(rm_cor$r, 2)
rm_ci_low <- round(rm_cor$CI[1],2)
rm_ci_hi <- round(rm_cor$CI[2],2)
rm_p <- if(rm_cor$p < 0.001){"< 0.001"} else {rm_cor$p} 


combined_scatter <- data_subset %>%
  ggplot(aes(lag1_diet_carb_g_kg, AM_zq_score))+
  geom_point(size = 1.5, colour = "grey70", alpha = 0.5) +  #6f6fde
  annotate("text", x = 17, y = 98, label = paste0("r = ", rm_r, ", 95% CI ", rm_ci_low, " to ", rm_ci_hi),
           size = 4, hjust = 1) +
  stat_smooth(method = "lm", aes(),  formula =my.formula, se = TRUE, linetype="solid", size=0.7, alpha = 0.1, color = "black") +
  labs(x = "Prior day CHO (g/kg)", y = "AM PRS score")

combined_scatter



# * individ correlations ----------------------------------------------------------
all_corr_plot <- nested_corr_function(data_subset, "lag1_diet_carb_g_kg", "AM_zq_score") +
  labs(x = "Correlation between prior day CHO\nand AM PRS score", y = "Participant ID") 

all_corr_plot



# *facet highlights -------------------------------------------------------
facet_6_colors <- c("#801d1d", "#801d1d", "#801d1d", "#1d801d", "#1d801d", "#1d801d")

facet_6_plot <- data_subset %>% 
  filter(subject_id %in% c(29,23,12,39,2,18)) %>% 
  mutate(subject_id = factor(subject_id, levels = c("29","23","12","39","2","18"))) %>% 
  ggplot(aes(lag1_diet_carb_g_kg, AM_zq_score, color = subject_id))+
  geom_point(size = 1.5, alpha = .75) +
  gghighlight::gghighlight(use_direct_label = F, unhighlighted_params = list(color = "gray90")) +
  scale_color_manual(values = facet_6_colors)+
  stat_smooth(method = "lm", aes(),  formula =my.formula, se = F, linetype="solid", size=0.7, alpha = 0.1, color = "black") +
  
  labs(x = "Prior day CHO (g/kg)", y = "AM PRS score") +
  facet_wrap( ~ subject_id, ncol = 3) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(margin = margin(.9,0,.9,0, "mm"))
  )

facet_6_plot


# * patch -------------------------------------------------------------------
((combined_scatter /facet_6_plot)| all_corr_plot) + 
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(title = NULL, tag_levels = list(c('A', 'C', 'B'))) &
  theme(plot.tag = element_text(size = 17, face="bold"))



# group level analysis ----------------------------------------------------
group_daily_tbl <- data_subset %>% 
  group_by(subject_id) %>% 
  slice(8:47) %>% 
  ungroup() %>% 
  group_by(study_day) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(
    lm1 = map(data, ~lm(AM_zq_score ~ lag1_diet_carb_g_kg + lag1_roll_strain + lag1_exercise_load + AM_soreness + sleep_index, .x)),
    lm1_metrics = map(lm1, ~glance(.x)),
    lm1 = map(lm1, ~broom::tidy(.,conf.int = T)),
  ) %>% 
  unnest(lm1_metrics) %>% select(-c(adj.r.squared, sigma:nobs)) %>% 
  unnest(lm1) %>% select(-c( statistic)) %>%
  filter(term == "lag1_diet_carb_g_kg")

group_daily_tbl 

group_summary_tbl <- group_daily_tbl %>% 
  summarise(
    mean_coef = mean(estimate),
    mean_coef_SD = sd(estimate),
    mean_coef_ci_low = confintr::ci_mean(estimate)$interval[1],
    mean_coef_ci_high = confintr::ci_mean(estimate)$interval[2],
    mean_sd_ci_low = confintr::ci_sd(estimate)$interval[1],
    mean_sd_ci_high = confintr::ci_sd(estimate)$interval[2],
    mean_rsq = mean(r.squared),
    mean_rsq_SD = sd(r.squared),
    mean_rsq_ci_low = confintr::ci_mean(r.squared)$interval[1],
    mean_rsq_ci_high = confintr::ci_mean(r.squared)$interval[2],
    mean_rsq_sd_ci_low = confintr::ci_sd(r.squared)$interval[1],
    mean_rsq_sd_ci_high = confintr::ci_sd(r.squared)$interval[2],
  )  %>% 
  mutate(
    Level = "Group",
    .before = everything()
  )

group_summary_tbl


# individual level analysis --------------------------------------------------------

individ_arima_tbl <- data_subset %>%  
  group_by(subject_id) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(
    arima_model = map(data, ~arima_ci_fn(AM_zq_score, lag1_diet_carb_g_kg,  lag1_exercise_load,
                                         date, data = .x, AM_soreness, sleep_index, lag1_roll_strain)),  
    lm1 = map(data, ~lm(AM_zq_score ~ lag1_diet_carb_g_kg + lag1_exercise_load + AM_soreness + lag1_roll_strain + sleep_index, .x)),
    lm1 = map(lm1, ~broom::tidy(.,conf.int = T)),
  ) %>% 
  unnest(arima_model) 

individ_arima_tbl

individ_arima_tbl %>% 
  mutate(
    arima_sign = case_when(
      arima_ci_low * arima_ci_high < 0 ~ "NS",
      arima_ci_low * arima_ci_high >0 & arima_ci_high >0 ~ "+",
      arima_ci_low * arima_ci_high >0 & arima_ci_high < 0  ~ "-"
    )) %>%
  
  ggplot(aes(arima_estimate, fct_rev(fct_reorder(factor(subject_id), arima_estimate)), color = arima_sign))+
  geom_vline(xintercept=0, linetype = "dotted") +
  geom_point(aes(), size = 2.5) +
  geom_errorbarh(aes(xmin=arima_ci_low, xmax=arima_ci_high,height = 0.4)) +
  scale_color_manual(values = c("-" = "#801d1d",
                                "NS" = "grey82",
                                "+"= "#1d801d")) +
  labs(x = "CHO model coefficients", y = "Participant ID", title = NULL)+
  theme(
    legend.position = "none",
  )



individ_arima_tbl

individ_summary_tbl <-individ_arima_tbl %>% 
  summarise(
    mean_coef = mean(arima_estimate),
    mean_coef_SD = sd(arima_estimate),
    mean_coef_ci_low = confintr::ci_mean(arima_estimate)$interval[1],
    mean_coef_ci_high = confintr::ci_mean(arima_estimate)$interval[2],
    mean_sd_ci_low = confintr::ci_sd(arima_estimate)$interval[1],
    mean_sd_ci_high = confintr::ci_sd(arima_estimate)$interval[2],
    mean_rsq = mean(arima_r2),
    mean_rsq_SD = sd(arima_r2),
    mean_rsq_ci_low = confintr::ci_mean(arima_r2)$interval[1],
    mean_rsq_ci_high = confintr::ci_mean(arima_r2)$interval[2],
    mean_rsq_sd_ci_low = confintr::ci_sd(arima_r2)$interval[1],
    mean_rsq_sd_ci_high = confintr::ci_sd(arima_r2)$interval[2],
  ) %>% 
  mutate(
    Level = "Individual"
  )

individ_summary_tbl


# combine -----------------------------------------------------------------

combined_summary_tbl <- group_summary_tbl %>% 
  bind_rows(individ_summary_tbl) %>% 
  mutate(
    across(c(mean_coef_SD, mean_sd_ci_low, mean_sd_ci_high), ~round(., 1)),
    across(c(mean_coef_SD, mean_sd_ci_low, mean_sd_ci_high), ~format(., nsmall = 1)),
    across(c(mean_coef, mean_coef_ci_low, mean_coef_ci_high, mean_rsq:mean_rsq_sd_ci_high), ~round(., 2)),
    across(c(mean_coef, mean_coef_ci_low, mean_coef_ci_high, mean_rsq:mean_rsq_sd_ci_high), ~format(., nsmall = 2)),
    "CHO\nMean (95% CI)" = paste0(mean_coef, "\n(", mean_coef_ci_low, " to ", mean_coef_ci_high, ")"),
    "CHO\nSD (95% CI)" = paste0(mean_coef_SD, "\n(", mean_sd_ci_low, " to ", mean_sd_ci_high, ")"),
  ) %>% 
  select(-c(mean_coef:mean_rsq_sd_ci_high))

# combined_summary_tbl

df <- tibble(x = 1, y = 1, tb = list(combined_summary_tbl))


# plot --------------------------------------------------------------------
combined_tbl <-  individ_arima_tbl %>% select(arima_estimate) %>% 
  mutate(Level = "Individual") %>% 
  rename("estimate" = arima_estimate) %>% 
  bind_rows(
    group_daily_tbl %>% select(estimate) %>% mutate(Level = "Group")
  ) 

overlapped <- combined_tbl %>% 
  ggplot(aes(estimate, fill = Level))+
  geom_density(alpha = .5)+
  scale_fill_manual(values = c("#6f6fde", "#de6f6f")) +
  labs(x = "CHO model coefficient", y = "Density", fill = NULL)+
  guides(fill=guide_legend(nrow=2)) +
  theme(
    legend.position = c(0.1, .9),
  )


overlapped +
  geom_table_npc(data = df, aes(label = list(combined_summary_tbl)),
                 table.theme = ttheme_gtlight(base_size = 10, padding = unit(c(4, 3), "mm"), parse = T),
                 npcx = 0.9, npcy = 0.98
  )






# decision tree -----------------------------------------------------------
library(rpart.plot) 
library(finetune)

tbl_for_tree <- read_rds("synthetic_tree_data.rds")
  

set.seed(123)
tree_boots <- bootstraps(tbl_for_tree, times = 500)
tune_boots <- bootstraps(tbl_for_tree, times = 100)

class_bal_table <- table(tbl_for_tree$arima_sign)
class_bal_table

tree_spec <-
  decision_tree(
    cost_complexity = tune(),
    tree_depth = tune(),
    min_n = tune()
  ) %>%
  set_mode("classification") %>%
  set_engine("rpart")



tree_rec <- recipe(arima_sign ~ ., tbl_for_tree) %>% 
  step_other(all_nominal_predictors()) %>%
  themis::step_upsample(arima_sign) %>%
  step_nzv(all_numeric_predictors())

tree_juiced <- tree_rec %>% prep() %>% juice()
tree_juiced %>% glimpse()


set.seed(345)
tree_wf <- workflow() %>%
  add_model(tree_spec) %>%
  add_recipe(tree_rec)



metrics <- metric_set(kap)

entropy_grid <- grid_max_entropy(cost_complexity(),
                                 tree_depth(),
                                 min_n(),
                                 size = 50)

doParallel::registerDoParallel()
tree_res <- tune_race_anova(
  tree_wf,
  resamples = tune_boots,
  grid = entropy_grid,
  metrics = metrics,
  control = control_race(save_pred = TRUE,
                         parallel_over = "everything",
                         save_workflow = TRUE,
                         verbose_elim = TRUE)
)


tree_res


doParallel::stopImplicitCluster()


final_penalty <-
  tree_res %>%
  show_best("kap") %>% 
  slice(1)

final_penalty


final_rs <-
  tree_wf %>%
  finalize_workflow(final_penalty)

final_fit <- final_rs %>%
  fit(tbl_for_tree)


for_plot <- final_fit %>%
  extract_fit_engine()


for_plot %>%
  rpart.plot(
    branch = 1, #angle of branches
    roundint = FALSE,
    trace = 1,
    type = 1,
    extra = 9,
    tweak = 1,
    fallen.leaves = F,
    gap = 0,
    # space = 5,
    legend.y = 1,
    leaf.round = 9,
    branch.lty = 1,
    box.palette = list("#de6f6f", "grey70", "#6fde6f"),
    compress = F,
    ycompress = T, 
    main = "Predicting the response to prior day CHO\non Perceived Recovery Status", 
    cex.main=1.4,
    yesno=1
    
  )   #




resamp_fit <- final_rs %>% fit_resamples(tree_boots, control = control_resamples(save_pred = TRUE))

collect_predictions(resamp_fit, summarize = T) %>% filter(!is.na(arima_sign)) %>%  metrics(truth = arima_sign, estimate = .pred_class)
collect_predictions(resamp_fit, summarize = T) %>% conf_mat(arima_sign, .pred_class) 


tibble(
  "Negative" = c(5,0,0),
  "Non-\nsignificant" = c(1,27,3),
  "Positive" = c(4,0,0)
) %>% 
  mutate(y = c("Negative", "Non-\nsignificant", "Positive")) %>% 
  gather(x, value, c("Negative", "Non-\nsignificant", "Positive")) %>% 
  mutate(x = factor(x), # alphabetical order by default
         y = factor(y, levels = rev(unique(y))),
         z = c(1,2,2,2,1,2,2,2,1)) %>% 
  ggplot(aes(x,y, fill = z))+
  geom_tile(color="black") + 
  coord_equal() +
  scale_fill_gradient(low = "#afedaf", high = "#dff8df")+
  geom_text(aes(label=value), color="black") +
  guides(fill="none") +
  scale_x_discrete(position = "top", expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "True class", y = "Predicted class", fill = NULL)+  
  theme(
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "grey60", fill=NA, size=.8)
  )

