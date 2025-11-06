# main analysis doc

#Sam Code
required <- c("tidyverse")
new <- setdiff(required, rownames(installed.packages()))
if (length(new)) install.packages(new)
library(tidyverse)

raw <- read.csv("data/biomarker-raw.csv", check.names = FALSE)


# Peek at the columns
str(raw, max.level = 1)

# ---- Choose a few protein columns
# Heuristic: take numeric columns that look like protein intensities.
num_cols <- names(raw)[sapply(raw, is.numeric)]

# If you know specific protein names, replace sample_proteins with those names.
set.seed(1)
sample_proteins <- head(num_cols, 6)  # or sample(num_cols, 6)

# Keep ID/group columns if you have them (optional)
df <- raw[, sample_proteins, drop = FALSE] %>%
  mutate(row_id = row_number())

# ---- Long format for plotting
long_raw <- df %>%
  pivot_longer(-row_id, names_to = "protein", values_to = "value_raw")

# ---- Make a safe log-transform
# Use log1p (log(1+x)) to handle zeros; if negatives exist, shift by a small constant.
min_val <- min(long_raw$value_raw, na.rm = TRUE)
offset  <- ifelse(min_val <= 0, abs(min_val) + 1e-6, 0)

long <- long_raw %>%
  mutate(
    value_shifted = value_raw + offset,
    value_log = log1p(value_shifted)  # log(1 + shifted)
  )


# ---- Plot: RAW
p_raw <- ggplot(long, aes(x = value_raw)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ protein, scales = "free") +
  labs(
    title = "Raw protein distributions",
    x = "Raw intensity", y = "Count",
    caption = if (offset > 0) paste("Shifted by", signif(offset, 3), "for log check") else NULL
  ) +
  theme_minimal()

# ---- Plot: LOG
p_log <- ggplot(long, aes(x = value_log)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ protein, scales = "free") +
  labs(
    title = "Log-transformed protein distributions (log1p)",
    x = "log1p(intensity)", y = "Count"
  ) +
  theme_minimal()

# ---- Save figures

p_raw 
p_log


#Phillip Code
library(tidyverse)

# Get variable names
var_names <- read_csv('data/biomarker-raw.csv', 
                      col_names = F, 
                      n_max = 2, 
                      col_select = -(1:2)) %>%
  t() %>%
  as_tibble() %>%
  rename(name = V1, 
         abbreviation = V2) %>%
  na.omit()

# Read in data WITHOUT trimming
biomarker_no_trim <- read_csv('data/biomarker-raw.csv', 
                              skip = 2,
                              col_select = -2L,
                              col_names = c('group', 
                                            'empty',
                                            pull(var_names, abbreviation),
                                            'ados'),
                              na = c('-', '')) %>%
  filter(!is.na(group)) %>%
  # Log transform, center and scale (NO TRIMMING)
  mutate(across(.cols = -c(group, ados), 
                ~ scale(log10(.x))[, 1])) %>%
  select(group, ados, everything()) %>%
  # Add subject ID for tracking
  mutate(subject_id = row_number(), .before = 1)

outlier_threshold <- 3

# Count outliers per subject
outlier_counts <- biomarker_no_trim %>%
  pivot_longer(cols = -c(subject_id, group, ados),
               names_to = "biomarker",
               values_to = "z_score") %>%
  mutate(is_outlier = abs(z_score) > outlier_threshold) %>%
  group_by(subject_id, group) %>%
  summarise(
    n_outliers = sum(is_outlier, na.rm = TRUE),
    n_biomarkers_measured = sum(!is.na(z_score)),
    outlier_proportion = n_outliers / n_biomarkers_measured,
    .groups = 'drop'
  ) %>%
  arrange(desc(n_outliers))

outlier_counts %>% 
  head(10) %>%
  knitr::kable(caption = "Top 10 subjects with most outlying values")

group_summary <- outlier_counts %>%
  group_by(group) %>%
  summarise(
    n_subjects = n(),
    total_outliers = sum(n_outliers),
    mean_outliers_per_subject = mean(n_outliers),
    median_outliers_per_subject = median(n_outliers),
    subjects_with_outliers = sum(n_outliers > 0),
    .groups = 'drop'
  )

group_summary %>%
  knitr::kable(caption = "Outlier summary by group", digits = 2)

ggplot(outlier_counts, aes(x = n_outliers, fill = group)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
  labs(title = "Distribution of Outlier Counts per Subject",
       x = "Number of Outlying Biomarkers",
       y = "Number of Subjects",
       fill = "Group") +
  theme_minimal()

ggplot(outlier_counts, aes(x = group, y = n_outliers, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  labs(title = "Outlier Counts by Group",
       x = "Group",
       y = "Number of Outlying Biomarkers per Subject") +
  theme_minimal() +
  theme(legend.position = "none")


biomarker_outliers <- biomarker_no_trim %>%
  pivot_longer(cols = -c(subject_id, group, ados),
               names_to = "biomarker",
               values_to = "z_score") %>%
  mutate(is_outlier = abs(z_score) > outlier_threshold) %>%
  group_by(biomarker, group) %>%
  summarise(
    n_outliers = sum(is_outlier, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_outliers))

biomarker_outliers %>%
  head(15) %>%
  knitr::kable(caption = "Top 15 biomarkers with most outliers by group")


# Aidan Code
library(tidyverse)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)
library(pROC)

set.seed(133233)

clean_rdata <- "data/biomarker-clean.RData"
results_dir <- "results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

top_k <- 10
learn_threshold <- TRUE

# Loading our data
load(clean_rdata)
biomarker_clean <- biomarker_clean %>%
  mutate(group = factor(group, levels = c("ASD","TD")))

# Stratified train/test split ----
split_obj <- initial_split(biomarker_clean, prop = 0.8, strata = group)
train_df <- training(split_obj)
test_df <- testing(split_obj)

# Training

# t-tests
test_fn <- function(.df){
  t_test(.df,
         formula = level ~ group,
         order = c("ASD","TD"),
         alternative = "two-sided",
         var.equal = FALSE)
}

ttests_out_train <- train_df %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = purrr::map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(
    m = n(),
    hm = log(m) + 1/(2*m) - digamma(1),
    rank = row_number(),
    p.adj = m*hm*p_value/rank
  )

proteins_s1_train <- ttests_out_train %>%
  slice_min(p.adj, n = top_k) %>%
  pull(protein)

write_csv(ttests_out_train %>% select(protein, statistic, p_value, p.adj),
          file.path(results_dir, "train_ttests_all.csv"))

# Random Forest
predictors_train <- train_df %>% 
  select(-c(group, ados))
response_train   <- train_df %>% 
  pull(group) %>% factor()

rf_out_train <- randomForest(
  x = predictors_train,
  y = response_train,
  ntree = 1000,
  importance = TRUE
)

rf_imp_train <- rf_out_train$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out_train$importance)) %>%
  arrange(desc(MeanDecreaseGini))

proteins_s2_train <- rf_imp_train %>%
  slice_head(n = top_k) %>%
  pull(protein)

write_csv(rf_imp_train, file.path(results_dir, "train_rf_importance_all.csv"))

# Intersection information setup
proteins_panel <- intersect(proteins_s1_train, proteins_s2_train)
if (length(proteins_panel) < 2) {
  message("Intersection < 2; using union as a fallback.")
  proteins_panel <- union(proteins_s1_train, proteins_s2_train) %>% 
    unique()
}

tibble(source = c(rep("t_test_topk", length(proteins_s1_train)),
                  rep("rf_topk", length(proteins_s2_train))),
       protein = c(proteins_s1_train, proteins_s2_train)) %>%
  mutate(in_intersection = protein %in% proteins_panel) %>%
  write_csv(file.path(results_dir, "train_selected_proteins.csv"))

# Fitting our models with training data

train_panel <- train_df %>%
  select(group, any_of(proteins_panel)) %>%
  mutate(class = factor(group == "ASD", levels = c(FALSE, TRUE))) %>%
  select(-group)

logit_fit <- glm(class ~ ., data = train_panel, family = binomial())

# Learn threshold for our training data
if (learn_threshold) {
  train_probs <- predict(logit_fit, type = "response")
  roc_obj <- pROC::roc(response = train_panel$class, predictor = as.numeric(train_probs))
  best_coords <- pROC::coords(roc_obj, x = "best", best.method = "youden",
                              ret = c("threshold","sensitivity","specificity","auc"))
  thr <- as.numeric(best_coords["threshold"])
} else {
  thr <- 0.5
}

# Results for testing splits

test_panel <- test_df %>%
  select(group, any_of(proteins_panel)) %>%
  mutate(class = factor(group == "ASD", levels = c(FALSE, TRUE))) %>%
  select(-group)

test_probs <- predict(logit_fit, newdata = test_panel, type = "response")
test_pred <- factor(ifelse(test_probs > thr, TRUE, FALSE), levels = c(FALSE, TRUE))

# Metrics
test_metrics <- bind_rows(
  tibble(.metric = "accuracy", .est = accuracy_vec(truth = test_panel$class, estimate = test_pred)),
  tibble(.metric = "sensitivity", .est = sensitivity_vec(truth = test_panel$class, estimate = test_pred, event_level = "second")),
  tibble(.metric = "specificity", .est = specificity_vec(truth = test_panel$class, estimate = test_pred, event_level = "second")),
  tibble(.metric = "roc_auc", .est = roc_auc_vec(truth = test_panel$class, estimate = as.numeric(test_probs), event_level = "second"))
)

write_csv(test_metrics, file.path(results_dir, "test_metrics_panel.csv"))

# Confusion matrix of data
cm <- yardstick::conf_mat(
  tibble(truth = test_panel$class, estimate = test_pred),
  truth = truth, estimate = estimate
)
cm_tbl <- as_tibble(cm$table) %>% rename(Actual = Truth, Predicted = Prediction)
write_csv(cm_tbl, file.path(results_dir, "test_confusion_matrix.csv"))

print(test_metrics)
print(cm$table)

# Putting results into report
report_md <- file.path(results_dir, "bullet1_train_test_first_summary.md")
cat(
  "# Bullet 1: Train/Test-First Analysis\n\n",
  "## Panel (selected on training only)\n\n",
  paste0("- k_ttest = ", top_k, ", k_rf = ", top_k, "\n"),
  paste0("- Intersection size = ", length(proteins_panel), "\n"),
  paste0("- Proteins: ", paste(proteins_panel, collapse = ", "), "\n\n"),
  "## Test Set Metrics\n\n",
  paste0(readr::format_csv(test_metrics, na = "", eol = "\n")),
  "\n\n## Confusion Matrix (Test)\n\n",
  paste0(readr::format_csv(cm_tbl, na = "", eol = "\n")),
  file = report_md
)

# Best protein results

ttop <- ttests_out_train %>%
  arrange(p.adj) %>%
  slice_head(n = top_k) %>%
  transmute(rank_t = row_number(),
            protein,
            statistic, p_value, p.adj)

write_csv(ttop, file.path(results_dir, "proteins_ttest_topk.csv"))

cat("\nTop-k from t-tests (training):\n")
print(ttop %>% select(rank_t, protein))

# Top-k from RF importance (training)
rf_top <- rf_imp_train %>%
  slice_head(n = top_k) %>%
  transmute(rank_rf = row_number(),
            protein,
            MeanDecreaseGini)

write_csv(rf_top, file.path(results_dir, "proteins_rf_topk.csv"))

cat("\nTop-k from RF importance (training):\n")
print(rf_top %>% select(rank_rf, protein))

# Intersection panel
panel_tbl <- tibble(protein = proteins_panel) %>%
  left_join(ttop %>% select(protein, rank_t), by = "protein") %>%
  left_join(rf_top %>% select(protein, rank_rf), by = "protein") %>%
  arrange(coalesce(rank_t, Inf), coalesce(rank_rf, Inf))

write_csv(panel_tbl, file.path(results_dir, "proteins_panel_intersection.csv"))

cat("\nFinal panel (intersection) used for the logistic model:\n")
print(panel_tbl)


# Oscar Code
library(tidyverse)
library(broom)
library(pROC)
library(here)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)



set.seed(10292025)



rdata_path <- here::here("data", "biomarker-clean.RData")
objs <- load(rdata_path)

biomarker_clean



bio <- biomarker_clean


idx_td  <- which(bio$group == "TD")
idx_asd <- which(bio$group == "ASD")

split_idx <- function(ix, prop = 0.7) sample(ix, floor(length(ix)*prop))
tr_idx <- c(split_idx(idx_td), split_idx(idx_asd))
te_idx <- setdiff(seq_len(nrow(bio)), tr_idx)

bio_tr <- bio[tr_idx, , drop = FALSE]
bio_te <- bio[te_idx, , drop = FALSE]

table(Train = bio_tr$group)
table(Test  = bio_te$group)



## MULTIPLE TESTING
####################

# function to compute tests
test_fn <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out <- biomarker_clean %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1 <- ttests_out %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

proteins_s1_large <- ttests_out %>%
  slice_min(p.adj, n = 15) %>%
  pull(protein)

predictors <- biomarker_clean %>%
  select(-c(group, ados))

response <- biomarker_clean %>% pull(group) %>% factor()

# fit RF
set.seed(101422)
rf_out <- randomForest(x = predictors, 
                       y = response, 
                       ntree = 1000, 
                       importance = T)

# check errors
rf_out$confusion

# compute importance scores
proteins_s2 <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

proteins_s2_large <- rf_out$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 15) %>%
  pull(protein)


# select subset of interest
proteins_sstar <- intersect(proteins_s1, proteins_s2)

proteins_sstar_large <- intersect(proteins_s1_large, proteins_s2_large)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = factor(if_else(group == "ASD", "ASD", "TD"),
                        levels = c("TD","ASD"))) %>% #factorized group
  select(-group)

biomarker_sstar_large <- biomarker_clean %>%
  select(group, any_of(proteins_sstar_large)) %>%
  mutate(class = factor(if_else(group == "ASD", "ASD", "TD"),
                        levels = c("TD","ASD"))) %>% #factorized group
  select(-group)

# partition into training and test set
set.seed(101422)
biomarker_split <- biomarker_sstar %>%
  initial_split(prop = 0.8)

biomarker_split_large <- biomarker_sstar_large %>%
  initial_split(prop = 0.8)

#added formula to backtick column names because of column name errors when getting results
mk_formula <- function(df, resp = "class") {
  preds <- setdiff(names(df), resp)
  as.formula(paste(resp, "~", paste(sprintf("`%s`", preds), collapse = " + ")))
}

# fit logistic regression model to training set
fit <- glm(
  formula = mk_formula(training(biomarker_split), resp = "class"),
  data = training(biomarker_split), 
  family = 'binomial')

fit_large <- glm(
  formula = mk_formula(training(biomarker_split_large), resp = "class"),
  data = training(biomarker_split_large),
  family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

normal_results <- testing(biomarker_split) %>%
  add_predictions(fit, type = "response") %>%
  mutate(
    pred_prob  = pred,
    pred_class = factor(if_else(pred_prob > 0.5, "ASD", "TD"),
                        levels = c("TD","ASD"))
  ) # adjusted to allow for factorized class

large_n_results <- testing(biomarker_split_large) %>%
  add_predictions(fit_large, type = 'response') %>%
  mutate(
    pred_prob  = pred,
    pred_class = factor(if_else(pred_prob > 0.5, "ASD", "TD"),
                        levels = c("TD","ASD"))
  ) # adjusted to allow for factorized class

normal_results

large_n_results



asd_df <- biomarker_clean %>%
  filter(group == "ASD") %>%
  select(group, ados, everything())

prot_cols <- setdiff(names(asd_df), c("group","ados"))

cors <- map_dbl(prot_cols, ~ suppressWarnings(
  cor(asd_df[[.x]], asd_df$ados, use = "pairwise.complete.obs")
))

proteins_s3 <- tibble(protein = prot_cols, score = abs(cors)) %>%
  arrange(desc(score)) %>%
  slice_head(n = 10) %>%
  pull(protein)

proteins_s3


proteins_fuzzy <- union(
  union(intersect(proteins_s1, proteins_s2), intersect(proteins_s1, proteins_s3)),
  intersect(proteins_s2, proteins_s3)
)

proteins_fuzzy

# Eitan Code


## ---- part4-lasso-simple, message=FALSE, warning=FALSE ------------------
library(glmnet)
library(pROC)
library(caret)

# Load data (this file is one folder up, inside 'data/')
load("data/biomarker-clean.RData")

set.seed(133233)



idx  <- caret::createDataPartition(biomarker_clean$group, p = 0.8, list = FALSE)
train <- biomarker_clean[idx, ]
test  <- biomarker_clean[-idx, ]

train$group <- factor(train$group, levels = c("TD","ASD"))
test$group  <- factor(test$group,  levels = c("TD","ASD"))

x_train <- model.matrix(group ~ . - 1, data = subset(train, select = -ados))
y_train <- train$group
x_test  <- model.matrix(group ~ . - 1, data = subset(test,  select = -ados))
y_test  <- test$group


set.seed(133233)

cvfit <- cv.glmnet(
  x = x_train, y = y_train,
  family = "binomial",
  type.measure = "auc",
  nfolds = 8
)

coef_1se <- coef(cvfit, s = cvfit$lambda.1se)
selected <- rownames(coef_1se)[as.numeric(coef_1se) != 0]
selected <- setdiff(selected, "(Intercept)")

cat("Simpler panel size:", length(selected), "\n")
print(selected)


pred <- as.numeric(predict(cvfit, newx = x_test, s = cvfit$lambda.1se, type = "response"))
roc_obj <- roc(response = y_test, predictor = pred, levels = c("TD","ASD"), quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
cat(sprintf("Test AUROC (simpler λ_1SE model): %.3f\n", auc_val))
