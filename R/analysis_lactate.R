# ============================================================
# Logistic regression and ROC analysis for lactate and
# in-hospital mortality
#
# Default input:  data/SPSSlactato.csv
# Default output: outputs/
# ============================================================

# -----------------------------
# 0) Configuration
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- ifelse(length(args) >= 1, args[[1]], file.path("data", "SPSSlactato.csv"))
output_dir <- ifelse(length(args) >= 2, args[[2]], "outputs")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

required_packages <- c("dplyr", "broom", "pROC", "ggplot2")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(broom)
library(pROC)
library(ggplot2)

# -----------------------------
# 1) Helpers
# -----------------------------
to_numeric <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

standardize_text <- function(x) {
  toupper(trimws(as.character(x)))
}

format_p <- function(p_value) {
  ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
}

format_or_ci <- function(or, low, high) {
  sprintf("%.2f (%.2f--%.2f)", or, low, high)
}

extract_logit_results <- function(model, model_name, labels) {
  coefs <- summary(model)$coefficients
  ci_wald <- confint.default(model)

  out <- data.frame(
    term = rownames(coefs),
    B = coefs[, "Estimate"],
    SE = coefs[, "Std. Error"],
    OR = exp(coefs[, "Estimate"]),
    CI_low = exp(ci_wald[, 1]),
    CI_high = exp(ci_wald[, 2]),
    p_value = coefs[, "Pr(>|z|)"],
    row.names = NULL
  )

  out %>%
    filter(term != "(Intercept)") %>%
    mutate(
      Model = model_name,
      Predictor = recode(term, !!!labels),
      OR_CI = format_or_ci(OR, CI_low, CI_high),
      p_formatted = format_p(p_value)
    ) %>%
    select(Model, Predictor, B, SE, OR, CI_low, CI_high, p_value, OR_CI, p_formatted)
}

make_roc_df <- function(roc_object, curve_name) {
  out <- data.frame(
    fpr = 1 - roc_object$specificities,
    sensitivity = roc_object$sensitivities,
    curve = curve_name
  )

  out <- out[is.finite(out$fpr) & is.finite(out$sensitivity), ]
  out <- out[
    out$fpr >= 0 & out$fpr <= 1 &
      out$sensitivity >= 0 & out$sensitivity <= 1,
  ]
  unique(out)
}

write_latex_table <- function(results, path) {
  lines <- c(
    "\\begin{table}[ht]",
    "\\centering",
    "\\caption{Logistic regression for in-hospital mortality}",
    "\\begin{tabular}{llcc}",
    "\\hline",
    "Model & Predictor & OR (95\\% CI) & p value \\\\",
    "\\hline"
  )

  for (i in seq_len(nrow(results))) {
    line <- paste0(
      results$Model[i], " & ",
      results$Predictor[i], " & ",
      results$OR_CI[i], " & ",
      results$p_formatted[i], " \\\\"
    )
    lines <- c(lines, line)
  }

  lines <- c(
    lines,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )

  writeLines(lines, con = path)
}

# -----------------------------
# 2) Read and validate data
# -----------------------------
if (!file.exists(input_file)) {
  stop(
    "Input file not found: ", input_file, "\n",
    "Place the private dataset at data/SPSSlactato.csv or pass a custom path:\n",
    "Rscript R/analysis_lactate.R path/to/SPSSlactato.csv outputs"
  )
}

df <- read.csv2(input_file, stringsAsFactors = FALSE, check.names = FALSE)

required_columns <- c("Edad", "SEXO", "MORTALIDAD", "DIA1", "DIA2")
missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

cat("\nColumns found:\n")
print(names(df))

# -----------------------------
# 3) Cleaning and coding
# -----------------------------
df <- df %>%
  mutate(
    Edad = to_numeric(Edad),
    DIA1 = to_numeric(DIA1),
    DIA2 = to_numeric(DIA2),
    SEXO = standardize_text(SEXO),
    MORTALIDAD = standardize_text(MORTALIDAD),
    death = ifelse(MORTALIDAD == "SI", 1, 0),
    sex_male = ifelse(SEXO == "MASCULINO", 1, 0),
    T2_gt3 = ifelse(DIA2 > 3, 1, 0)
  )

df_model1 <- df %>% filter(complete.cases(death, DIA1, DIA2))
df_model2 <- df %>% filter(complete.cases(death, Edad, sex_male, T2_gt3))
df_roc <- df %>% filter(complete.cases(DIA1, DIA2, death))

cat("\nOutcome coding:\n")
print(table(df$MORTALIDAD, df$death, useNA = "ifany"))

cat("\nSex coding:\n")
print(table(df$SEXO, df$sex_male, useNA = "ifany"))

cat("\nT2 > 3 mmol/L vs mortality:\n")
print(table(Mortality = df$death, T2_gt3 = df$T2_gt3, useNA = "ifany"))

# -----------------------------
# 4) Logistic regression
# -----------------------------
model1 <- glm(
  death ~ DIA1 + DIA2,
  data = df_model1,
  family = binomial(link = "logit")
)

model2 <- glm(
  death ~ Edad + sex_male + T2_gt3,
  data = df_model2,
  family = binomial(link = "logit")
)

labels_m1 <- c(
  "DIA1" = "T1 lactate (continuous)",
  "DIA2" = "T2 lactate (continuous)"
)

labels_m2 <- c(
  "Edad" = "Age (years)",
  "sex_male" = "Sex (male vs female)",
  "T2_gt3" = "T2 lactate > 3 mmol/L"
)

res_m1 <- extract_logit_results(
  model1,
  "Model 1: Lactate at admission (T1) and 24 h (T2)",
  labels_m1
)

res_m2 <- extract_logit_results(
  model2,
  "Model 2: Age, sex, and T2 lactate > 3 mmol/L",
  labels_m2
)

logistic_results <- bind_rows(res_m1, res_m2)

write.csv(
  logistic_results,
  file = file.path(output_dir, "Logistic_regression_results.csv"),
  row.names = FALSE
)

write_latex_table(
  logistic_results,
  file.path(output_dir, "Logistic_regression_table_latex.txt")
)

cat("\nLogistic regression results:\n")
print(logistic_results)

# -----------------------------
# 5) ROC analysis
# -----------------------------
roc_d1 <- roc(
  response = df_roc$death,
  predictor = df_roc$DIA1,
  levels = c(0, 1),
  direction = "<",
  ci = TRUE,
  quiet = TRUE
)

roc_d2 <- roc(
  response = df_roc$death,
  predictor = df_roc$DIA2,
  levels = c(0, 1),
  direction = "<",
  ci = TRUE,
  quiet = TRUE
)

auc_d1 <- auc(roc_d1)
ci_d1 <- ci.auc(roc_d1)

auc_d2 <- auc(roc_d2)
ci_d2 <- ci.auc(roc_d2)

delong_test <- roc.test(roc_d2, roc_d1, method = "delong")

roc_summary <- data.frame(
  Curve = c("T1", "T2", "DeLong T2 vs T1"),
  AUC = c(as.numeric(auc_d1), as.numeric(auc_d2), NA),
  CI_low = c(as.numeric(ci_d1[1]), as.numeric(ci_d2[1]), NA),
  CI_high = c(as.numeric(ci_d1[3]), as.numeric(ci_d2[3]), NA),
  p_value = c(NA, NA, delong_test$p.value)
)

write.csv(
  roc_summary,
  file = file.path(output_dir, "ROC_summary.csv"),
  row.names = FALSE
)

cat("\nROC summary:\n")
print(roc_summary)

# -----------------------------
# 6) DIA2 cutoff > 3 mmol/L
# -----------------------------
cutoff <- 3

coords_d2 <- coords(
  roc_d2,
  x = cutoff,
  input = "threshold",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

pred_d2 <- ifelse(df_roc$DIA2 > cutoff, 1, 0)

TP <- sum(pred_d2 == 1 & df_roc$death == 1)
FN <- sum(pred_d2 == 0 & df_roc$death == 1)
FP <- sum(pred_d2 == 1 & df_roc$death == 0)
TN <- sum(pred_d2 == 0 & df_roc$death == 0)

cutoff_metrics <- data.frame(
  cutoff = cutoff,
  rule = "DIA2 > 3 mmol/L",
  TP = TP,
  FN = FN,
  FP = FP,
  TN = TN,
  sensitivity = TP / (TP + FN),
  specificity = TN / (TN + FP),
  accuracy = (TP + TN) / nrow(df_roc)
)

write.csv(
  cutoff_metrics,
  file = file.path(output_dir, "DIA2_cutoff_3_confusion_matrix.csv"),
  row.names = FALSE
)

cat("\nDIA2 cutoff > 3 mmol/L:\n")
print(cutoff_metrics)

# -----------------------------
# 7) Publication ROC plot
# -----------------------------
roc_df <- rbind(
  make_roc_df(roc_d2, "T2"),
  make_roc_df(roc_d1, "T1")
)

roc_df$curve <- factor(roc_df$curve, levels = c("T2", "T1"))
roc_df <- roc_df[order(roc_df$curve, roc_df$fpr, roc_df$sensitivity), ]

stats_text <- paste0(
  "AUC (95% CI)\n",
  "T2: ", sprintf("%.3f", as.numeric(auc_d2)),
  " (", sprintf("%.3f", as.numeric(ci_d2[1])),
  "-", sprintf("%.3f", as.numeric(ci_d2[3])), ")\n",
  "T1: ", sprintf("%.3f", as.numeric(auc_d1)),
  " (", sprintf("%.3f", as.numeric(ci_d1[1])),
  "-", sprintf("%.3f", as.numeric(ci_d1[3])), ")\n",
  ifelse(
    delong_test$p.value < 0.001,
    "DeLong test: p < 0.001",
    paste0("DeLong test: p = ", sprintf("%.3f", delong_test$p.value))
  )
)

cutoff_df <- data.frame(
  fpr = 1 - as.numeric(coords_d2$specificity),
  sensitivity = as.numeric(coords_d2$sensitivity)
)

cutoff_label <- paste0(
  "T2 cutoff > ", cutoff, " mmol/L\n",
  "Sensitivity = ", sprintf("%.2f", as.numeric(coords_d2$sensitivity)), "\n",
  "Specificity = ", sprintf("%.2f", as.numeric(coords_d2$specificity))
)

roc_plot <- ggplot(
  roc_df,
  aes(x = fpr, y = sensitivity, linetype = curve)
) +
  geom_step(
    linewidth = 1.0,
    color = "black",
    direction = "hv"
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    linewidth = 0.6,
    color = "black"
  ) +
  geom_point(
    data = cutoff_df,
    aes(x = fpr, y = sensitivity),
    inherit.aes = FALSE,
    size = 2.6,
    color = "black"
  ) +
  annotate(
    "text",
    x = 0.98,
    y = 0.22,
    label = stats_text,
    hjust = 1,
    vjust = 1,
    size = 3.0,
    lineheight = 0.95
  ) +
  annotate(
    "text",
    x = min(cutoff_df$fpr + 0.12, 0.75),
    y = cutoff_df$sensitivity,
    label = cutoff_label,
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    lineheight = 0.95
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("T2" = "solid", "T1" = "longdash"),
    labels = c("T2" = "T2", "T1" = "T1"),
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "1 - Specificity",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = "Sensitivity",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_equal(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.83, 0.56),
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.25
    ),
    legend.key.width = grid::unit(1.0, "cm"),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.5, color = "black"),
    plot.margin = margin(8, 8, 8, 8)
  )

print(roc_plot)

ggsave(
  filename = file.path(output_dir, "ROC_lactate_T1_T2_publication.pdf"),
  plot = roc_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "ROC_lactate_T1_T2_publication.tiff"),
  plot = roc_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(output_dir, "ROC_lactate_T1_T2_publication.png"),
  plot = roc_plot,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 600
)

cat("\nDone. Outputs saved to: ", normalizePath(output_dir), "\n", sep = "")
