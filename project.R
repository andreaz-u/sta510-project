# STA510 project: do ER status and positive nodes relate to death within 5 years
# in the same way in SEER and METABRIC? (notes 6.2: data -> priors -> fit -> check)
#
#   logit(mu_i) = alpha_j(i) + beta1_j(i) X1 + beta2_j(i) X2 + beta3 X3 + beta4 X4 + beta5 X5
#   X1 = ER-1, X2 = log1p(nodes)-log1p(2), X3 = age-54, X4 = size-25, X5 = grade-2
#   j(i) = dataset. M1: beta3-5 shared (9 parameters). M2: beta3-5 per dataset (12).
#   y_i = 1 if the patient died within 60 months of diagnosis.
#   delta = beta_METABRIC - beta_SEER; exp(delta) = ratio of odds ratios.
# Needs: readxl, brms, cmdstanr, ggplot2.
# Sections: 0 setup | 1 data | 2 priors and models | 3 prior predictive | 4 fit
#           5 convergence | 6 estimands | 7 posterior predictive | 8 comparison
#           9 next-100 prediction (to do)

# ---- 0. Setup ---------------------------------------------------------------
library(brms)
library(ggplot2)

PROJECT <- Sys.getenv("PROJECT", getwd())       # repository root (run from there, or set PROJECT)
dir.create(file.path(PROJECT, "figures"), showWarnings = FALSE)
HORIZON <- 60                                    # months
BACKEND <- "cmdstanr"
centre  <- c(er = 1, nodes_log = log1p(2), age = 54, size = 25, grade = 2)

odds   <- function(p) p / (1 - p)
log_or <- function(p_ref, p) log(odds(p) / odds(p_ref))
prior_from_or_range <- function(or_range, per = 1) {   # 95% range -> N(mean, sd)
  b <- log(or_range) / per
  c(mean = mean(b), sd = diff(b) / 4)
}

# ---- 1. Data ----------------------------------------------------------------
load_seer <- function() {
  raw <- read.csv(file.path(PROJECT, "data", "SEER Breast Cancer Dataset .csv"),
                  sep = ";", check.names = FALSE, strip.white = TRUE)
  names(raw) <- trimws(names(raw))
  grade <- c("Grade I" = 1, "Grade II" = 2, "Grade III" = 3, "anaplastic" = 3)
  data.frame(age = raw$Age, size = raw[["Tumor Size"]],
             grade  = unname(grade[trimws(raw$Grade)]),
             er     = c(Positive = 1, Negative = 0)[raw[["Estrogen Status"]]],
             nodes  = raw[["Reginol Node Positive"]],            # sic
             months = raw[["Survival Months"]],
             died   = c(Dead = 1, Alive = 0)[raw$Status], row.names = NULL)
}

load_metabric <- function() {
  dir <- file.path(PROJECT, "data")
  pat <- as.data.frame(readxl::read_excel(file.path(dir, "data_clinical_patient.xlsx"),
                                          na = c("", "NA")))[-(1:4), ]
  smp <- read.delim(file.path(dir, "data_clinical_sample.txt"),
                    check.names = FALSE, quote = "", comment.char = "")[-(1:4), ]
  raw <- merge(pat, smp, by = "#Patient Identifier")
  num <- function(x) suppressWarnings(as.numeric(x))
  st  <- raw[["Overall Survival Status"]]
  data.frame(age = num(raw[["Age at Diagnosis"]]), size = num(raw[["Tumor Size"]]),
             grade  = num(raw[["Neoplasm Histologic Grade"]]),
             er     = c(Positive = 1, Negative = 0)[raw[["ER Status"]]],
             nodes  = num(raw[["Lymph nodes examined positive"]]),
             months = num(raw[["Overall Survival (Months)"]]),
             died   = ifelse(is.na(st), NA, as.numeric(startsWith(as.character(st), "1"))),
             row.names = NULL)
}

# Alive with < HORIZON months of follow-up has unknown 5-year status: excluded
# (assumed missing completely at random). METABRIC is restricted to SEER-like
# patients (age 30-69, >= 1 positive node). Every exclusion is counted.
prepare <- function(df, name, restrict = FALSE) {
  n <- c(raw = nrow(df))
  df <- df[!is.na(df$months) & !is.na(df$died), ];             n["no_outcome"] <- n["raw"] - nrow(df)
  known <- df$died == 1 | df$months >= HORIZON
  n["unknown_5y"] <- sum(!known);  df <- df[known, ]
  df$y <- as.integer(df$died == 1 & df$months < HORIZON)
  k <- nrow(df);  df <- df[complete.cases(df[c("age", "size", "grade", "er", "nodes")]), ]
  n["missing_cov"] <- k - nrow(df)
  if (restrict) {
    k <- nrow(df);  df <- df[df$age >= 30 & df$age <= 69 & df$nodes >= 1, ]
    n["restricted_out"] <- k - nrow(df)
  }
  df$dataset <- name
  list(data = df, counts = c(n, analysed = nrow(df), deaths_5y = sum(df$y)))
}

seer <- prepare(load_seer(), "SEER")
mb   <- prepare(load_metabric(), "METABRIC", restrict = TRUE)
dat  <- rbind(seer$data, mb$data)
dat$dataset <- factor(dat$dataset, levels = c("SEER", "METABRIC"))
dat$X1 <- dat$er - centre["er"]                          # centred covariates, same
dat$X2 <- log1p(dat$nodes) - centre["nodes_log"]         # constants in both datasets
dat$X3 <- dat$age - centre["age"]
dat$X4 <- dat$size - centre["size"]
dat$X5 <- dat$grade - centre["grade"]

cat("Exclusions\n")
print(t(sapply(list(SEER = seer, METABRIC = mb), function(p)
  p$counts[c("raw", "no_outcome", "unknown_5y", "missing_cov", "analysed")])))
cat("\nRestricted out of METABRIC:", mb$counts["restricted_out"], "\n")
cat("\nBy dataset\n")
print(do.call(rbind, lapply(split(dat, dat$dataset), function(d)
  c(n = nrow(d), death_rate = round(mean(d$y), 3), ER_pos = round(mean(d$er), 3),
    med_nodes = median(d$nodes), med_age = median(d$age),
    med_size = median(d$size), med_grade = median(d$grade)))))
cat("\n5-year death rate by dataset x ER (0 = ER-, 1 = ER+)\n")
print(round(with(dat, tapply(y, list(dataset, er), mean)), 3))

# ---- 2. Priors and models -----------------------------------------------------
# Inputs: evidence next to each choice. Judgement (shrinkage, sd) is marked.
# alpha: NCI regional 5-yr survival 87.5% -> death 12.5%, raised to 14% (other-cause deaths)
alpha <- c(mean = qlogis(0.14), sd = 1)

# beta1 (ER+ vs ER-): Danish HR 2.08, Moffat "weak" HR < 1.5, signs flipped
er_evidence <- -log(c(danish = 2.08, moffat_weak = 1.5))
beta1 <- c(mean = -0.57, sd = 0.5)               # midpoint of the two sources

# beta2 (per log1p(nodes) unit): SEER Fig 28-3, negative vs 1-3 nodes within size bins
fig283 <- data.frame(neg = c(95, 90.5, 87.5, 88, 82.5) / 100,
                     n13 = c(87.5, 78.5, 75.5, 72, 59) / 100)
nodes_evidence <- log_or(1 - fig283$neg, 1 - fig283$n13) / log1p(2)
beta2 <- c(mean = 0.7, sd = 0.4)                 # figure mean ~0.9, shrunk

# beta3 (age/yr): OR per decade from PeerJ 7252 (40 vs 50) and PMC8693310 (60-75 vs 40-60),
# rounded to 1.0-1.5; a linear term cannot capture the U-shape at young ages
age_evidence <- c(peerj = exp(log_or(0.0762, 0.1147) / 1.7),
                  pmc_60_75 = exp(log(1.3) / 1.7), pmc_over75 = exp(log(2.5) / 3))
beta3 <- prior_from_or_range(c(1.0, 1.5), per = 10)

# beta4 (size/mm): SEER Fig 28-3, log-odds rise per mm between size bins (mean OR 1.49 per
# 10 mm); shrunk to 1.0-1.7 (unadjusted for grade/ER/age, older era)
fig_size <- 1 - rbind(neg = c(.95, .905, .875), n13 = c(.875, .785, .755), n4 = c(.725, .635, .60))
size_evidence <- cbind(log_or(fig_size[, 1], fig_size[, 2]), log_or(fig_size[, 2], fig_size[, 3])) / 10
beta4 <- prior_from_or_range(c(1.0, 1.7), per = 10)

# beta5 (grade/step): adjusted HR 2.16 (CI 1.04-4.49), BMC Cancer 2006;6:194, per step = log(HR)/2.
# sd 0.25 (CI alone gives 0.18): that study selected young women. Unadjusted SEER (1.38/step) = ceiling
grade_evidence <- c(adjusted = log(2.161) / 2, adjusted_lo = log(1.04) / 2, adjusted_hi = log(4.49) / 2,
                    unadjusted_2to3 = log_or(1 - 0.932, 1 - 0.776))
beta5 <- c(mean = 0.4, sd = 0.25)

priors <- data.frame(
  parameter = c("alpha", "beta1 ER+ vs ER-", "beta2 log1p(nodes)", "beta3 age/yr", "beta4 size/mm", "beta5 grade/step"),
  rbind(alpha, beta1, beta2, beta3, beta4, beta5), row.names = NULL)
priors$lo95 <- priors$mean - 2 * priors$sd
priors$hi95 <- priors$mean + 2 * priors$sd
priors$P_below_0 <- pnorm(0, priors$mean, priors$sd)
cat("\nPriors (alpha, beta1, beta2 per dataset; beta3-5 shared in M1)\n")
print(cbind(priors[1], round(priors[-1], 3)), row.names = FALSE)

or95 <- function(p, per = 1) round(exp(per * (p["mean"] + c(-2, 2) * p["sd"])), 2)
cat("\nalpha reference-patient risk, 95%:", round(plogis(alpha["mean"] + c(-2, 2) * alpha["sd"]), 3),
    "\nOR 95% ranges: ER+ vs ER-", or95(beta1), "| nodes/unit", or95(beta2),
    "| age/decade", or95(beta3, 10), "| size/10mm", or95(beta4, 10), "| grade/step", or95(beta5),
    "\nSd of prior dataset differences: ER", round(beta1["sd"] * sqrt(2), 2),
    "nodes", round(beta2["sd"] * sqrt(2), 2), "\n")

# One place that defines a model: formula AND brms priors, built from the objects above.
# `shared` = covariates with one slope for both datasets (M1: X3-X5; M2: none).
slope_priors <- list(X1 = beta1, X2 = beta2, X3 = beta3, X4 = beta4, X5 = beta5)

make_model <- function(shared) {
  coefs <- list(datasetSEER = alpha, datasetMETABRIC = alpha)
  terms <- "dataset"
  for (x in names(slope_priors)) {
    if (x %in% shared) {
      coefs[[x]] <- slope_priors[[x]];  terms <- c(terms, x)
    } else {
      for (d in c("SEER", "METABRIC")) coefs[[paste0("dataset", d, ":", x)]] <- slope_priors[[x]]
      terms <- c(terms, paste0("dataset:", x))
    }
  }
  list(formula = reformulate(c("0", terms), response = "y"),
       prior   = do.call(c, lapply(names(coefs), function(nm)
         prior_string(sprintf("normal(%.4f, %.4f)", coefs[[nm]]["mean"], coefs[[nm]]["sd"]),
                      class = "b", coef = nm))))
}
M1 <- make_model(shared = c("X3", "X4", "X5"))
M2 <- make_model(shared = character(0))
M0 <- make_model(shared = c("X1", "X2", "X3", "X4", "X5"))   # ER and nodes shared too: delta = 0 exactly
print(M1$formula);  print(M1$prior)                      # check the nine lines
print(M2$formula)
print(M0$formula)

# MH: partial pooling of ER and nodes. beta_j = beta_0 + u_j, u_j ~ N(0, sigma),
# sigma ~ exponential(rate). X1, X2 (population level) are the common means, with the beta1/beta2
# priors; alpha stays separate per dataset; X3-X5 shared. rate 2 (mean sigma 0.5, so delta has
# a typical sd 0.5*sqrt(2) = 0.71) matches the scale of the delta prior in M1.
make_model_h <- function(rate = 2) {
  m <- make_model(shared = c("X1", "X2", "X3", "X4", "X5"))
  list(formula = as.formula(paste(paste(deparse(m$formula), collapse = ""), "+ (0 + X1 + X2 || dataset)")),
       prior   = c(m$prior, prior_string(sprintf("exponential(%g)", rate), class = "sd")))
}
MH <- make_model_h(rate = 2)
print(MH$formula)

# ---- 3. Prior predictive check --------------------------------------
# brms with the likelihood switched off: the draws come from the priors alone.
rate_by_dataset <- function(mu) cbind(   # mu: draws x patients -> overall rate per draw
  SEER     = rowMeans(mu[, dat$dataset == "SEER"]),
  METABRIC = rowMeans(mu[, dat$dataset == "METABRIC"]))

fit_model <- function(m, file, prior_only = FALSE, ...)
  brm(m$formula, data = dat, family = bernoulli(), prior = m$prior,
      sample_prior = if (prior_only) "only" else "no",
      backend = BACKEND, chains = 4, iter = 2000, seed = 1, file = file, ...)

fit_prior <- fit_model(M1, "fit_prior", prior_only = TRUE)
cat("\nPrior predictive with brms: overall 5-year death rate (2.5%, 50%, 97.5%)\n")
print(round(apply(rate_by_dataset(posterior_epred(fit_prior)), 2, quantile, c(0.025, 0.5, 0.975)), 3))

# ---- 4. Fit -------------------------------------------------------------------------------
fit_m1 <- fit_model(M1, "fit_m1")        
fit_m2 <- fit_model(M2, "fit_m2")
fit_m0 <- fit_model(M0, "fit_m0")        # only alpha differs by dataset
fit_mh <- fit_model(MH, "fit_mh_099", control = list(adapt_delta = 0.99))   # 0.95 gave 9 divergences


# ---- 5. Convergence ------------------------------------------
check_fit <- function(fit) {
  print(summary(fit))                    # Rhat < 1.01, Bulk_ESS and Tail_ESS in the hundreds+
  np <- nuts_params(fit)
  cat("divergent transitions:", sum(np$Value[np$Parameter == "divergent__"]), "\n")
}
save_traces <- function(fit, file) {     # one PDF with all trace plots
  pdf(file.path(PROJECT, "figures", file), width = 7, height = 9);  plot(fit, ask = FALSE, newpage = FALSE);  dev.off()
}
check_fit(fit_m1)
check_fit(fit_m2)
check_fit(fit_mh)

save_traces(fit_m1, "traces_m1.pdf")
save_traces(fit_mh, "traces_mh.pdf")

# ---- 6. Estimands: dataset differences delta = beta_METABRIC - beta_SEER -------------------
delta_draws <- function(fit, x) {
  d <- as_draws_df(fit)
  d[[paste0("b_datasetMETABRIC:", x)]] - d[[paste0("b_datasetSEER:", x)]]
}
summarise_delta <- function(delta) {
  q <- quantile(delta, c(0.025, 0.5, 0.975), names = FALSE)
  c(q2.5 = q[1], median = q[2], q97.5 = q[3], P_gt_0 = mean(delta > 0),
    OR_ratio = exp(q[2]), OR_lo = exp(q[1]), OR_hi = exp(q[3]))
}
delta_draws_h <- function(fit, x) {      # MH: difference of the two group-level deviations
  d <- as_draws_df(fit)
  d[[paste0("r_dataset[METABRIC,", x, "]")]] - d[[paste0("r_dataset[SEER,", x, "]")]]
}
sigma_check <- function(fit, rate = 2) { # does the data inform sigma? prior vs posterior
  d <- as_draws_df(fit);  qs <- c(0.025, 0.5, 0.975)
  round(rbind(prior = quantile(rexp(nrow(d), rate), qs),
              post_ER = quantile(d$sd_dataset__X1, qs), post_nodes = quantile(d$sd_dataset__X2, qs)), 3)
}
delta_table <- function(fit, delta_fun = delta_draws)   # ER (X1) and nodes (X2), per model
  round(t(sapply(c(ER = "X1", nodes = "X2"), function(x) summarise_delta(delta_fun(fit, x)))), 3)

cat("\nM1: delta (METABRIC - SEER) and ratio of odds ratios\n");  print(delta_table(fit_m1))
cat("\nM2:\n");                                                    print(delta_table(fit_m2))
cat("\nMH (partial pooling, sigma ~ exponential(2)):\n");          print(delta_table(fit_mh, delta_draws_h))
cat("\nsigma: prior vs posterior (does the data inform it?)\n");    print(sigma_check(fit_mh))

# ---- 7. Posterior predictive checks (observed vs replicated death rates by group) -------------
# A failure = observed rate in the far tail of the replications (p_upper < 0.025 or > 0.975)
# in a coherent pattern. Check on the covariates whose slopes are shared in M1.
ppc_rates <- function(yrep, y, group) {
  groups <- sort(unique(as.character(group)))
  list(rep = sapply(groups, function(g) rowMeans(yrep[, group == g, drop = FALSE])),
       obs = sapply(groups, function(g) mean(y[group == g])))
}
summarise_ppc <- function(res) data.frame(
  group    = names(res$obs),  observed = res$obs,
  rep_2.5  = apply(res$rep, 2, quantile, 0.025),
  rep_50   = apply(res$rep, 2, quantile, 0.5),
  rep_97.5 = apply(res$rep, 2, quantile, 0.975),
  p_upper  = sapply(names(res$obs), function(g) mean(res$rep[, g] >= res$obs[g])),  # replications >= observed
  row.names = NULL)
plot_ppc <- function(res, title = "") {
  long <- data.frame(group = rep(colnames(res$rep), each = nrow(res$rep)), rate = as.vector(res$rep))
  obs  <- data.frame(group = names(res$obs), rate = res$obs)
  ggplot(long, aes(rate)) + geom_histogram(bins = 40) +
    geom_vline(data = obs, aes(xintercept = rate), colour = "red", linewidth = 0.8) +
    facet_wrap(~ group, scales = "free") +
    labs(x = "5-year death rate", y = "replications", title = paste("Posterior predictive check:", title, "(red = observed)"))
}
ppc_check <- function(yrep, group, title) {
  cat("\n=====", title, "=====\n");  print(table(group))
  res <- ppc_rates(yrep, dat$y, group)
  print(summarise_ppc(res), digits = 3)
  plot_ppc(res, title)                   # returns the plot: print it, or save with ggsave()
}

groupings <- list(
  ER    = paste(dat$dataset, ifelse(dat$er == 1, "ER+", "ER-")),
  nodes = paste(dat$dataset, cut(dat$nodes, c(0, 1, 3, 6, Inf), labels = c("1", "2-3", "4-6", "7+"))),
  grade = paste(dat$dataset, "grade", dat$grade),
  size  = paste(dat$dataset, cut(dat$size, c(0, 15, 25, 40, Inf), labels = c("<=15", "16-25", "26-40", ">40"))),
  age   = paste(dat$dataset, cut(dat$age, c(29, 44, 54, 64, 69), labels = c("30-44", "45-54", "55-64", "65-69"))))

yrep_m1 <- posterior_predict(fit_m1)     # draws x patients, simulated 0/1 outcomes
yrep_m2 <- posterior_predict(fit_m2)
# one at a time, e.g.:   ppc_check(yrep_m1, groupings$age, "age, M1")
#                        ppc_check(yrep_m2, groupings$age, "age, M2")

# Figures for the slides: readable font size, one PNG per check
save_ppc <- function(p, file, width = 10, height = 5) {
  ggsave(file.path(PROJECT, "figures", file), p + theme_minimal(base_size = 16), width = width, height = height, dpi = 200)
}
save_ppc(ppc_check(yrep_m1, groupings$size, "size, M1"), "ppc_size.png")
save_ppc(ppc_check(yrep_m1, groupings$age,  "age, M1"),  "ppc_age.png")

# Same check for the covariates that passed (ER, nodes, grade): prints the table, saves one PNG each
for (v in c("ER", "nodes", "grade")) {
  save_ppc(ppc_check(yrep_m1, groupings[[v]], paste(v, "M1")), paste0("ppc_", tolower(v), ".png"))
}

# M0 (ER and nodes shared across datasets): does the "no difference" model still reproduce the rates?
yrep_m0 <- posterior_predict(fit_m0)
save_ppc(ppc_check(yrep_m0, groupings$ER,    "ER, M0"),    "ppc_er_m0.png")
save_ppc(ppc_check(yrep_m0, groupings$nodes, "nodes, M0"), "ppc_nodes_m0.png")

# ---- 8. Model comparison ---------------------------------------------------------------------
# elpd_diff vs se_diff: ~2 se or more = a real difference; < 1 se = prefer the simpler M1
fit_m1 <- add_criterion(fit_m1, "loo")
fit_m2 <- add_criterion(fit_m2, "loo")
fit_m0 <- add_criterion(fit_m0, "loo")
fit_mh <- add_criterion(fit_mh, "loo")
loo_compare(fit_m1, fit_m2)

# Do ER and nodes need dataset-specific effects at all? M0 (shared, delta = 0) vs M1:
loo_compare(fit_m0, fit_m1)

# Complete pooling (M0) vs partial (MH) vs no pooling (M1) for ER and nodes:
loo_compare(fit_m0, fit_mh, fit_m1)