# ==============================================================================
# UCI Cryotherapy | MixedDeepMM 0.1.0 | six model configurations
# Gaussian / Student-t / skew-t x single / two heads
# 20 independent initializations per model; final paper report is Mean +/- SE.
# Main clustering: Single=common first hidden layer; Two=shared-tail layer.
# Final output: one Mean+SE CSV and a fixed-coordinate 7-panel Gower PCoA figure.
# Monte Carlo sampling is retained for estimation; no sampling-efficiency diagnostics reported.
# ==============================================================================

# 1. SETTINGS -----------------------------------------------------------------
RUN_MODE <- "fit"              # "fit" or "plot_only"
NSTART <- 20L                  # Smoke test: set to 1L
FIRST_SEED <- 42L
MAX_ITER <- 100L               # Smoke test: set to 10L
MIN_ITER <- 8L                 # Smoke test: set to 3L
MC_SAMPLES <- 48L
MC_EVAL <- 192L
TOL <- 0.002
NU_INITIAL <- 8
ESTIMATE_NU <- TRUE
K_HEAD <- 2L
K_DISCRETE <- 2L
K_TAIL <- 2L
RIDGE <- 0.001
PSI_FLOOR <- 0.03

get_script_directory <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", args, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f[[1L]]), mustWork = FALSE)))
  sources <- lapply(sys.frames(), function(fr) fr$ofile)
  sources <- Filter(function(z) is.character(z) && length(z) == 1L && nzchar(z), sources)
  if (length(sources)) return(dirname(normalizePath(tail(sources, 1L)[[1L]], mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    z <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (is.character(z) && length(z) == 1L && nzchar(z))
      return(dirname(normalizePath(z, mustWork = FALSE)))
  }
  getwd()
}
SCRIPT_DIR <- get_script_directory()
DATA_FILE <- file.path(SCRIPT_DIR, "Cryotherapy_90.csv")
if (!file.exists(DATA_FILE)) {
  alt <- file.path(getwd(), "Cryotherapy_90.csv")
  if (file.exists(alt)) {
    DATA_FILE <- normalizePath(alt)
  } else {
    stop("Cannot find Cryotherapy_90.csv; extract this ZIP into one directory.")
  }
}
OUTPUT_ROOT <- dirname(DATA_FILE)
RESULT_DIR <- file.path(OUTPUT_ROOT, "Cryotherapy_PaperAligned_6Models_YYYYMMDD_HHMMSS")
if (!(RUN_MODE %in% c("fit", "plot_only"))) stop("RUN_MODE must be 'fit' or 'plot_only'.")
OUT_DIR <- if (RUN_MODE == "fit") {
  file.path(OUTPUT_ROOT, paste0("Cryotherapy_PaperAligned_6Models_",
                                format(Sys.time(), "%Y%m%d_%H%M%S")))
} else {
  RESULT_DIR
}
if (RUN_MODE == "fit") {
  if (!dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output directory.")
  for (sub in c("models", "starts", "assignments", "posteriors", "reports", "figures", "traces")) {
    if (!dir.create(file.path(OUT_DIR, sub), showWarnings = FALSE)) stop("Cannot create subdirectory: ", sub)
  }
} else if (!dir.exists(OUT_DIR)) stop("RESULT_DIR does not exist: ", OUT_DIR)
out <- function(...) file.path(OUT_DIR, ...)
save_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "")

# 2. DATA ---------------------------------------------------------------------
continuous <- c("age", "time", "area")
types <- c(age = "continuous", time = "continuous", area = "continuous", n_warts = "count", sex_male = "binary", wart_type = "nominal")

if (RUN_MODE == "fit") {
  if (!requireNamespace("MixedDeepMM", quietly = TRUE) ||
      utils::packageVersion("MixedDeepMM") != "0.1.0")
    stop("This analysis requires MixedDeepMM 0.1.0. Install the versioned source package MixedDeepMM_0.1.0.tar.gz before fitting.")
  if (!requireNamespace("GIGrvg", quietly = TRUE))
    stop("Install GIGrvg for skew-t models: install.packages('GIGrvg').")
  raw <- utils::read.csv(DATA_FILE, check.names = FALSE, stringsAsFactors = FALSE,
                         na.strings = c("NA", ""))
  if (nrow(raw) != 90L || anyNA(raw) || anyDuplicated(raw$Source_Row))
    stop("Input CSV must contain 90 complete, uniquely indexed observations.")
  if (!("Result_of_Treatment" %in% names(raw))) stop("Missing external evaluation label.")
  source_rows <- as.integer(raw$Source_Row)
  truth <- factor(raw[["Result_of_Treatment"]], levels = c(0, 1),
                  labels = c("Class0", "Class1"))
  if (anyNA(truth) || length(unique(truth)) != 2L)
    stop("External label must contain classes 0 and 1.")
  X <- data.frame(
    age = as.numeric(raw$age),
    time = as.numeric(raw$Time),
    area = as.numeric(raw$Area),
    n_warts = as.numeric(raw$Number_of_Warts),
    sex_male = as.integer(raw$sex == 1L),
    wart_type = factor(raw$Type, levels = c(1, 2, 3))
  )
  if (!identical(names(X), names(types)) || anyNA(X)) stop("Predictor preparation failed.")
  for (nm in continuous) {
    if (!is.finite(stats::sd(X[[nm]])) || stats::sd(X[[nm]]) <= 0)
      stop("Continuous feature has zero/invalid variance: ", nm)
    X[[nm]] <- as.numeric(scale(X[[nm]]))
  }
  save_csv(data.frame(Source_Row = source_rows, Truth = truth, X, check.names = FALSE),
           out("Cryotherapy_Processed.csv"))
  save_csv(X, out("Cryotherapy_Predictors_Standardized.csv"))
  writeLines(c(
    "UCI Cryotherapy", "90 cryotherapy-treated observations", "Observations: 90",
    "Input: Cryotherapy_90.csv", "Original source: https://archive.ics.uci.edu/dataset/429/cryotherapy%2Bdataset",
    "External label: Result_of_Treatment; classes 0=Failure / 1=Success",
    "The external label is excluded from model fitting and model selection.",
    "Continuous features are z-standardized on the full evaluation subset.",
    "Input measurement types: age=continuous, time=continuous, area=continuous, n_warts=count, sex_male=binary, wart_type=nominal",
    "Cryotherapy subset of the original UCI wart data (90 patients); sex and wart type enter as categorical predictors; treatment outcome is withheld.",
    "Single-head primary clustering: common first mixture layer.",
    "Two-head primary clustering: shared-tail mixture layer.",
    "Best fit: highest average Gower silhouette width of an individual completed restart (ties by LogLik). External labels never select models.",
    "Observed treatment outcomes:", capture.output(print(table(truth)))
  ), out("Data_Description.txt"))
  rm(raw)
} else {
  processed <- out("Cryotherapy_Processed.csv")
  if (!file.exists(processed)) stop("No processed data in RESULT_DIR: ", processed)
  z <- utils::read.csv(processed, check.names = FALSE, stringsAsFactors = FALSE)
  source_rows <- as.integer(z$Source_Row)
  truth <- factor(z$Truth, levels = c("Class0", "Class1"))
  X <- z[, names(types), drop = FALSE]
  X$wart_type <- factor(X$wart_type, levels = c(1, 2, 3))
  if (anyNA(X) || anyNA(truth)) stop("Invalid processed data.")
}
stopifnot(nrow(X) == 90L, identical(names(X), names(types)))
# Gower distance is computed ONCE from predictors only: continuous/count
# numeric range-scaled by daisy, binary/nominal symmetric factors, ordinal ordered.
# Do not pass the externally withheld labels to this function.
if (!requireNamespace("cluster", quietly = TRUE)) {
  stop("Install 'cluster' to compute the mandatory Gower silhouette: install.packages('cluster').")
}
if (!requireNamespace("MixedDeepMM", quietly = TRUE) ||
    utils::packageVersion("MixedDeepMM") != "0.1.0") {
  stop("This analysis requires MixedDeepMM 0.1.0. Install the versioned source package MixedDeepMM_0.1.0.tar.gz before running this script.")
}
gower_dist <- MixedDeepMM::mixeddeepmm_gower_distance(X, types)

# 3. SIX CONFIGURATIONS ------------------------------------------------------
models <- list(
  Gaussian_Single = c(family = "gaussian", heads = "single"),
  Gaussian_Two    = c(family = "gaussian", heads = "two"),
  StudentT_Single = c(family = "t", heads = "single"),
  StudentT_Two    = c(family = "t", heads = "two"),
  SkewT_Single    = c(family = "skew-t", heads = "single"),
  SkewT_Two       = c(family = "skew-t", heads = "two")
)
model_titles <- c(
  Gaussian_Single = "DGMM (Single)", Gaussian_Two = "DGMM (Double)",
  StudentT_Single = "RDMM (Single)", StudentT_Two = "RDMM (Double)",
  SkewT_Single = "DStMM (Single)", SkewT_Two = "DStMM (Double)")
single_layers <- list(k = c(K_HEAD, K_TAIL), r = c(2L, 1L))
two_layers <- list(k = c(K_HEAD, K_DISCRETE, K_TAIL),
                   r = list(continuous = 2L, discrete = 2L, tail = 1L))

# 4. HELPERS -----------------------------------------------------------------
num <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  y <- suppressWarnings(as.numeric(x[1L]))
  if (length(y) == 1L && is.finite(y)) y else NA_real_
}
ari <- function(g) as.numeric(MixedDeepMM::mixeddeepmm_ari(truth, g))
mr <- function(g) as.numeric(MixedDeepMM::mixeddeepmm_mr(truth, g))
fmt <- function(x, digits = 3L) {
  if (length(x) == 0L || !is.finite(x[1L])) return("NA")
  sprintf(paste0("%.", digits, "f"), x[1L])
}
empty_init <- function(name, j, seed) {
  data.frame(Model = name, Initialization = j, Seed = seed,
             Status = "Fit failed", Gower_Silhouette = NA_real_,
             ARI = NA_real_, MR = NA_real_,
             Cluster1_N = NA_integer_, Cluster2_N = NA_integer_,
             Nonempty_Clusters = NA_integer_, BIC = NA_real_,
             LogLik = NA_real_, Nu_Hat = NA_real_, Iterations = NA_real_,
             Fit_Sec = NA_real_, Predict_Sec = NA_real_, Total_Sec = NA_real_,
             Error = "", stringsAsFactors = FALSE)
}
print_init <- function(r) {
  cat(sprintf("%-18s | Init %02d/%02d | %-23s | Sil %s | ARI %s | MR %s | BIC %s | nu %s | Clusters %s | Iter %s | Time %ss\n",
              r$Model, r$Initialization, NSTART, r$Status,
              fmt(r$Gower_Silhouette, 4), fmt(r$ARI, 4), fmt(r$MR, 4),
              fmt(r$BIC, 2), fmt(r$Nu_Hat, 3), fmt(r$Nonempty_Clusters, 0),
              fmt(r$Iterations, 0), fmt(r$Total_Sec, 1)))
  if (nzchar(r$Error)) cat("    Error: ", r$Error, "\n", sep = "")
  utils::flush.console()
}

# 5. INDEPENDENT FITS PER MODEL ----------------------------------------------
if (RUN_MODE == "fit") {
  all_inits <- empty_init("", 0L, 0L)[FALSE, , drop = FALSE]
  model_results <- data.frame()
  cat("\nInput:", DATA_FILE, "\nOutput:", OUT_DIR, "\n")
  cat("Six configurations x", NSTART, "single-start fits. Select an individual best fit by Gower silhouette (LogLik tie-break only); paper summaries use all valid fits.\n")
  for (name in names(models)) {
    cfg <- models[[name]]
    family <- cfg[["family"]]
    heads <- cfg[["heads"]]
    layers <- if (heads == "single") single_layers else two_layers
    per_model <- empty_init("", 0L, 0L)[FALSE, , drop = FALSE]
    best_sil <- -Inf
    best_loglik <- -Inf
    best_row <- NULL
    best_pred <- NULL
    cat("\n==================== ", name, " ====================\n", sep = "")
    for (j in seq_len(NSTART)) {
      seed <- as.integer(FIRST_SEED + j - 1L)
      r <- empty_init(name, j, seed)
      tick <- proc.time()[["elapsed"]]
      ctrl <- list(seed = seed, nstart = 1L, max_iter = MAX_ITER,
                   min_iter = MIN_ITER, tol = TOL,
                   mc_samples = MC_SAMPLES, mc_eval = MC_EVAL,
                   nu = NU_INITIAL, estimate_nu = ESTIMATE_NU,
                   ridge = RIDGE, psi_floor = PSI_FLOOR, verbose = FALSE)
      fit <- tryCatch(
        MixedDeepMM::mixeddeepmm(data = X, types = types,
                                 family = family, heads = heads,
                                 layers = layers, control = ctrl),
        error = function(e) e)
      r$Fit_Sec <- proc.time()[["elapsed"]] - tick
      if (inherits(fit, "error")) {
        r$Error <- conditionMessage(fit)
      } else {
        r$LogLik <- num(fit$loglik)
        r$BIC <- num(fit$bic)
        r$Nu_Hat <- if (family == "gaussian") NA_real_ else num(fit$params$nu)
        r$Iterations <- num(fit$iterations)
        valid <- is.finite(r$LogLik) &&
          (family == "gaussian" || (is.finite(r$Nu_Hat) && r$Nu_Hat > 2))
        if (!valid) {
          r$Status <- "Invalid fit"
          r$Error <- "Nonfinite likelihood or invalid degrees of freedom."
        } else {
          tpred <- proc.time()[["elapsed"]]
          pred <- tryCatch(stats::predict(fit, type = "layer", cluster_level = "paper"),
                           error = function(e) e)
          r$Predict_Sec <- proc.time()[["elapsed"]] - tpred
          if (inherits(pred, "error")) {
            r$Status <- "Prediction failed"
            r$Error <- conditionMessage(pred)
          } else {
            metrics <- tryCatch({
              if (!is.matrix(pred) || nrow(pred) != nrow(X) ||
                  ncol(pred) != K_HEAD || any(!is.finite(pred)) ||
                  any(abs(rowSums(pred) - 1) > 1e-3))
                stop("Invalid primary-layer posterior probabilities")
              g <- max.col(pred, ties.method = "first")
              sizes <- tabulate(g, nbins = K_HEAD)
              c(Gower_Silhouette = MixedDeepMM::mixeddeepmm_gower_silhouette(g, gower_dist),
                ARI = ari(g), MR = mr(g), Cluster1_N = sizes[1L],
                Cluster2_N = sizes[2L], Nonempty_Clusters = sum(sizes > 0L))
            }, error = function(e) e)
            if (inherits(metrics, "error")) {
              r$Status <- "Metrics failed"
              r$Error <- conditionMessage(metrics)
            } else {
              for (nm in names(metrics)) r[[nm]] <- metrics[[nm]]
              if (!is.finite(r$Gower_Silhouette)) {
                r$Status <- "Undefined silhouette"
                r$Error <- "Fewer than two occupied MAP clusters (or invalid silhouette); cannot rank this start."
              } else {
                r$Status <- "OK"
              }
            }
          }
          # Save EACH valid run's MAP partition for label-invariant 20-run consensus.
          # These checkpoint files are removed after final outputs are created.
          if (r$Status == "OK") {
            save_csv(
              data.frame(Source_Row = source_rows,
                         PaperCluster = max.col(pred, ties.method = "first")),
              out("assignments", sprintf("%s_init_%02d.csv", name, j))
            )
          }
          # Selection uses only the predictor-based Gower silhouette.
          # External truth / ARI / MR are NEVER used in selection.
          # Likelihood breaks exact numerical silhouette ties ONLY.
          if (r$Status == "OK" &&
              (r$Gower_Silhouette > best_sil + 1e-10 ||
               (abs(r$Gower_Silhouette - best_sil) <= 1e-10 &&
                r$LogLik > best_loglik))) {
            best_sil <- r$Gower_Silhouette
            best_loglik <- r$LogLik
            best_row <- r
            best_pred <- pred
            saveRDS(fit, out("models", paste0(name, "_best.rds")))
          }
        }
      }
      r$Total_Sec <- proc.time()[["elapsed"]] - tick
      if (!is.null(best_row) && best_row$Initialization == j)
        best_row$Total_Sec <- r$Total_Sec
      per_model <- rbind(per_model, r)
      all_inits <- rbind(all_inits, r)
      save_csv(per_model, out("starts", paste0(name, "_initializations.csv")))
      save_csv(all_inits, out("all_initializations.csv"))
      print_init(r)
      rm(fit)
    }
    if (is.null(best_row)) {
      selected <- empty_init(name, NA_integer_, NA_integer_)
      selected$Status <- "No valid silhouette"
    } else {
      selected <- best_row
      best_fit <- readRDS(out("models", paste0(name, "_best.rds")))
      if (length(best_fit$trace))
        save_csv(data.frame(Iteration = seq_along(best_fit$trace),
                            LogLik = as.numeric(best_fit$trace)),
                 out("traces", paste0(name, "_best.csv")))
      g <- max.col(best_pred, ties.method = "first")
      assignments <- data.frame(Source_Row = source_rows,
                                Sample_ID = seq_len(nrow(X)),
                                PaperCluster = g, True_Class = truth)
      save_csv(assignments, out("assignments", paste0(name, "_best.csv")))
      posterior <- cbind(data.frame(Source_Row = source_rows,
                                    Sample_ID = seq_len(nrow(X))),
                         stats::setNames(as.data.frame(best_pred),
                                         paste0("PaperLayer_", seq_len(ncol(best_pred)))))
      save_csv(posterior, out("posteriors", paste0(name, "_best.csv")))
      tab <- table(truth, g)
      save_csv(data.frame(True_Class = rownames(tab), as.data.frame.matrix(tab),
                          check.names = FALSE),
               out("reports", paste0(name, "_cluster_contingency.csv")))
      writeLines(c(paste("Model:", name),
                   paste("Selected initialization:", best_row$Initialization),
                   paste("Seed:", best_row$Seed),
                   paste("Selection: highest Gower silhouette of an individual completed start; LogLik for ties only."),
                   paste("Gower Silhouette:", fmt(best_row$Gower_Silhouette, 5)),
                   paste("BIC (secondary):", fmt(best_row$BIC, 5)),
                   paste("ARI:", fmt(best_row$ARI, 5)),
                   paste("MR:", fmt(best_row$MR, 5)),
                   paste("Number of nonempty clusters:", best_row$Nonempty_Clusters),
                   paste("Estimated nu:", fmt(best_row$Nu_Hat, 5)),
                   paste("Iterations:", fmt(best_row$Iterations, 0))),
                 out("reports", paste0(name, "_best.txt")))
      rm(best_fit)
    }
    selected$Successful_Inits <- sum(per_model$Status %in% c("OK", "Undefined silhouette"))
    selected$Valid_Silhouette_Inits <- sum(per_model$Status == "OK")
    selected$Undefined_Silhouette_Inits <- sum(per_model$Status == "Undefined silhouette")
    selected$Failed_Inits <- NSTART - selected$Successful_Inits
    model_results <- rbind(model_results, selected)
    save_csv(model_results, out("model_results.csv"))
    cat("BEST-GOWER ", name, " | Init ", selected$Initialization,
        " | Sil ", fmt(selected$Gower_Silhouette, 4),
        " | BIC (secondary) ", fmt(selected$BIC, 2),
        " | ARI ", fmt(selected$ARI, 4),
        " | MR ", fmt(selected$MR, 4), "\n", sep = "")
    gc(verbose = FALSE)
  }
  # Compare the six selected individual fits by Gower silhouette (descriptive only).
  # This never uses the external truth and does not refit or alter the MCEM objective.
  valid_model_rows <- which(model_results$Status == "OK" &
                            is.finite(model_results$Gower_Silhouette))
  if (length(valid_model_rows)) {
    ord_models <- valid_model_rows[order(-model_results$Gower_Silhouette[valid_model_rows])]
    leaderboard <- model_results[ord_models, c("Model", "Initialization", "Seed",
                        "Gower_Silhouette", "ARI", "MR", "BIC", "Nu_Hat",
                        "Nonempty_Clusters"), drop = FALSE]
    leaderboard$Silhouette_Rank <- seq_len(nrow(leaderboard))
    save_csv(leaderboard, out("reports", "six_model_gower_ranking.csv"))
    writeLines(c("Across-model ranking of selected individual fits: Gower silhouette, descending.",
                 "Exact ties have no statistical interpretation; no external label was used.",
                 "The first row has the maximum silhouette among the six selected fitted models.",
                 capture.output(print(leaderboard, row.names = FALSE))),
               out("reports", "six_model_gower_ranking.txt"))
  }
  writeLines(c("UCI Cryotherapy | six mixed-response model configurations",
               paste("Number of observations:", nrow(X)),
               paste("Runs per model:", NSTART),
               "Single-head paper clustering: common first latent layer.",
               "Two-head paper clustering: shared-tail latent layer.",
               "Individual best-fit selection: highest Gower silhouette among completed fits; only ties by LogLik. Paper summaries use all valid fits.",
               "This script does not select individual EM iterations by silhouette.",
               "BIC is secondary diagnostic only; neither BIC nor external labels select fits.",
               "ARI/MR use external labels for evaluation only, never for fit selection.",
               "Multiple starts are restarts on one dataset, not independent samples.",
               "Selected fits:",
               capture.output(print(model_results[, c("Model", "Status", "Initialization",
                                                     "Gower_Silhouette", "ARI", "MR", "BIC", "Nu_Hat",
                                                     "Nonempty_Clusters")], row.names = FALSE)),
               "", capture.output(utils::sessionInfo())),
             out("Final_Report.txt"))
  cat("\n===== FITTING COMPLETE =====\nResults:", OUT_DIR, "\n")
  print(model_results[, c("Model", "Status", "Initialization", "Seed",
                          "Gower_Silhouette", "ARI", "MR", "BIC", "Nu_Hat", "Nonempty_Clusters")],
        row.names = FALSE)
}

# 6. TWO-DIMENSIONAL PLOTS ----------------------------------------------------
resfile <- out("model_results.csv")
if (!file.exists(resfile)) stop("No model_results.csv: ", resfile)
model_results <- utils::read.csv(resfile, stringsAsFactors = FALSE)
labels_for_plot <- lapply(names(model_titles), function(name) {
  path <- out("assignments", paste0(name, "_best.csv"))
  if (!file.exists(path)) return(NULL)
  a <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (nrow(a) != nrow(X) ||
      !identical(as.integer(a$Source_Row), as.integer(source_rows)))
    stop("Cluster assignment row alignment error: ", path)
  a$PaperCluster
})
names(labels_for_plot) <- names(model_titles)
plot_groups <- function(xy, groups, heading, detail = NULL, x_lab, y_lab) {
  if (is.null(groups)) {
    graphics::plot.new(); graphics::title(main = paste(heading, "(unavailable)"))
    return(invisible(NULL))
  }
  g <- factor(groups)
  cols <- grDevices::hcl.colors(max(2L, nlevels(g)), "Dark 3")
  graphics::plot(xy[, 1], xy[, 2], pch = 19, cex = 0.9,
                 col = cols[as.integer(g)], xlab = x_lab, ylab = y_lab,
                 main = if (is.null(detail)) heading else paste(heading, detail, sep = "\n"))
}
draw_comparison <- function(xy, x_lab, y_lab) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(4, 2), mar = c(4.5, 5, 4.3, 1),
                mgp = c(2.6, 0.8, 0), cex = 1)
  plot_groups(xy, truth, "Observed treatment outcomes", x_lab = x_lab, y_lab = y_lab)
  for (name in names(model_titles)) {
    sel <- model_results[model_results$Model == name, , drop = FALSE]
    det <- if (nrow(sel) == 1L && sel$Status == "OK" &&
               is.finite(sel$ARI) && is.finite(sel$MR)) {
      sprintf("Sil %.3f; ARI %.3f; MR %.3f",
              sel$Gower_Silhouette, sel$ARI, sel$MR)
    } else {
      NULL
    }
    plot_groups(xy, labels_for_plot[[name]], unname(model_titles[[name]]),
                det, x_lab, y_lab)
  }
  graphics::plot.new()
  graphics::text(0.5, 0.67, "UCI Cryotherapy | n = 90 | two external classes")
  graphics::text(0.5, 0.47, paste0("Six models x ", NSTART, " separate fits"))
  graphics::text(0.5, 0.27, "Paper-layer clusters; both continuous and discrete data fitted")
}
save_figure_pair <- function(xy, prefix, x_lab, y_lab) {
  pngf <- out("figures", paste0(prefix, ".png"))
  pdff <- out("figures", paste0(prefix, ".pdf"))
  grDevices::png(pngf, width = 2600, height = 3600, res = 180, pointsize = 18)
  tryCatch(draw_comparison(xy, x_lab, y_lab), finally = grDevices::dev.off())
  grDevices::pdf(pdff, width = 2600/180, height = 3600/180, pointsize = 18)
  tryCatch(draw_comparison(xy, x_lab, y_lab), finally = grDevices::dev.off())
  cat("Saved ", pngf, " and ", pdff, "\n", sep = "")
}
xy <- stats::prcomp(X[, continuous, drop = FALSE], center = FALSE,
                    scale. = FALSE)$x[, 1:2, drop = FALSE]
save_figure_pair(xy, "PCA_6Models_Comparison", "Continuous PC1", "Continuous PC2")
# Identical cached Gower distances underpin both silhouette and PCoA.
ord <- tryCatch(stats::cmdscale(gower_dist, k = 2, eig = TRUE, add = TRUE),
                error = function(e) e)
if (!inherits(ord, "error") && is.matrix(ord$points) &&
    ncol(ord$points) >= 2L) {
  save_figure_pair(ord$points[, 1:2, drop = FALSE],
                   "Gower_PCoA_6Models_Comparison",
                   "Gower PCoA 1", "Gower PCoA 2")
} else {
  warning("Gower PCoA unavailable for these variables.")
}
cat("\n===== COMPLETE =====\nResults:", OUT_DIR, "\n")



# 7. FINAL OUTPUTS ONLY: 20-RUN MEAN +/- SE AND GOWER PCoA CONSENSUS -----------
# APPENDED TO THE ORIGINAL EXPERIMENT. No change to the six fits/MCEM objective.
# The CSV summarises ALL valid fits, not the best-Silhouette initialization.
# The plot uses label-invariant co-assignment across runs, not an arithmetic
# average of the arbitrary cluster numbers. The underlying 2D Gower PCoA
# coordinates are the SAME for the true-label panel and all six model panels.
# Intermediate artifacts are kept if the script fails; upon complete success
# in RUN_MODE='fit', they are cleaned so the output folder has ONLY 3 files.

final_metric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  c(Mean = if (n) mean(x) else NA_real_,
    SE = if (n >= 2L) stats::sd(x) / sqrt(n) else NA_real_)
}

format_mean_se <- function(mu, se, digits = 4L) {
  if (!is.finite(mu)) return(NA_character_)
  if (!is.finite(se)) return(paste0(formatC(mu, format = "f", digits = digits),
                                    " +/- NA"))
  sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), mu, se)
}

# Load one assignment per valid start, preserving the exact subject order.
load_valid_partitions <- function(model_name, selected_runs) {
  if (!nrow(selected_runs)) return(NULL)
  labels <- vector("list", nrow(selected_runs))
  for (i in seq_len(nrow(selected_runs))) {
    j <- as.integer(selected_runs$Initialization[i])
    f <- out("assignments", sprintf("%s_init_%02d.csv", model_name, j))
    if (!file.exists(f)) {
      warning("Missing per-run labels for ", model_name, " run ", j,
              "; full 20-run consensus requires the updated fit script.")
      return(NULL)
    }
    dat <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (!all(c("Source_Row", "PaperCluster") %in% names(dat)) ||
        nrow(dat) != nrow(X) ||
        !identical(as.integer(dat$Source_Row), as.integer(source_rows))) {
      stop("Per-run assignment row mismatch: ", f)
    }
    g <- as.integer(dat$PaperCluster)
    if (anyNA(g) || !all(g %in% seq_len(K_HEAD)) ||
        length(unique(g)) < 2L) stop("Invalid assignment: ", f)
    labels[[i]] <- g
  }
  do.call(cbind, labels)
}

consensus_labels <- function(M, xy) {
  if (is.null(M) || !is.matrix(M) || ncol(M) < 2L) return(NULL)
  n <- nrow(M)
  C <- matrix(0, n, n)
  for (j in seq_len(ncol(M))) {
    g <- M[, j]
    C <- C + outer(g, g, "==")
  }
  C <- C / ncol(M)
  diag(C) <- 1
  D <- 1 - C
  D[D < 0] <- 0
  D[D > 1] <- 1
  D <- (D + t(D)) / 2
  diag(D) <- 0
  
  tree <- stats::hclust(
    stats::as.dist(D),
    method = "average"
  )
  g <- as.integer(stats::cutree(tree, k = K_HEAD))
  # Canonical cluster color order (NOT an external-label alignment).
  centroid <- tapply(xy[, 1], g, mean)
  g <- match(g, order(centroid))
  as.integer(g)
}

# Determine shared Gower PCoA coordinates ONCE; values of labels do not affect xy.
final_ord <- stats::cmdscale(gower_dist, k = 2L, eig = TRUE, add = TRUE)
if (!is.matrix(final_ord$points) || ncol(final_ord$points) < 2L) {
  stop("Cannot form a two-dimensional Gower PCoA: final plot was not created.")
}
final_xy <- final_ord$points[, 1:2, drop = FALSE]
final_xy_range <- list(x = range(final_xy[, 1]), y = range(final_xy[, 2]))

# In plot_only mode, only existing result files can be used. Folders
# containing only BEST labels cannot reconstruct the all-run consensus.
final_runs <- if (RUN_MODE == "fit") {
  all_inits
} else {
  utils::read.csv(out("all_initializations.csv"), stringsAsFactors = FALSE)
}

metric_names <- c("Gower_Silhouette", "ARI", "MR", "BIC", "Nu_Hat", "Total_Sec")
need <- c("Model", "Initialization", "Status", metric_names)
if (!all(need %in% names(final_runs))) {
  stop("Missing metric columns: ", paste(setdiff(need, names(final_runs)), collapse = ", "))
}

final_partitions <- setNames(vector("list", length(models)), names(models))
final_rows <- vector("list", length(models))

for (m in seq_along(models)) {
  model_name <- names(models)[m]
  one <- final_runs[final_runs$Model == model_name, , drop = FALSE]
  valid <- one[one$Status == "OK" & is.finite(one$Gower_Silhouette), , drop = FALSE]
  rr <- data.frame(
    Model = model_name,
    Attempted = nrow(one), Valid = nrow(valid),
    Collapsed = sum(one$Status == "Undefined silhouette"),
    Failed = sum(!one$Status %in% c("OK", "Undefined silhouette")),
    stringsAsFactors = FALSE
  )
  for (nm in metric_names) {
    z <- final_metric(valid[[nm]])
    rr[[paste0(nm, "_Mean")]] <- unname(z["Mean"])
    rr[[paste0(nm, "_SE")]] <- unname(z["SE"])
    rr[[paste0(nm, "_Mean_SE")]] <- format_mean_se(z["Mean"], z["SE"],
                                                     if (nm == "BIC") 2L else 4L)
  }
  partitions <- load_valid_partitions(model_name, valid)
  g <- consensus_labels(partitions, final_xy)
  final_partitions[[model_name]] <- g
  rr$Consensus_Runs <- if (is.null(g)) 0L else ncol(partitions)
  rr$Consensus_Gower_Silhouette <- if (is.null(g)) NA_real_ else
    as.numeric(MixedDeepMM::mixeddeepmm_gower_silhouette(g, gower_dist))
  rr$Consensus_ARI <- if (is.null(g)) NA_real_ else ari(g)
  rr$Consensus_MR <- if (is.null(g)) NA_real_ else mr(g)
  rr$Consensus_Cluster1_N <- if (is.null(g)) NA_integer_ else sum(g == 1L)
  rr$Consensus_Cluster2_N <- if (is.null(g)) NA_integer_ else sum(g == 2L)
  final_rows[[m]] <- rr
}
final_table <- do.call(rbind, final_rows)
rownames(final_table) <- NULL

# Export a single comprehensive CSV, preserving both the mean run performance
# and the consensus-partition performance as DIFFERENT quantities.
final_csv_name <- "Cryotherapy_20Runs_Mean_SE_and_Consensus.csv"
final_csv <- out(final_csv_name)
save_csv(final_table, final_csv)

# Do not silently label a best-start plot as a 20-run consensus plot.
# Result folders without per-start MAP assignments cannot reconstruct consensus.
if (RUN_MODE == "plot_only" &&
    any(vapply(final_partitions, is.null, logical(1L)))) {
  stop("Mean +/- SE CSV saved, but genuine 20-run consensus figure requires ",
       "all per-start assignment CSVs. Refit once using this updated script.")
}

# One 7-panel figure: true labels + 6 model consensus partitions.
# Circle vs triangle denotes TRUE class; fill color denotes PREDICTED cluster.
# For the true-label panel, fill color instead denotes the observed class.
# Neither external labels nor ARI participate in the consensus construction.
final_colors <- c("#236EAB", "#DF883A")
final_shapes <- ifelse(truth == levels(truth)[1L], 21L, 24L)

plot_final_panels <- function() {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mfrow = c(4L, 2L), mar = c(3.7, 3.5, 3.1, 0.7),
                oma = c(0.3, 0.5, 2.0, 0.5), mgp = c(2.2, 0.7, 0),
                cex = 0.9, cex.main = 0.83)
  panel <- function(g, caption) {
    if (is.null(g)) {
      graphics::plot.new()
      graphics::title(main = paste0(caption, "\nConsensus unavailable"))
      return(invisible(NULL))
    }
    graphics::plot(final_xy[, 1L], final_xy[, 2L],
                   xlim = final_xy_range$x, ylim = final_xy_range$y,
                   pch = final_shapes, bg = final_colors[g], col = "#2F3539",
                   cex = 0.85, lwd = 0.5,
                   xlab = "Gower PCoA 1", ylab = "Gower PCoA 2", main = caption)
  }
  panel(as.integer(truth), "Observed treatment outcome (external labels)")
  for (name in names(models)) {
    z <- final_table[final_table$Model == name, , drop = FALSE]
    title <- paste0(unname(model_titles[[name]]), " | runs ", z$Valid,
                    "\nSil ", formatC(z$Gower_Silhouette_Mean, digits = 3L, format = "f"),
                    " +/- ", formatC(z$Gower_Silhouette_SE, digits = 3L, format = "f"),
                    " | ARI ", formatC(z$ARI_Mean, digits = 3L, format = "f"),
                    " +/- ", formatC(z$ARI_SE, digits = 3L, format = "f"))
    panel(final_partitions[[name]], title)
  }
  graphics::plot.new()
  graphics::legend("center",
    legend = c("Cluster / observed class 1", "Cluster / observed class 2",
               "Observed Class0: circle", "Observed Class1: triangle"),
    pt.bg = c(final_colors, "white", "white"),
    pch = c(21L, 21L, 21L, 24L),
    col = "#2F3539", pt.cex = 1.4, cex = 1.0, bty = "n")
  graphics::mtext("Cryotherapy (n=90): fixed Gower PCoA + co-assignment consensus",
                  outer = TRUE, side = 3, cex = 1.05, line = 0.4)
}

final_png_name <- "Gower_PCoA_20Runs_Consensus.png"
final_pdf_name <- "Gower_PCoA_20Runs_Consensus.pdf"
final_png <- out(final_png_name)
final_pdf <- out(final_pdf_name)
grDevices::png(filename = final_png, width = 2800, height = 3400,
               res = 200, pointsize = 17)
tryCatch(plot_final_panels(), finally = grDevices::dev.off())
grDevices::pdf(file = final_pdf, width = 14, height = 17, pointsize = 17)
tryCatch(plot_final_panels(), finally = grDevices::dev.off())

cat("\n===== FINAL MEAN +/- SE (ALL VALID RUNS) =====\n")
print(final_table[, c("Model", "Attempted", "Valid", "Collapsed", "Failed",
                       "Gower_Silhouette_Mean_SE", "ARI_Mean_SE", "MR_Mean_SE",
                       "Consensus_Runs")], row.names = FALSE)
cat("\nFinal files:\n", final_csv, "\n", final_png, "\n", final_pdf, "\n", sep = "")

# Only delete temporary files from a newly created fit folder, after the final
# CSV and figures were successfully generated. Never clean user's plot_only folder.
if (RUN_MODE == "fit" &&
    grepl("^Cryotherapy_PaperAligned_6Models_[0-9]{8}_[0-9]{6}$",
          basename(OUT_DIR)) &&
    all(file.exists(c(final_csv, final_png, final_pdf)))) {
  keep <- c(final_csv_name, final_png_name, final_pdf_name)
  extra <- list.files(OUT_DIR, all.files = TRUE, full.names = TRUE,
                      no.. = TRUE)
  extra <- extra[!(basename(extra) %in% keep)]
  if (length(extra)) unlink(extra, recursive = TRUE, force = TRUE)
  cat("Output directory cleaned: only final CSV + PNG + PDF remain.\n")
}

# ============================================================
# FINAL SIMPLE FIGURE
# True labels + 6 predicted consensus partitions
# ============================================================

stopifnot(
  exists("final_xy"),
  exists("final_partitions"),
  exists("final_png"),
  exists("final_pdf")
)

cols <- c("#2878B5", "#E69F40")

panel_names <- c(
  "Observed Treatment Outcomes",
  "DGMM Single",
  "DGMM Double",
  "RDMM Single",
  "RDMM Double",
  "DStMM Single",
  "DStMM Double"
)

panel_groups <- c(
  list(as.integer(truth)),
  unname(final_partitions[names(models)])
)

draw_simple_figure <- function() {
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  par(
    mfrow = c(4, 2),
    mar = c(3.8, 3.8, 2.8, 1),
    oma = c(0, 0, 1.5, 0),
    cex.main = 1
  )
  
  xlim <- range(final_xy[, 1])
  ylim <- range(final_xy[, 2])
  
  for (i in seq_along(panel_names)) {
    
    g <- panel_groups[[i]]
    
    if (is.null(g)) {
      plot.new()
      title(main = paste(panel_names[i], "(unavailable)"))
      next
    }
    
    plot(
      final_xy[, 1],
      final_xy[, 2],
      col = cols[as.integer(g)],
      pch = 19,
      cex = 1.5,
      xlim = xlim,
      ylim = ylim,
      xlab = "Gower PCoA 1",
      ylab = "Gower PCoA 2",
      main = panel_names[i]
    )
  }
  
  plot.new()
  legend(
    "center",
    legend = c("Group 1", "Group 2"),
    col = cols,
    pch = 19,
    pt.cex = 2, 
    cex = 1.5,     
    text.font = 2,  
    bty = "n"
  )
  
  mtext(
    "Cryotherapy: Observed Outcomes and Six Model Results",
    outer = TRUE,
    side = 3,
    line = 0,
    font = 2
  )
}

# Overwrite the original complex PNG figure
png(final_png, width = 2400, height = 3000, res = 200)
tryCatch(
  draw_simple_figure(),
  finally = dev.off()
)

# Overwrite the original complex PDF figure
pdf(final_pdf, width = 12, height = 15)
tryCatch(
  draw_simple_figure(),
  finally = dev.off()
)

cat("Simple figures saved:\n")
cat(final_png, "\n")
cat(final_pdf, "\n")
