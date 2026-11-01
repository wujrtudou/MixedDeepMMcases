# ============================================================================
# UCI Wine | MixedDeepMM 0.2.1 | continuous-data DGMM / RDMM / DStMM
# 20 repeated fits PER model, ONE initialization per fit (nstart = 1).
# No K-means comparison. No Gower distance/silhouette (all predictors numeric).
# Final result folder contains ONLY:
#   Wine_20Runs_Mean_SE.csv
#   Wine_PCA_20Runs_Consensus.png
#   Wine_PCA_20Runs_Consensus.pdf
# PCA plots use 20-run CO-ASSIGNMENT consensus, never averaged cluster IDs.
# ============================================================================

# 1. SETTINGS -----------------------------------------------------------------
N_RUNS <- 20L
FIRST_SEED <- 42L
NSTART <- 1L       # Fixed: ONE internal start per run (20 outer repetitions)
MAX_ITER <- 100L
MIN_ITER <- 8L
TOL <- 0.002
MC_SAMPLES <- 48L # Used for skew-t estimation
NU_INITIAL <- 8
ESTIMATE_NU <- TRUE
K_LAYERS <- c(3L, 2L)
R_LAYERS <- c(4L, 2L)
SHOW_PLOTS_IN_RSTUDIO <- TRUE

# 2. FILES AND DEPENDENCIES ---------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  zz <- grep("^--file=", args, value = TRUE)
  if (length(zz)) {
    return(dirname(normalizePath(sub("^--file=", "", zz[[1L]]),
                                 mustWork = FALSE)))
  }
  sources <- lapply(sys.frames(), function(fr) fr$ofile)
  sources <- Filter(function(z) is.character(z) && length(z) == 1L &&
                      nzchar(z), sources)
  if (length(sources)) {
    return(dirname(normalizePath(tail(sources, 1L)[[1L]],
                                 mustWork = FALSE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    z <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (is.character(z) && length(z) == 1L && nzchar(z)) {
      return(dirname(normalizePath(z, mustWork = FALSE)))
    }
  }
  getwd()
}

if (!requireNamespace("MixedDeepMM", quietly = TRUE) ||
    utils::packageVersion("MixedDeepMM") < "0.2.1") {
  stop("Please install the bundled MixedDeepMM_0.2.1.tar.gz package.")
}

for (pkg in c("mvtnorm", "corpcor", "GIGrvg")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing dependency: ", pkg,
         ". Install using install.packages('", pkg, "').")
  }
}

SCRIPT_DIR <- get_script_dir()
DATA_FILE <- file.path(SCRIPT_DIR, "wine_uci.csv")
if (!file.exists(DATA_FILE)) {
  alt <- file.path(getwd(), "wine_uci.csv")
  if (file.exists(alt)) {
    DATA_FILE <- alt
  } else {
    stop("Cannot find wine_uci.csv; extract the entire ZIP into one folder.")
  }
}

if (NSTART != 1L) {
  stop("NSTART must remain 1: the experiment uses 20 separate single-start fits.")
}
if (N_RUNS != 20L) warning("Formal Wine analysis was specified for 20 runs.")

raw <- utils::read.csv(DATA_FILE, check.names = FALSE)
if (!"Class" %in% names(raw)) stop("CSV has no external Class label.")
truth <- factor(raw$Class, levels = c(1, 2, 3))
features <- raw[, setdiff(names(raw), "Class"), drop = FALSE]
if (nrow(features) != 178L || ncol(features) != 13L ||
    anyNA(raw) || anyNA(truth) ||
    !all(vapply(features, is.numeric, logical(1))) ||
    any(!is.finite(as.matrix(features)))) {
  stop("UCI Wine must contain 178 complete rows, 13 numerical inputs and Class.")
}
if (any(vapply(features, stats::sd, numeric(1)) <= 0)) {
  stop("A numerical input has zero variance.")
}

X <- as.data.frame(scale(features))
stopifnot(nrow(X) == 178L, ncol(X) == 13L)
layers <- list(k = K_LAYERS, r = R_LAYERS)
n <- nrow(X)
proj <- stats::prcomp(X, center = FALSE, scale. = FALSE)$x[, 1:2, drop = FALSE]

root <- dirname(normalizePath(DATA_FILE))
OUT_DIR <- file.path(root, paste0("Wine_20Runs_Final_",
                                  format(Sys.time(), "%Y%m%d_%H%M%S")))
if (!dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)) {
  stop("Cannot create output directory: ", OUT_DIR)
}

# Output folder gets only the 1 CSV + 2 figure files at completion.
# All per-run estimates and cluster vectors are retained in this R session.
models <- c(DGMM = "gaussian", RDMM = "t", DStMM = "skew-t")
run_rows <- list()
assignments <- setNames(vector("list", length(models)), names(models))

safe_stat <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  nn <- length(x)
  c(n = nn,
    mean = if (nn) mean(x) else NA_real_,
    sd = if (nn > 1L) stats::sd(x) else NA_real_,
    se = if (nn > 1L) stats::sd(x) / sqrt(nn) else NA_real_)
}

get_num <- function(x) {
  if (is.null(x) || !length(x)) return(NA_real_)
  x <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(x) == 1L && is.finite(x)) x else NA_real_
}

# 3. 20 RUNS x 3 MODELS --------------------------------------------------------
cat("\nUCI Wine: n=", n, ", p=", ncol(X), ", k=", K_LAYERS[[1L]],
    ", 20 runs per model, one initialization per run.\n", sep = "")
for (model in names(models)) {
  family <- unname(models[[model]])
  assignments[[model]] <- list()
  for (r in seq_len(N_RUNS)) {
    seed <- as.integer(FIRST_SEED + r - 1L)
    row <- data.frame(Model = model, Run = r, Seed = seed,
                      Status = "Fit failed", BIC = NA_real_,
                      LogLik = NA_real_, ARI = NA_real_, MR = NA_real_,
                      Nu_Hat = NA_real_, Nonempty = NA_integer_,
                      Fit_Seconds = NA_real_, Iterations = NA_real_,
                      Error = "", stringsAsFactors = FALSE)
    start <- proc.time()[["elapsed"]]
    control <- list(seed = seed, nstart = 1L,
                    max_iter = MAX_ITER, min_iter = MIN_ITER, tol = TOL,
                    nu = NU_INITIAL, estimate_nu = ESTIMATE_NU,
                    mc_samples = MC_SAMPLES, psi_floor = 0.03,
                    verbose = FALSE)
    fit <- tryCatch(MixedDeepMM::mixeddeepmm(
      data = X, family = family, heads = "single",
      layers = layers, control = control), error = function(e) e)
    row$Fit_Seconds <- proc.time()[["elapsed"]] - start
    if (inherits(fit, "error")) {
      row$Error <- conditionMessage(fit)
    } else {
      row$BIC <- get_num(fit$bic)
      row$LogLik <- get_num(fit$loglik)
      row$Iterations <- get_num(fit$iterations)
      if (family != "gaussian") row$Nu_Hat <- get_num(fit$model$nu)
      valid <- is.finite(row$BIC) && is.finite(row$LogLik) &&
        (family == "gaussian" ||
         (is.finite(row$Nu_Hat) && row$Nu_Hat > 2))
      if (!valid) {
        row$Status <- "Invalid fit"
        row$Error <- "Nonfinite BIC/logLik or invalid fitted nu."
      } else {
        pred <- tryCatch({
          set.seed(as.integer(seed + 100000L))
          p <- stats::predict(fit, type = "layer", layer = 1L)
          if (!is.matrix(p) || nrow(p) != n ||
              ncol(p) != K_LAYERS[[1L]] || any(!is.finite(p)) ||
              any(abs(rowSums(p) - 1) > 1e-3)) {
            stop("Invalid first-layer posterior matrix.")
          }
          max.col(p, ties.method = "first")
        }, error = function(e) e)
        if (inherits(pred, "error")) {
          row$Status <- "Prediction failed"
          row$Error <- conditionMessage(pred)
        } else {
          row$ARI <- get_num(MixedDeepMM::mixeddeepmm_ari(truth, pred))
          row$MR <- get_num(MixedDeepMM::mixeddeepmm_mr(truth, pred))
          row$Nonempty <- length(unique(pred))
          if (!is.finite(row$ARI) || !is.finite(row$MR)) {
            row$Status <- "Invalid metrics"
            row$Error <- "ARI or MR is nonfinite."
          } else {
            row$Status <- "OK"
            assignments[[model]][[as.character(r)]] <- as.integer(pred)
          }
        }
      }
    }
    run_rows[[length(run_rows) + 1L]] <- row
    cat(sprintf("%-6s | Run %02d/%02d | %-17s | BIC %s | ARI %s | MR %s | Clusters %s\n",
                model, r, N_RUNS, row$Status,
                if (is.finite(row$BIC)) sprintf("%.2f", row$BIC) else "NA",
                if (is.finite(row$ARI)) sprintf("%.4f", row$ARI) else "NA",
                if (is.finite(row$MR)) sprintf("%.4f", row$MR) else "NA",
                if (is.finite(row$Nonempty)) as.character(row$Nonempty) else "NA"))
    if (nzchar(row$Error)) cat("  Reason: ", row$Error, "\n", sep = "")
    rm(fit)
  }
}
all_runs <- do.call(rbind, run_rows)

# 4. ONE SUMMARY CSV: 20-RUN MEAN / SE ----------------------------------------
model_summary <- function(model) {
  z <- all_runs[all_runs$Model == model, , drop = FALSE]
  good <- z[z$Status == "OK", , drop = FALSE]
  bic <- safe_stat(good$BIC)
  ari <- safe_stat(good$ARI)
  mr <- safe_stat(good$MR)
  nu <- safe_stat(good$Nu_Hat)
  secs <- safe_stat(good$Fit_Seconds)
  data.frame(Model = model, Runs_Requested = N_RUNS,
             Runs_OK = nrow(good), Runs_Failed = N_RUNS - nrow(good),
             Cluster_Collapse_Runs = sum(good$Nonempty < K_LAYERS[[1L]]),
             BIC_Mean = bic["mean"], BIC_SD = bic["sd"], BIC_SE = bic["se"],
             ARI_Mean = ari["mean"], ARI_SD = ari["sd"], ARI_SE = ari["se"],
             MR_Mean = mr["mean"], MR_SD = mr["sd"], MR_SE = mr["se"],
             Nu_Mean = nu["mean"], Nu_SD = nu["sd"], Nu_SE = nu["se"],
             FitSec_Mean = secs["mean"], FitSec_SE = secs["se"],
             stringsAsFactors = FALSE, row.names = NULL)
}
summary_table <- do.call(rbind, lapply(names(models), model_summary))
row.names(summary_table) <- NULL
csv_file <- file.path(OUT_DIR, "Wine_20Runs_Mean_SE.csv")
utils::write.csv(summary_table, csv_file, row.names = FALSE, na = "")
cat("\n========== FINAL 20-RUN MEAN +/- SE ==========\n")
print(summary_table, row.names = FALSE)

# 5. 20-RUN CO-ASSIGNMENT CONSENSUS, THEN PCA FIGURE --------------------------
# Labels 1/2/3 may permute between fits, so averaging numeric cluster IDs
# would be meaningless. Coassignment is invariant to label switching.
# The number of clusters in each consensus partition is the MODAL number of
# occupied clusters across successful runs (NOT automatically forced to 3).
consensus_partition <- function(parts, n) {
  R <- length(parts)
  if (R == 0L) return(rep(NA_integer_, n))
  C <- matrix(0, n, n)
  occupancies <- integer(R)
  for (j in seq_len(R)) {
    g <- parts[[j]]
    if (length(g) != n || anyNA(g)) stop("Bad partition length in consensus.")
    occupancies[[j]] <- length(unique(g))
    C <- C + outer(g, g, "==")
  }
  C <- C / R
  freq <- table(occupancies)
  k <- as.integer(names(freq)[which.max(as.numeric(freq))])
  if (k < 2L) return(rep(1L, n))
  distance <- stats::as.dist(pmax(0, 1 - C))
  D <- 1 - C
  D[D < 0] <- 0
  D[D > 1] <- 1
  D <- (D + t(D)) / 2
  diag(D) <- 0
  
  tree <- stats::hclust(
    stats::as.dist(D),
    method = "average"
  )
  as.integer(stats::cutree(tree, k = min(k, n)))
}
consensus <- lapply(assignments, consensus_partition, n = n)

# Every panel shows the SAME PCA coordinates and axis range; colours denote
# the true classes in panel A or the consensus group in model panels B-D.
colors <- c("#2C76B8", "#DB8F27", "#2B9C82")
xy_lim <- range(proj[, 1], finite = TRUE)
y_lim <- range(proj[, 2], finite = TRUE)
expand <- function(z) {
  delta <- diff(z)
  if (delta <= 0) delta <- 1
  z + c(-1, 1) * delta * 0.045
}
xy_lim <- expand(xy_lim)
y_lim <- expand(y_lim)

plot_panel <- function(groups, heading, subheading = "", true_panel = FALSE) {
  if (all(is.na(groups))) {
    graphics::plot.new()
    graphics::title(main = paste(heading, "(no successful fits)"))
    return(invisible(NULL))
  }
  g <- factor(groups)
  ids <- as.integer(g)
  graphics::plot(proj[, 1], proj[, 2], type = "n", xlim = xy_lim,
                 ylim = y_lim, xlab = "PC1 (standardized features)",
                 ylab = "PC2", main = heading)
  graphics::abline(h = 0, v = 0, col = "#E8E8E8", lty = 3)
  graphics::points(proj[, 1], proj[, 2], col = colors[ids], pch = 16, cex = 0.77)
  graphics::legend("topright", legend = paste(if (true_panel) "Class" else "Cluster",
                                          levels(g)),
                   pch = 16, col = colors[seq_along(levels(g))],
                   bty = "n", cex = 0.76)
  if (nzchar(subheading)) graphics::mtext(subheading, side = 3, line = 0.15, cex = 0.80)
}

plot_all <- function() {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4.6, 4.5, 3.4, 1.3),
                oma = c(0.7, 0.7, 1.5, 0.7), cex = 0.98)
  plot_panel(truth, "(A) True Wine classes", "n = 178; 13 continuous inputs", TRUE)
  for (j in seq_along(models)) {
    nm <- names(models)[[j]]
    zz <- summary_table[summary_table$Model == nm, , drop = FALSE]
    note <- sprintf("ARI %.3f (%.3f)",
                    zz$BIC_Mean, zz$BIC_SE, zz$ARI_Mean, zz$ARI_SE)
    plot_panel(consensus[[nm]], paste0("(", LETTERS[j + 1L], ") ", nm,
                                              " consensus"), note)
  }
  graphics::mtext("UCI Wine",
                  side = 3, outer = TRUE, line = 0.0, cex = 1.0, font = 2)
}

png_file <- file.path(OUT_DIR, "Wine_PCA_20Runs_Consensus.png")
pdf_file <- file.path(OUT_DIR, "Wine_PCA_20Runs_Consensus.pdf")
grDevices::png(png_file, width = 2400L, height = 1800L,
               res = 200L, pointsize = 15L)
tryCatch(plot_all(), finally = grDevices::dev.off())
grDevices::pdf(pdf_file, width = 12, height = 9, pointsize = 15,
               useDingbats = FALSE)
tryCatch(plot_all(), finally = grDevices::dev.off())
if (interactive() && isTRUE(SHOW_PLOTS_IN_RSTUDIO)) plot_all()
cat("\nFinished. Final output folder: ", OUT_DIR, "\n", sep = "")
cat("Files: Wine_20Runs_Mean_SE.csv, Wine_PCA_20Runs_Consensus.png/pdf\n")
cat("Note: SE describes across-seed variability on ONE fixed dataset, not a data-sampling confidence interval.\n")
