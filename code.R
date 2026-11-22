# Cancer Classification with High-Dimensional and Imbalanced Gene Expression Data
# Experiment 1: six sampling conditions after FR with K = 500.
# Experiment 2: Fisher, ReliefF and FR; no sampling and the retained sampler.
# Experiment 3: FR at K = 100, 200, ..., 1000 with the retained sampler.
EXPERIMENT <- 1L
DATA_DIR <- "data"
OUTPUT_DIR <- "results"
DATASET_IDS <- c("D1", "D2", "D3")
CPU_BUDGET <- 100L
PARALLEL_FOLDS <- 5L

if (.Platform$OS.type == "windows") suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
suppressWarnings(try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE))
options(stringsAsFactors = FALSE, timeout = 7200, repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

DATASET_SPECS <- data.frame(
  id = c("D1", "D2", "D3"),
  dataset = c("D1_SCANB_HER2", "D2_GSE106817_LUNG", "D3_GSE211692_SARCOMA"),
  expected_n = c(3151L, 4046L, 12284L), expected_positive = c(420L, 115L, 187L),
  expected_p = c(30865L, 2550L, 2550L),
  prefix = c("symbol:", "mirbase:", "mirna:"),
  sampler = c("Borderline_SMOTE1", "ADASYN", "SMOTE"),
  ir_min = c(5, 30, 50), ir_max = c(10, 50, 100)
)
ALL_SAMPLERS <- c("None", "RUS", "ROS", "SMOTE", "ADASYN", "Borderline_SMOTE1")
DATASETS <- DIR <- DOWNLOAD_FILES <- CFG <- NULL
EXPERIMENT_ID <- DATASET_INDEX <- WORKERS_PER_FOLD <- NULL
DATASET <- FEATURE_PREFIX <- EXPECTED_FILE <- EXPECTED_N <- EXPECTED_POSITIVE <- EXPECTED_P <- NULL
K_VALUES <- 500L



need_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) {
    user_lib <- strsplit(path.expand(Sys.getenv("R_LIBS_USER")), .Platform$path.sep, fixed=TRUE)[[1]][1]
    if(!is.na(user_lib) && nzchar(user_lib)) {
      dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
      .libPaths(unique(c(user_lib, .libPaths())))
    }
    install.packages(missing, dependencies=NA)
  }
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop("Install missing packages: ", paste(missing,collapse=", "))
}



initialize_dependencies <- function(analysis = TRUE) {
  pkgs <- if (analysis) c("data.table", "matrixStats", "RANN", "ranger", "gbm", "xgboost", "pROC") else "data.table"
  need_packages(pkgs)
  suppressPackageStartupMessages(library(data.table))
  DATASETS <<- as.data.table(DATASET_SPECS[, c("dataset", "expected_n", "expected_positive", "expected_p", "ir_min", "ir_max")])
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores)) cores <- 1L
  CPU_BUDGET <<- max(1L, min(as.integer(CPU_BUDGET), cores))
  PARALLEL_FOLDS <<- max(1L, min(5L, as.integer(PARALLEL_FOLDS), CPU_BUDGET))
  WORKERS_PER_FOLD <<- max(1L, as.integer(CPU_BUDGET %/% PARALLEL_FOLDS))
}

configure_experiment <- function(experiment, dataset_index) {
  spec <- DATASET_SPECS[dataset_index, ]
  EXPERIMENT_ID <<- as.integer(experiment)
  DATASET_INDEX <<- as.integer(dataset_index)
  DATASET <<- spec$dataset
  FEATURE_PREFIX <<- spec$prefix
  EXPECTED_FILE <<- paste0(DATASET, ".csv")
  EXPECTED_N <<- spec$expected_n
  EXPECTED_POSITIVE <<- spec$expected_positive
  EXPECTED_P <<- spec$expected_p
  K_VALUES <<- if (experiment == 3L) seq.int(100L, 1000L, by = 100L) else 500L
  selectors <- switch(as.character(experiment),
    "1" = "FR_intersection_top500",
    "2" = c("Fisher", "ReliefF", "FR_intersection_top500"),
    "3" = paste0("FR_K", K_VALUES), stop("Experiment must be 1, 2 or 3."))
  methods <- switch(as.character(experiment), "1" = ALL_SAMPLERS,
                    "2" = c("None", spec$sampler), "3" = spec$sampler)
  CFG <<- list(seed = 20261121L, train_fraction = 0.7, folds = 5L, top_k = max(K_VALUES),
    relief_reference_per_class = 250L, relief_batch_size = 16L, relief_neighbors = 10L,
    threads = 1L, max_missing_fraction = 0.2, sampling_neighbors = 5L,
    borderline_neighbors = 10L, methods = methods, selectors = selectors,
    models = c("RF", "GBM", "XGBoost"), classification_threshold = 0.5)
}

# Training preprocessing, independent rankings and resampling.


safe_seed <- function(...) as.integer((sum(as.double(unlist(list(...))))%%2147483000) + 1)


column_blocks <- function(p, size = 256L) split(seq_len(p), ceiling(seq_len(p)/size))


fit_preprocessor <- function(x, max_missing = CFG$max_missing_fraction) {
    missing <- colMeans(!is.finite(x))
    keep <- which(missing <= max_missing)
    if (length(keep) < CFG$top_k) 
        stop("Fewer than K eligible features")
    z <- x[, keep, drop = FALSE]
    z[!is.finite(z)] <- NA_real_
    med <- matrixStats::colMedians(z, na.rm = TRUE)
    for (j in which(colSums(is.na(z)) > 0)) z[is.na(z[, j]), j] <- med[j]
    mn <- matrixStats::colMins(z)
    mx <- matrixStats::colMaxs(z)
    good <- is.finite(mx - mn) & mx > mn
    z <- z[, good, drop = FALSE]
    keep <- keep[good]
    med <- med[good]
    mn <- mn[good]
    mx <- mx[good]
    if (ncol(z) < CFG$top_k) 
        stop("Fewer than K nonconstant training genes")
    for (cols in column_blocks(ncol(z))) z[, cols] <- sweep(sweep(z[, cols, drop = FALSE], 2, mn[cols], "-"), 2, mx[cols] - 
        mn[cols], "/")
    list(x = z, features = colnames(x)[keep], median = med, min = mn, range = mx - mn)
}


apply_preprocessor <- function(x, pp, features = pp$features) {
    if (!length(features)) 
        return(matrix(numeric(), nrow(x), 0L, dimnames = list(rownames(x), character())))
    j <- match(features, pp$features)
    if (anyNA(j)) 
        stop("Unknown preprocessing feature")
    z <- x[, features, drop = FALSE]
    z[!is.finite(z)] <- NA_real_
    for (k in which(colSums(is.na(z)) > 0)) z[is.na(z[, k]), k] <- pp$median[j[k]]
    sweep(sweep(z, 2, pp$min[j], "-"), 2, pp$range[j], "/")
}


fisher_scores <- function(x, y) {
    n0 <- sum(y == 0L)
    n1 <- sum(y == 1L)
    m0 <- colMeans(x[y == 0, , drop = FALSE])
    m1 <- colMeans(x[y == 1, , drop = FALSE])
    m <- (n0 * m0 + n1 * m1)/(n0 + n1)
    numerator <- n0 * (m0 - m)^2 + n1 * (m1 - m)^2
    denominator <- (n0 - 1) * matrixStats::colVars(x[y == 0, , drop = FALSE]) + (n1 - 1) * matrixStats::colVars(x[y == 1, 
        , drop = FALSE])
    out <- numerator/denominator
    out[denominator == 0 & numerator > 0] <- Inf
    out[denominator == 0 & numerator == 0] <- 0
    names(out) <- colnames(x)
    out
}


relieff_scores <- function(x, y, seed, k = CFG$relief_neighbors, per_class = CFG$relief_reference_per_class, batch_size = CFG$relief_batch_size) {
    stopifnot(all(y %in% 0:1), all(is.finite(x)))
    idx <- list(which(y == 0L), which(y == 1L))
    if (min(lengths(idx)) < 2) 
        stop("ReliefF requires >=2 rows per class")
    set.seed(seed)
    refs <- unlist(lapply(idx, function(z) z[sample.int(length(z), min(per_class, length(z)))]), use.names = FALSE)
    nref <- table(factor(y[refs], levels = 0:1))
    prior <- lengths(idx)/length(y)
    weights <- prior[y[refs] + 1L]/as.numeric(nref[y[refs] + 1L])
    norm2 <- numeric(nrow(x))
    for (cols in column_blocks(ncol(x))) norm2 <- norm2 + rowSums(x[, cols, drop = FALSE]^2)
    scores <- numeric(ncol(x))
    k <- min(k, min(lengths(idx)) - 1L)
    for (batch in split(seq_along(refs), ceiling(seq_along(refs)/batch_size))) {
        anchors <- refs[batch]
        d <- outer(norm2, norm2[anchors], "+") - 2 * tcrossprod(x, x[anchors, , drop = FALSE])
        d <- pmax(d, 0)
        for (b in seq_along(anchors)) {
            r <- anchors[b]
            hit <- idx[[y[r] + 1L]]
            hit <- hit[hit != r]
            miss <- idx[[2L - y[r]]]
            h <- hit[order(d[hit, b], hit)[seq_len(k)]]
            v <- miss[order(d[miss, b], miss)[seq_len(k)]]
            delta <- colMeans(abs(sweep(x[v, , drop = FALSE], 2, x[r, ], "-"))) - colMeans(abs(sweep(x[h, , drop = FALSE], 
                2, x[r, ], "-")))
            scores <- scores + weights[batch[b]] * delta
        }
    }
    names(scores) <- colnames(x)
    list(scores = scores, references = refs, reference_weights = weights, k = k)
}


rank_genes <- function(x, y, seed) {
    f <- fisher_scores(x, y)
    r <- relieff_scores(x, y, seed)
    fn <- names(f)[order(-f, names(f))]
    rn <- names(r$scores)[order(-r$scores, names(r$scores))]
    if (ncol(x) < max(K_VALUES)) 
        stop("Fewer eligible training features than max(K_VALUES).")
    list(fisher = fn, relieff = rn)
}


nearest <- function(ref, query, k, self_index = NULL) {
    ask <- min(nrow(ref), k + as.integer(!is.null(self_index)))
    a <- RANN::nn2(ref, query, k = ask, eps = 0)$nn.idx
    if (is.null(dim(a))) 
        a <- matrix(a, nrow = nrow(query))
    out <- matrix(NA_integer_, nrow(query), k)
    for (i in seq_len(nrow(query))) {
        v <- a[i, ]
        if (!is.null(self_index)) 
            v <- v[v != self_index[i]]
        if (length(v) < k) 
            stop("Insufficient nonself neighbors")
        out[i, ] <- v[seq_len(k)]
    }
    out
}


largest_remainder <- function(n, weights) {
    v <- n * weights/sum(weights)
    out <- floor(v)
    left <- n - sum(out)
    if (left > 0) 
        out[order(-(v - out), seq_along(v))[seq_len(left)]] <- out[order(-(v - out), seq_along(v))[seq_len(left)]] + 1L
    as.integer(out)
}


resample_training <- function(x, y, method, seed) {
    set.seed(seed)
    i0 <- which(y == 0L)
    i1 <- which(y == 1L)
    n0 <- length(i0)
    n1 <- length(i1)
    if (n1 >= n0) 
        stop("Expected positive minority class")
    extra <- 0L
    danger_n <- NA_integer_
    k <- min(CFG$sampling_neighbors, n1 - 1L)
    m <- min(CFG$borderline_neighbors, nrow(x) - 1L)
    if (method == "None") 
        ix <- seq_along(y)
    else if (method == "RUS") 
        ix <- c(i0[sample.int(n0, n1)], i1)
    else if (method == "ROS") 
        ix <- c(seq_along(y), i1[sample.int(n1, n0 - n1, replace = TRUE)])
    else {
        if (k < 1) 
            stop("Synthetic sampling requires >=2 minority rows")
        xm <- x[i1, , drop = FALSE]
        kn <- nearest(xm, xm, k, self_index = seq_len(n1))
        target <- n0 - n1
        if (method == "SMOTE") 
            bases <- sample.int(n1, target, replace = TRUE)
        else if (method == "ADASYN") {
            alln <- nearest(x, xm, min(CFG$sampling_neighbors, nrow(x) - 1L), self_index = i1)
            ratio <- rowMeans(matrix(y[alln] == 0L, nrow = n1))
            if (sum(ratio) == 0) 
                stop("NOT_APPLICABLE: ADASYN has zero majority-neighbor density; no fallback sampler used")
            counts <- largest_remainder(target, ratio)
            bases <- rep(seq_len(n1), counts)
        }
        else if (method == "Borderline_SMOTE1") {
            alln <- nearest(x, xm, m, self_index = i1)
            maj <- rowSums(matrix(y[alln] == 0L, nrow = n1))
            danger <- which(maj >= m/2 & maj < m)
            danger_n <- length(danger)
            if (!danger_n) 
                stop("NOT_APPLICABLE: no Borderline-SMOTE1 danger samples; no fallback sampler used")
            bases <- danger[sample.int(danger_n, target, replace = TRUE)]
        }
        else stop("Unknown sampler: ", method)
        partner <- kn[cbind(bases, sample.int(k, length(bases), replace = TRUE))]
        lam <- runif(length(bases))
        syn <- xm[bases, , drop = FALSE] + lam * (xm[partner, , drop = FALSE] - xm[bases, , drop = FALSE])
        x <- rbind(x, syn)
        y <- c(y, rep(1L, nrow(syn)))
        ix <- seq_along(y)
        extra <- nrow(syn)
    }
    x <- x[ix, , drop = FALSE]
    y <- y[ix]
    list(x = x, y = y, audit = data.table(method = method, real_negative = n0, real_positive = n1, sampled_negative = sum(y == 
        0L), sampled_positive = sum(y == 1L), synthetic_rows = extra, minority_k = if (method %in% c("SMOTE", "ADASYN", "Borderline_SMOTE1")) k else NA_integer_, 
        danger_samples = danger_n, seed = seed, status = "OK"))
}


make_grids <- function(p) {
    rf <- CJ(mtry_rule = c("sqrt", "one_fifth", "one_third"), min_node_size = c(1L, 5L))
    rf[, `:=`(num_trees, 300L)]
    gb <- CJ(n_trees = c(200L, 500L), interaction_depth = c(1L, 3L), shrinkage = c(0.05, 0.1))
    gb[, `:=`(n_minobsinnode = 10L, bag_fraction = 0.8)]
    xg <- CJ(nrounds = c(200L, 500L), max_depth = c(2L, 4L, 6L), eta = c(0.05, 0.1))
    xg[, `:=`(min_child_weight = 1, subsample = 0.8, colsample_bytree = 0.8, lambda = 1)]
    list(RF = rf, GBM = gb, XGBoost = xg)
}


fit_predict <- function(model, x, y, newx, par, seed) {
    set.seed(seed)
    if (model == "RF") {
        fit <- ranger::ranger(x = as.data.frame(x), y = factor(y, levels = 0:1), probability = TRUE, num.trees = par$num_trees, 
            mtry = rf_mtry(par$mtry_rule, ncol(x)), min.node.size = par$min_node_size, seed = seed, num.threads = CFG$threads, 
            write.forest = TRUE)
        p <- predict(fit, data = as.data.frame(newx), num.threads = CFG$threads)$predictions[, "1"]
    }
    else if (model == "GBM") {
        fit <- gbm::gbm.fit(x = x, y = y, distribution = "bernoulli", n.trees = par$n_trees, interaction.depth = par$interaction_depth, 
            shrinkage = par$shrinkage, n.minobsinnode = par$n_minobsinnode, bag.fraction = par$bag_fraction, nTrain = nrow(x), 
            keep.data = FALSE, verbose = FALSE)
        p <- predict(fit, newdata = newx, n.trees = par$n_trees, type = "response")
    }
    else if (model == "XGBoost") {
        params <- list(objective = "binary:logistic", eval_metric = "logloss", max_depth = par$max_depth, eta = par$eta, 
            min_child_weight = par$min_child_weight, subsample = par$subsample, colsample_bytree = par$colsample_bytree, 
            lambda = par$lambda, nthread = CFG$threads, seed = seed, tree_method = "hist")
        fit <- xgboost::xgb.train(params = params, data = xgboost::xgb.DMatrix(x, label = y), nrounds = par$nrounds, verbose = 0)
        p <- predict(fit, xgboost::xgb.DMatrix(newx))
    }
    else stop("Unknown model")
    as.numeric(p)
}


rf_mtry <- function(rule, p) max(1L, min(p, as.integer(switch(rule, sqrt = sqrt(p), one_fifth = p/5, one_third = p/3))))


split_dataset <- function(d, dataset_index) {
    stopifnot(!anyDuplicated(d$sample_id), !anyDuplicated(d$patient_id))
    set.seed(safe_seed(CFG$seed, 100 * dataset_index))
    train <- integer()
    for (v in 0:1) {
        ix <- which(d$y == v)
        train <- c(train, ix[sample.int(length(ix), floor(CFG$train_fraction * length(ix)))])
    }
    train <- sort(train)
    test <- setdiff(seq_along(d$y), train)
    fold <- rep(NA_integer_, length(d$y))
    for (v in 0:1) {
        ix <- train[d$y[train] == v]
        ix <- ix[sample.int(length(ix))]
        fold[ix] <- rep(seq_len(CFG$folds), length.out = length(ix))
    }
    if (min(table(d$y[train])) < CFG$folds) 
        stop("Too few positives for stratified CV")
    data.table(sample_id = d$sample_id, patient_id = d$patient_id, label = d$y, partition = ifelse(seq_along(d$y) %in% train, 
        "train", "test"), validation_fold = fold)
}


parameter_text <- function(model, par, p) {
    if (model == "RF") 
        return(paste0("num.trees=", par$num_trees, "; mtry=", rf_mtry(par$mtry_rule, p), "; mtry_rule=", par$mtry_rule, "; min.node.size=", 
            par$min_node_size, "; probability=TRUE"))
    if (model == "GBM") 
        return(paste0("n.trees=", par$n_trees, "; interaction.depth=", par$interaction_depth, "; shrinkage=", par$shrinkage, 
            "; n.minobsinnode=", par$n_minobsinnode, "; bag.fraction=", par$bag_fraction, "; distribution=bernoulli"))
    paste0("nrounds=", par$nrounds, "; max_depth=", par$max_depth, "; eta=", par$eta, "; min_child_weight=", par$min_child_weight, 
        "; subsample=", par$subsample, "; colsample_bytree=", par$colsample_bytree, "; lambda=", par$lambda, "; alpha=0; gamma=0; max_bin=256; objective=binary:logistic; tree_method=hist")
}


read_input <- function(path) {
    z <- fread(path, na.strings = c("NA", "NaN", ""), check.names = FALSE)
    required <- c("sample_id", "patient_id", "label")
    if (!all(required %in% names(z))) 
        stop("Use a processed CSV generated by code.R --download-only.")
    if (anyNA(z$sample_id) || anyNA(z$patient_id) || anyDuplicated(z$sample_id) || anyDuplicated(z$patient_id)) 
        stop("Missing or repeated sample/patient identifiers.")
    if (anyNA(z$label) || !all(z$label %in% 0:1)) 
        stop("Labels must be 0 or 1.")
    if (sum(z$label == 1L) >= sum(z$label == 0L) || sum(z$label == 1L) < 30L) 
        stop("Invalid minority class size.")
    features <- setdiff(names(z), required)
    if (!length(features) || anyDuplicated(features) || !all(startsWith(features, FEATURE_PREFIX))) 
        stop("Wrong dataset CSV: expected ", DATASET, " with feature prefix ", FEATURE_PREFIX)
    if (nrow(z) != EXPECTED_N || sum(z$label == 1L) != EXPECTED_POSITIVE || length(features) != EXPECTED_P) 
        stop("Unexpected dimensions/class counts for ", DATASET, ". Use the unmodified prepared CSV.")
    if (!all(vapply(z[, ..features], is.numeric, logical(1)))) 
        stop("Non-numeric expression columns.")
    d <- list(y = as.integer(z$label), sample_id = as.character(z$sample_id), patient_id = as.character(z$patient_id), x = as.matrix(z[, 
        ..features]))
    storage.mode(d$x) <- "double"
    rm(z)
    invisible(gc())
    d
}


init_worker <- function(paths) {
    .libPaths(paths)
    if (.Platform$OS.type == "windows") 
        suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
    suppressWarnings(try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE))
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    suppressPackageStartupMessages(library(data.table))
    data.table::setDTthreads(1L)
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) 
        RhpcBLASctl::blas_set_num_threads(1L)
    invisible(TRUE)
}


condition_key <- function(selector, sampler, model = NULL) paste(c(selector, sampler, model), collapse = "|")


make_jobs <- function(grids, final = FALSE, chosen = NULL) {
    jobs <- list()
    for (selector in CFG$selectors) for (sampler in CFG$methods) for (mi in seq_along(CFG$models)) {
        model <- CFG$models[mi]
        key <- condition_key(selector, sampler, model)
        gs <- if (final) {
            if (is.null(chosen[[key]])) 
                integer()
            else chosen[[key]]$grid_id
        }
        else seq_len(nrow(grids[[model]]))
        for (g in gs) jobs[[length(jobs) + 1L]] <- list(selector = selector, sampler = sampler, model = model, mi = mi, g = g)
    }
    jobs
}


fit_job <- function(job) {
    a <- WORK_PAYLOAD$sampled[[condition_key(job$selector, job$sampler)]]
    if (inherits(a, "error")) 
        return(c(job, list(error = conditionMessage(a))))
    par <- as.list(WORK_PAYLOAD$grids[[job$model]][job$g])
    offset <- if (WORK_PAYLOAD$final) 
        9990000L
    else WORK_PAYLOAD$fold * 10000L
    seed <- safe_seed(CFG$seed, DATASET_INDEX * 1e+06, offset, job$mi * 100, job$g)
    newx <- WORK_PAYLOAD$newx[[job$selector]]
    p <- tryCatch(fit_predict(job$model, a$x, a$y, newx, par, seed), error = function(e) e)
    if (inherits(p, "error")) 
        return(c(job, list(error = conditionMessage(p))))
    if (length(p) != nrow(newx) || any(!is.finite(p))) 
        return(c(job, list(error = "Invalid probabilities")))
    c(job, list(p = p, error = ""))
}


run_job_pool <- function(payload, jobs, n) {
    if (!length(jobs)) 
        return(list())
    WORK_PAYLOAD <- payload
    if (n <= 1L) {
        f <- fit_job
        environment(f) <- environment()
        return(lapply(jobs, f))
    }
    cl <- parallel::makePSOCKcluster(min(n, length(jobs)), outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, init_worker, .libPaths())
    parallel::clusterExport(cl, c("fit_predict", "rf_mtry", "safe_seed", "condition_key", "CFG", "DATASET_INDEX", "fit_job", 
        "WORK_PAYLOAD"), envir = environment())
    parallel::parLapplyLB(cl, jobs, fit_job)
}


make_sampled_sets <- function(x, y, f, final = FALSE) {
    out <- list()
    for (selector in CFG$selectors) for (sampler in CFG$methods) {
        si <- match(sampler, c("None", "RUS", "ROS", "SMOTE", "ADASYN", "Borderline_SMOTE1"))
        out[[condition_key(selector, sampler)]] <- tryCatch({
            if (ncol(x[[selector]]) == 0L) 
                stop("EMPTY_FR_INTERSECTION: no substitute features added")
            resample_training(x[[selector]], y, sampler, safe_seed(CFG$seed, DATASET_INDEX * 1e+05, if (final) 
                999000L
            else f * 1000L, si * 10))
        }, error = function(e) e)
    }
    out
}


run_fold_parallel <- function(f) {
    message("CV fold ", f, ": training-only preprocessing and both rankings.")
    ia <- itrain[fold != f]
    iv <- itrain[fold == f]
    pp <- fit_preprocessor(d$x[ia, , drop = FALSE])
    sets <- feature_sets(pp$x, d$y[ia], safe_seed(CFG$seed, DATASET_INDEX * 10000, f))
    xa <- lapply(sets, function(sel) pp$x[, sel, drop = FALSE])
    xv <- lapply(sets, function(sel) apply_preprocessor(d$x[iv, , drop = FALSE], pp, sel))
    rm(pp)
    invisible(gc())
    sampled <- make_sampled_sets(xa, d$y[ia], f)
    message("CV fold ", f, ": selected features ", paste(names(sets), lengths(sets), collapse = "; "), "; model workers=", 
        WORKERS_PER_FOLD)
    payload <- list(sampled = sampled, newx = xv, grids = grids, fold = f, final = FALSE)
    ans <- list(fold = f, selected_counts = lengths(sets), predictions = run_job_pool(payload, make_jobs(grids), WORKERS_PER_FOLD))
    message("CV fold ", f, ": complete.")
    ans
}



metric_values <- function(y, p, threshold = 0.5) {
  stopifnot(length(y) == length(p), all(is.finite(p)), all(y %in% 0:1))
  pred <- as.integer(p >= threshold)
  tp <- sum(pred == 1L & y == 1L); fp <- sum(pred == 1L & y == 0L)
  fn <- sum(pred == 0L & y == 1L); tn <- sum(pred == 0L & y == 0L)
  auc <- if (length(unique(y)) == 2L) as.numeric(pROC::auc(pROC::roc(
    y, p, levels = c(0, 1), direction = "<", quiet = TRUE))) else NA_real_
  c(Accuracy = (tp + tn) / length(y), Recall = if (tp + fn) tp / (tp + fn) else NA_real_,
    F1 = if (2 * tp + fp + fn) 2 * tp / (2 * tp + fp + fn) else 0, AUC = auc)
}

choose_parameters <- function(prob, y, fold, grid) {
  scores <- lapply(seq_len(nrow(grid)), function(g) {
    if (any(!is.finite(prob[, g]))) return(NULL)
    f1 <- vapply(seq_len(CFG$folds), function(f) {
      ix <- which(fold == f)
      unname(metric_values(y[ix], prob[ix, g], CFG$classification_threshold)["F1"])
    }, numeric(1))
    if (any(!is.finite(f1))) return(NULL)
    data.table(grid_id = g, F1 = mean(f1))
  })
  scores <- rbindlist(scores)
  if (!nrow(scores)) return(NULL)
  setorderv(scores, c("F1", "grid_id"), c(-1L, 1L))
  scores[1L]
}

feature_sets <- function(x, y, seed) {
  ranks <- rank_genes(x, y, seed)
  if (EXPERIMENT_ID == 2L) return(list(Fisher = ranks$fisher[seq_len(500L)],
    ReliefF = ranks$relieff[seq_len(500L)],
    FR_intersection_top500 = intersect(ranks$fisher[seq_len(500L)], ranks$relieff[seq_len(500L)])))
  setNames(lapply(K_VALUES, function(k) intersect(ranks$fisher[seq_len(k)],
                                                ranks$relieff[seq_len(k)])), CFG$selectors)
}

run_experiment <- function(d) {
  splits <- split_dataset(d, DATASET_INDEX)
  itrain <- which(splits$partition == "train"); itest <- which(splits$partition == "test")
  ytrain <- d$y[itrain]; fold <- splits$validation_fold[itrain]
  grids <- make_grids(CFG$top_k)
  message(DATASET, ": training n=", length(itrain), " (positive=", sum(ytrain == 1L),
          "), test n=", length(itest), " (positive=", sum(d$y[itest] == 1L), ").")
  if (PARALLEL_FOLDS == 1L) {
    ffun <- run_fold_parallel; environment(ffun) <- environment()
    fold_results <- lapply(seq_len(CFG$folds), ffun)
  } else {
    cl <- parallel::makePSOCKcluster(PARALLEL_FOLDS, outfile = "")
    tryCatch({
      parallel::clusterCall(cl, init_worker, .libPaths())
      parallel::clusterExport(cl, WORKER_EXPORTS, envir = environment(run_experiment))
      parallel::clusterExport(cl, c("d", "itrain", "fold", "grids"), envir = environment())
      fold_results <- parallel::parLapplyLB(cl, seq_len(CFG$folds), run_fold_parallel)
    }, finally = parallel::stopCluster(cl))
  }
  pred <- list()
  for (selector in CFG$selectors) for (sampler in CFG$methods) for (model in CFG$models)
    pred[[condition_key(selector, sampler, model)]] <- matrix(NA_real_, length(itrain), nrow(grids[[model]]))
  for (z in fold_results) for (p in z$predictions) if (!nzchar(p$error))
    pred[[condition_key(p$selector, p$sampler, p$model)]][fold == z$fold, p$g] <- p$p
  chosen <- lapply(names(pred), function(key) {
    model <- strsplit(key, "|", fixed = TRUE)[[1L]][3L]
    choose_parameters(pred[[key]], ytrain, fold, grids[[model]])
  }); names(chosen) <- names(pred)
  fold_results <- fold_results[order(vapply(fold_results, function(z) z$fold, integer(1)))]
  cv_counts <- setNames(lapply(CFG$selectors, function(sel) vapply(fold_results,
    function(z) as.integer(z$selected_counts[[sel]]), integer(1))), CFG$selectors)
  rm(pred, fold_results); invisible(gc())
  message("CV complete. Refit on all training observations and evaluate the held-out test set.")
  pp <- fit_preprocessor(d$x[itrain, , drop = FALSE])
  sets <- feature_sets(pp$x, ytrain, safe_seed(CFG$seed, DATASET_INDEX * 10000, 999L))
  xtrain <- lapply(sets, function(sel) pp$x[, sel, drop = FALSE])
  xtest <- lapply(sets, function(sel) apply_preprocessor(d$x[itest, , drop = FALSE], pp, sel))
  rm(pp); invisible(gc())
  sampled <- make_sampled_sets(xtrain, ytrain, 999L, final = TRUE)
  payload <- list(sampled = sampled, newx = xtest, grids = grids, fold = 999L, final = TRUE)
  jobs <- make_jobs(grids, final = TRUE, chosen = chosen)
  fits <- run_job_pool(payload, jobs, min(length(jobs), CPU_BUDGET))
  assign("FR_LAST_TEST_PREDICTIONS", list(fits = fits, truth = d$y[itest],
    sample_id = d$sample_id[itest]), envir = .GlobalEnv)
  rows <- list()
  for (selector in CFG$selectors) for (sampler in CFG$methods) for (model in CFG$models) {
    key <- condition_key(selector, sampler, model); b <- chosen[[key]]
    a <- sampled[[condition_key(selector, sampler)]]
    k <- if (EXPERIMENT_ID == 3L) K_VALUES[match(selector, CFG$selectors)] else 500L
    row <- data.table(Experiment = EXPERIMENT_ID, Dataset = DATASET,
      Feature_selection = if (EXPERIMENT_ID == 3L) "FR_intersection" else selector,
      Candidate_K = k, Classifier = model, Sampling = sampler,
      Accuracy = NA_real_, Recall = NA_real_, F1 = NA_real_, AUC = NA_real_, Threshold = 0.5,
      Selected_features = length(sets[[selector]]), Training_N = length(itrain),
      Training_positive = sum(ytrain == 1L), Test_N = length(itest), Test_positive = sum(d$y[itest] == 1L),
      Fitted_negative = NA_integer_, Fitted_positive = NA_integer_, Final_hyperparameters = "", Status = "OK")
    for (f in seq_len(CFG$folds)) row[[paste0("CV_fold", f, "_selected_features")]] <- cv_counts[[selector]][f]
    if (!length(sets[[selector]])) row$Status <- "EMPTY_FR_INTERSECTION"
    else if (is.null(b)) row$Status <- "CV_FAILED: no candidate completed all five folds"
    else {
      row$Final_hyperparameters <- parameter_text(model, as.list(grids[[model]][b$grid_id]), ncol(xtrain[[selector]]))
      fit <- Filter(function(z) identical(z$selector, selector) && identical(z$sampler, sampler) &&
                      identical(z$model, model), fits)[[1L]]
      if (nzchar(fit$error)) row$Status <- fit$error else {
        row$Fitted_negative <- sum(a$y == 0L); row$Fitted_positive <- sum(a$y == 1L)
        values <- metric_values(d$y[itest], fit$p, 0.5)
        for (m in names(values)) row[[m]] <- unname(values[m])
      }
    }
    rows[[length(rows) + 1L]] <- row
  }
  ans <- rbindlist(rows)
  stopifnot(nrow(ans) == length(CFG$selectors) * length(CFG$methods) * length(CFG$models))
  if (EXPERIMENT_ID == 3L) ans <- ans[order(Candidate_K, match(Classifier, CFG$models))]
  else ans <- ans[order(match(Classifier, CFG$models), match(Feature_selection, CFG$selectors),
                        match(Sampling, CFG$methods))]
  ans
}

write_summary <- function(ans, path) {
  assign("FR_LAST_RESULTS", ans, envir = .GlobalEnv)
  tryCatch({ fwrite(ans, path, bom = TRUE, na = "NA"); path }, error = function(e) {
    alternate <- file.path(dirname(path), paste0(tools::file_path_sans_ext(basename(path)),
      "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
    message("Output is locked or unavailable; saving to ", alternate)
    tryCatch({ fwrite(ans, alternate, bom = TRUE, na = "NA"); alternate },
      error = function(e2) stop("Saving failed. Results remain in FR_LAST_RESULTS. ", conditionMessage(e2)))
  })
}

run_study <- function(experiment, data_dir = DATA_DIR, output_dir = OUTPUT_DIR, datasets = DATASET_IDS) {
  initialize_dependencies()
  ids <- match(datasets, DATASET_SPECS$id)
  if (anyNA(ids) || anyDuplicated(ids)) stop("Use distinct dataset IDs: D1, D2, D3.")
  files <- file.path(data_dir, paste0(DATASET_SPECS$dataset[ids], ".csv"))
  if (any(!file.exists(files))) stop("Missing processed CSV: ", paste(files[!file.exists(files)], collapse = ", "),
    "\nRun code.R --download-only --data-dir <folder> first, or set --data-dir to your existing CSV folder.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(output_dir, paste0("Experiment", experiment, "_test_results.csv"))
  results <- list()
  for (i in ids) {
    configure_experiment(experiment, i)
    d <- read_input(file.path(data_dir, EXPECTED_FILE))
    results[[length(results) + 1L]] <- run_experiment(d)
    saved <- write_summary(rbindlist(results), path)
    message("Saved ", nrow(rbindlist(results)), " conditions: ", saved)
    rm(d); invisible(gc())
  }
  ans <- rbindlist(results)
  if (any(ans$Status != "OK")) warning("Some conditions failed; inspect Status.")
  invisible(ans)
}

experiment1 <- function(data_dir = DATA_DIR, output_dir = OUTPUT_DIR, datasets = DATASET_IDS)
  run_study(1L, data_dir, output_dir, datasets)
experiment2 <- function(data_dir = DATA_DIR, output_dir = OUTPUT_DIR, datasets = DATASET_IDS)
  run_study(2L, data_dir, output_dir, datasets)
experiment3 <- function(data_dir = DATA_DIR, output_dir = OUTPUT_DIR, datasets = DATASET_IDS)
  run_study(3L, data_dir, output_dir, datasets)

WORKER_EXPORTS <- c("EXPERIMENT_ID", "K_VALUES", "CFG", "DATASET_INDEX", "WORKERS_PER_FOLD", "init_worker",
  "safe_seed", "column_blocks", "fit_preprocessor", "apply_preprocessor", "fisher_scores", "relieff_scores",
  "rank_genes", "nearest", "largest_remainder", "resample_training", "rf_mtry", "fit_predict", "condition_key",
  "feature_sets", "make_jobs", "fit_job", "run_job_pool", "make_sampled_sets", "run_fold_parallel")

# GEO acquisition. The source matrices are processed once before experiments.


download_files_manifest <- function() {
DOWNLOAD_FILES <- data.table(name=c('GSE96058_family.soft.gz','GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz','GSE106817_series_matrix.txt.gz','GSE211692_series_matrix.txt.gz','GSE211692_processed_data_matrix.txt.gz'),
bytes=c(442471,591676211,22986259,475860,230724363),
md5=c('df5df20772738e86e8d3f41f61185b4a','4362d1d9b5e80a6e1409cd7b9dc4a0ab','0bcec83204211163a21a07c1c9dc885e','96e929a40291294241e0615698a1abb4','21469a5f50690183ebadd5e63ab86007'))


DOWNLOAD_FILES
}


fetch <- function(name, url) {
    target <- file.path(DIR$raw, name)
    part <- paste0(target, ".part")
    file_name <- name
    info <- DOWNLOAD_FILES[DOWNLOAD_FILES$name == file_name]
    if (nrow(info) != 1L) 
        stop("Unknown download: ", name)
    verified <- function(path) {
        if (!file.exists(path) || file.size(path) != info$bytes) 
            return(FALSE)
        actual <- unname(tools::md5sum(path))
        length(actual) == 1L && !is.na(actual) && identical(actual, info$md5)
    }
    if (verified(target)) {
        message("Already complete; reusing ", name)
        return(target)
    }
    if (file.exists(target)) 
        stop("Existing final file failed verification: ", target, "\nIt has been preserved. Move it aside or replace it with the complete original public file.")
    curl_exe <- unname(Sys.which(if (.Platform$OS.type == "windows") 
        "curl.exe"
    else "curl"))
    if (!nzchar(curl_exe)) 
        stop("This version needs the curl command-line downloader for safe resuming.", "\nWindows normally includes curl.exe. Alternatively download the URL in a browser and save it to: ", 
            target, "\nURL: ", url)
    for (attempt in 1:4) {
        current <- if (file.exists(part)) 
            file.size(part)
        else 0
        if (current > info$bytes) 
            stop("Partial file is larger than the original; preserved for inspection: ", part)
        if (current == info$bytes) {
            if (!verified(part)) 
                stop("Complete-length partial file failed MD5 verification; it has been preserved: ", part)
            if (!file.rename(part, target)) 
                stop("Cannot finalize download: ", target)
            return(target)
        }
        message("\n", name, ": ", round(current/1024^2, 1), " / ", round(info$bytes/1024^2, 1), " MiB already saved; ", if (current > 
            0) 
            "RESUMING"
        else "STARTING", " (attempt ", attempt, "/4).")
        args <- c("--location", "--fail", "--connect-timeout", "30", "--speed-limit", "1024", "--speed-time", "120", "--max-time", 
            "0", "--continue-at", "-", "--progress-bar", "--output", shQuote(part), shQuote(url))
        status <- tryCatch(suppressWarnings(system2(curl_exe, args = args)), error = function(e) {
            message(conditionMessage(e))
            127L
        })
        if (verified(part)) {
            if (!file.rename(part, target)) 
                stop("Cannot finalize download: ", target)
            message("Download complete and verified: ", name)
            return(target)
        }
        if (file.exists(part) && file.size(part) == info$bytes) 
            stop("Downloaded file failed MD5 verification; preserved: ", part)
        if (status == 33L) 
            stop("This server did not allow resuming this request. Partial file preserved: ", part, "\nUse the complete prepared CSV, or download this original file in a browser: ", 
                url)
        if (status %in% c(22L, 23L, 60L, 127L)) 
            stop("Download stopped (curl exit ", status, "). Partial file preserved: ", part, "\nRead the curl error above. Public URL: ", 
                url)
        if (attempt < 4L) {
            message("Transfer interrupted. Keeping all saved bytes; retrying in 3 seconds.")
            Sys.sleep(3)
        }
    }
    stop("Download interrupted repeatedly. The partial file has NOT been deleted or restarted.", "\nRerun this script and select the SAME output folder to resume: ", 
        part, "\nURL: ", url)
}


plain_file <- function(path) {
    if (!grepl("\\.gz$", path)) 
        return(path)
    out <- file.path(DIR$cache, sub("\\.gz$", "", basename(path)))
    if (file.exists(out)) 
        return(out)
    message("Download is complete. Decompressing LOCALLY: ", basename(path))
    part <- paste0(out, ".part")
    input <- gzfile(path, "rb")
    output <- file(part, "wb")
    tryCatch({
        repeat {
            b <- readBin(input, "raw", n = 1024L * 1024L)
            if (!length(b)) 
                break
            writeBin(b, output)
        }
    }, finally = {
        close(input)
        close(output)
    })
    if (!file.rename(part, out)) 
        stop("Cannot complete decompression: ", out)
    out
}


finalize <- function(id, x, meta, feature_map, details) {
    stopifnot(nrow(x) == nrow(meta), ncol(x) == nrow(feature_map), identical(rownames(x), meta$sample_id))
    if (anyNA(meta$sample_id) || anyNA(meta$patient_id) || anyDuplicated(meta$sample_id) || anyDuplicated(meta$patient_id)) 
        stop("Missing or repeated sample/patient IDs: ", id)
    if (anyNA(meta$label) || !all(meta$label %in% 0:1)) 
        stop("Invalid labels")
    if (anyDuplicated(feature_map$feature_id)) 
        stop("Duplicated feature identifiers")
    spec <- DATASETS[dataset == id]
    n1 <- sum(meta$label == 1L)
    n0 <- sum(meta$label == 0L)
    if (nrow(x) != spec$expected_n || ncol(x) != spec$expected_p || n1 != spec$expected_positive || n0/n1 < spec$ir_min || 
        n0/n1 > spec$ir_max) 
        stop("Unexpected source dimensions/class counts for ", id)
    colnames(x) <- feature_map$feature_id
    output <- file.path(DIR$processed, paste0(id, ".csv"))
    message("Writing processed CSV LOCALLY: ", basename(output))
    fwrite(cbind(meta[, .(sample_id, patient_id, label)], as.data.table(x)), output, na = "NA", bom = TRUE)
    message(id, ": n=", nrow(x), "; positive=", n1, "; negative=", n0, "; IR=", round(n0/n1, 4), "; features=", ncol(x), 
        "\nSaved: ", normalizePath(output, winslash = "/"))
    invisible(output)
}


parse_soft_samples <- function(path) {
    con <- gzfile(path, "rt", encoding = "UTF-8")
    on.exit(close(con))
    rows <- list()
    cur <- NULL
    processing <- character()
    repeat {
        a <- readLines(con, n = 1000, warn = FALSE)
        if (!length(a)) 
            break
        for (line in a) {
            if (startsWith(line, "^SAMPLE = ")) {
                if (!is.null(cur)) 
                  rows[[length(rows) + 1L]] <- cur
                cur <- list(sample_id = sub("^\\^SAMPLE = ", "", line))
            }
            else if (!is.null(cur) && startsWith(line, "!Sample_title = ")) 
                cur$title <- sub("^!Sample_title = ", "", line)
            else if (!is.null(cur) && startsWith(line, "!Sample_characteristics_ch1 = ")) {
                s <- sub("^!Sample_characteristics_ch1 = ", "", line)
                at <- regexpr(": ", s, fixed = TRUE)[1]
                if (at > 0) 
                  cur[[substr(s, 1, at - 1)]] <- substr(s, at + 2, nchar(s))
            }
            else if (startsWith(line, "!Sample_data_processing = ")) 
                processing <- c(processing, sub("^!Sample_data_processing = ", "", line))
        }
    }
    if (!is.null(cur)) 
        rows[[length(rows) + 1L]] <- cur
    list(meta = rbindlist(rows, fill = TRUE), processing = unique(processing))
}


series_header <- function(path) {
    con <- gzfile(path, "rt", encoding = "UTF-8")
    on.exit(close(con))
    a <- character()
    repeat {
        line <- readLines(con, n = 1, warn = FALSE)
        if (!length(line)) 
            stop("Missing matrix marker")
        if (startsWith(line, "!series_matrix_table_begin")) 
            break
        if (startsWith(line, "!Sample_")) 
            a <- c(a, line)
    }
    ids <- NULL
    chars <- list()
    proc <- character()
    for (line in a) {
        z <- scan(text = line, what = character(), sep = "\t", quote = "\"", quiet = TRUE)
        if (z[1] == "!Sample_geo_accession") 
            ids <- z[-1]
        if (z[1] == "!Sample_title") 
            titles <- z[-1]
        if (z[1] == "!Sample_characteristics_ch1") 
            chars[[length(chars) + 1L]] <- z[-1]
        if (z[1] == "!Sample_data_processing") 
            proc <- unique(c(proc, z[-1]))
    }
    out <- data.table(sample_id = ids, title = titles)
    for (v in chars) {
        k <- sub(":.*$", "", v[1])
        out[[k]] <- sub("^[^:]*: ?", "", v)
    }
    list(meta = out, processing = proc)
}


read_series_expression <- function(path, keep_samples) {
    con <- gzfile(path, "rt", encoding = "UTF-8")
    on.exit(close(con))
    repeat {
        a <- readLines(con, n = 1, warn = FALSE)
        if (!length(a)) 
            stop("No matrix")
        if (startsWith(a, "!series_matrix_table_begin")) 
            break
    }
    header <- scan(text = readLines(con, n = 1), what = character(), sep = "\t", quote = "\"", quiet = TRUE)
    chunks <- list()
    found_end <- FALSE
    repeat {
        lines <- readLines(con, n = 1000, warn = FALSE)
        if (!length(lines)) 
            break
        if (any(startsWith(lines, "!series_matrix_table_end"))) 
            found_end <- TRUE
        lines <- lines[!startsWith(lines, "!")]
        if (!length(lines)) 
            break
        z <- fread(text = paste(lines, collapse = "\n"), header = FALSE, sep = "\t", quote = "\"", showProgress = FALSE, 
            select = c(1, match(keep_samples, header)))
        setnames(z, c("feature_id", keep_samples))
        chunks[[length(chunks) + 1L]] <- z
    }
    if (!found_end) 
        stop("Incomplete GEO matrix: table_end marker missing. Redownload the complete original file.")
    out <- rbindlist(chunks)
    if (nrow(out) != 2565L) 
        stop("Unexpected miRNA row count: ", nrow(out), " (expected 2565).")
    out
}


prepare_scanb <- function() {
    sf <- fetch("GSE96058_family.soft.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/soft/GSE96058_family.soft.gz")
    ef <- fetch("GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz")
    a <- parse_soft_samples(sf)
    m <- copy(a$meta)
    m[, `:=`(patient_id, sub("\\..*$", "", get("scan-b external id")))]
    m[, `:=`(original_class, paste0("Clinical HER2 status ", get("her2 status")))]
    m[, `:=`(label, match(get("her2 status"), c("0", "1")) - 1L)]
    m[, `:=`(exclusion_reason, ifelse(grepl("repl$", title, ignore.case = TRUE), "technical replicate", ifelse(is.na(label), 
        "missing clinical HER2 status", "")))]
    m[!nzchar(exclusion_reason) & duplicated(patient_id), `:=`(exclusion_reason, "repeated patient")]
    keep <- m[exclusion_reason == ""]
    setorder(keep, title)
    plain <- plain_file(ef)
    z <- fread(plain, select = c(1L, match(keep$title, names(fread(plain, nrows = 0)))), showProgress = TRUE)
    genes <- as.character(z[[1]])
    z[[1]] <- NULL
    if (anyDuplicated(genes) || any(!nzchar(genes))) 
        stop("Unexpected gene symbol representation")
    x <- t(as.matrix(z))
    rownames(x) <- keep$sample_id
    fmap <- data.table(feature_id = paste0("symbol:", genes), gene_symbol = genes, identifier_type = "depositor gene symbol")
    finalize("D1_SCANB_HER2", x, keep, fmap, list(platform = "GPL11154/GPL18573; Illumina RNA sequencing", source_scale = "log2(FPKM + 0.1), already transformed by depositor", 
        additional_log_transform = FALSE, reference_label = "clinical her2 status; NOT PAM50 or expression-derived HER2 predictions", 
        preprocessing = a$processing, exclusions = m[exclusion_reason != "", .N, by = exclusion_reason]))
}


prepare_serum <- function() {
    sf <- fetch("GSE106817_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106817/matrix/GSE106817_series_matrix.txt.gz")
    a <- series_header(sf)
    m <- copy(a$meta)
    m[, `:=`(original_class, sub(" [^ ]+$", "", title))]
    m[, `:=`(label, as.integer(original_class == "Lung Cancer"))]
    m[, `:=`(patient_id, title)]
    setorder(m, sample_id)
    z <- read_series_expression(sf, m$sample_id)
    if (anyNA(z$feature_id)) 
        stop("Missing miRBase identifiers")
    z <- z[grepl("^MIMAT[0-9]+$", feature_id)]
    ids <- as.character(z$feature_id)
    z$feature_id <- NULL
    x <- t(as.matrix(z))
    storage.mode(x) <- "double"
    rownames(x) <- m$sample_id
    fmap <- data.table(feature_id = paste0("mirbase:", ids), mirbase_mature_accession = ids)
    finalize("D2_GSE106817_LUNG", x, m, fmap, list(source = "GSE106817", source_scale = "depositor background-corrected log2 miRNA microarray signals", 
        additional_log_transform = FALSE, reference_label = "Lung Cancer title versus all other specimens", subject_id_limit = "unique specimen title is a sample-level proxy; independent patient IDs unavailable"))
}


prepare_sarcoma <- function() {
    sf <- fetch("GSE211692_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211692/matrix/GSE211692_series_matrix.txt.gz")
    ef <- fetch("GSE211692_processed_data_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211692/suppl/GSE211692_processed_data_matrix.txt.gz")
    df <- fetch("GSE106817_series_matrix.txt.gz", "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106817/matrix/GSE106817_series_matrix.txt.gz")
    m <- series_header(sf)$meta
    d2 <- series_header(df)$meta
    excluded <- sub("^.* ", "", d2$title)
    m <- m[!title %in% excluded]
    m[, `:=`(patient_id, title)]
    m[, `:=`(original_class, get("disease state"))]
    m[, `:=`(label, as.integer(original_class == "bone and soft tissue sarcomas"))]
    setorder(m, sample_id)
    plain <- plain_file(ef)
    header <- names(fread(plain, nrows = 0, check.names = FALSE))
    if (anyNA(match(m$title, header))) 
        stop("Missing expression columns for retained source specimens")
    z <- fread(plain, select = c(1L, match(m$title, header)), check.names = FALSE, showProgress = TRUE)
    ids <- as.character(z[[1]])
    keep <- grepl("^hsa-", ids) & !grepl(",", ids, fixed = TRUE)
    z <- z[keep]
    ids <- ids[keep]
    z[[1]] <- NULL
    x <- t(as.matrix(z))
    storage.mode(x) <- "double"
    rownames(x) <- m$sample_id
    if (anyDuplicated(as.data.table(x))) 
        stop("Duplicated retained expression profiles; inspect sample identity before proceeding")
    message("D3: excluded ", sum(series_header(sf)$meta$title %in% excluded), " known D2 source IDs.")
    finalize("D3_GSE211692_SARCOMA", x, m, data.table(feature_id = paste0("mirna:", ids)), list(source = "GSE211692", label = "bone and soft tissue sarcomas versus other retained specimens", 
        known_D2_source_overlap_removed = 3906L, source_scale = "depositor internally normalized signal; no extra transform"))
}



prepare_data <- function(data_dir = DATA_DIR, datasets = DATASET_IDS) {
  initialize_dependencies(analysis = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  root <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  DIR <<- list(raw = file.path(root, "raw"), cache = file.path(root, "cache"), processed = root)
  invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))
  DOWNLOAD_FILES <<- download_files_manifest()
  prepare <- list(D1 = prepare_scanb, D2 = prepare_serum, D3 = prepare_sarcoma)
  if (any(!datasets %in% names(prepare))) stop("Use dataset IDs D1, D2 and/or D3.")
  for (id in datasets) {
    path <- file.path(root, paste0(DATASET_SPECS$dataset[match(id, DATASET_SPECS$id)], ".csv"))
    if (file.exists(path)) message("Processed CSV already exists: ", path)
    else prepare[[id]]()
    invisible(gc())
  }
  message("Processed CSV files are in ", root)
  invisible(root)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  experiment <- as.character(EXPERIMENT); data_dir <- DATA_DIR; output_dir <- OUTPUT_DIR
  datasets <- DATASET_IDS; download <- FALSE; download_only <- FALSE
  i <- 1L
  while (i <= length(args)) {
    flag <- args[i]
    if (flag == "--help") {
      cat("Rscript code.R --experiment <1|2|3|all> --data-dir <folder> --output-dir <folder>\n",
          "Optional: --cpus <n> --datasets D1,D2,D3 --download --download-only\n", sep = "")
      return(invisible(NULL))
    }
    if (flag %in% c("--download", "--download-only")) {
      download <- TRUE; download_only <- flag == "--download-only" || download_only
    } else {
      if (i == length(args)) stop("Missing value for ", flag)
      value <- args[i + 1L]; i <- i + 1L
      if (flag == "--experiment") experiment <- value
      else if (flag == "--data-dir") data_dir <- value
      else if (flag == "--output-dir") output_dir <- value
      else if (flag == "--datasets") datasets <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
      else if (flag == "--cpus") {
        budget <- suppressWarnings(as.integer(value))
        if (is.na(budget) || budget < 1L) stop("--cpus must be a positive integer.")
        CPU_BUDGET <<- budget
      } else stop("Unknown option: ", flag)
    }
    i <- i + 1L
  }
  if (!experiment %in% c("1", "2", "3", "all")) stop("--experiment must be 1, 2, 3 or all.")
  if (download) prepare_data(data_dir, datasets)
  if (download_only) return(invisible(NULL))
  experiments <- if (experiment == "all") 1:3 else as.integer(experiment)
  for (e in experiments) run_study(e, data_dir, output_dir, datasets)
  invisible(NULL)
}

if (!isTRUE(getOption("fr.portable.load_only", FALSE))) main()
