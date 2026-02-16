#' Compute a weighted score from betas (File A) and abundances (File B) of significant targets
#' Targets are standardized before weighting and combining
#' A mean of the weighted targets is used (with people with missing values included via na.rm=TRUE)
#' THIS SCRIPT DOES NOT SELECT BY P-VALUE THRESHOLDS
#' 


#----------------------Build function -----------------------------------------#
#path_betas = "/media/Analyses/CARDIA-ADPQS-replication/Data/diet-MESA-lasso-coefs-2026-02-06.csv"
#abund_df    = tar_read(Proteins_long)
#beta_col    = "lasso_APDQS",
#id_col      = "idno",
#time_col    = "Exam",            
#score_name  = "APDQS_protein_score",
#verbose     = TRUE,
#na_rm       = TRUE,            # <-- robust to NA by exam
#min_non_missing = 3,
#metabolite_col   = "OlinkID"
    
    
    
#---------------------------------------------------------------------------#   
#---------------------------------------------------------------------------#     
    
    

build_weighted_score_from_LASSO <- function(
    path_betas,
    abund_df,
    id_col,
    metabolite_col,
    beta_col,
    time_col,
    score_name,
    normalize_names  = FALSE,
    verbose = TRUE,
    na_rm = TRUE,              # ignore NA abundances (treat as zero contribution)
    min_non_missing       # require at least this many non-missing metabolites per row
) {
  
  
  betas_df <- read.csv(path_betas)
  #----------------Helper -----------------------------#

  
  msg <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  
  # checks
  stopifnot(id_col %in% names(abund_df))
  if (!is.null(time_col)) stopifnot(time_col %in% names(abund_df))
  for (nm in c(metabolite_col, beta_col)) {
    if (!nm %in% names(betas_df)) stop(sprintf("Column '%s' not found in betas_df.", nm))
  }
  
  
  
  
  
  #----------------Select betas and Ps -----------------------------#

  
  df_b <- betas_df[, c(metabolite_col, beta_col)]
  names(df_b) <- c("target", "beta")
  
  
  #----------------Select hits -----------------------------#
  df_b <- df_b[df_b$beta !=0 & !is.na(df_b$beta), , drop = FALSE]
  if (nrow(df_b) == 0L) stop("No metabolites pass the p-value threshold.")
  
  
  if (nrow(df_b) == 0L) stop("No metabolites selected by LASSO")
  # QC info
  Sig_targets_n <- dim(df_b)[1]
  Sig_targets_names <- as.data.frame(df_b$target)
  
  #----------------Check for duplicates in raw data -----------------------------#
  
  dupe_counts <- as.data.frame(table(df_b$target), stringsAsFactors = FALSE)
  names(dupe_counts) <- c("target", "n_rows"); dupe_counts$n_rows <- as.integer(dupe_counts$n_rows)
  n_with_duplicates <- sum(dupe_counts$n_rows > 1L)
  n_total_dupe_rows <- sum(pmax(dupe_counts$n_rows - 1L, 0L))
  if (n_with_duplicates > 0L) {
    dup_preview <- dupe_counts[dupe_counts$n_rows > 1L, ]
    dup_preview <- dup_preview[order(-dup_preview$n_rows, dup_preview$target), ]
    preview_str <- paste(utils::head(sprintf("%s (n=%d)", dup_preview$target, dup_preview$n_rows), 10L), collapse = "; ")
    msg("Found %d metabolites with duplicate rows",
        n_with_duplicates, n_total_dupe_rows)
    msg("Duplicates (up to 10): %s%s", preview_str, if (nrow(dup_preview) > 10L) " ..." else "")
  }
  
  # collapse dupes
  # Removed from this analysis
  #df_b <- aggregate(beta ~ target, data = df_b, FUN = function(z) mean(z, na.rm = TRUE))
  
  #----------------Make scores -----------------------------#
  target_cols_in_abund <- setdiff(names(abund_df), c(id_col, time_col))
  common_targets <- intersect(df_b$target, target_cols_in_abund)
  
  
  missing_targets <- setdiff(df_b$target, target_cols_in_abund)
  nonsig_targets   <- setdiff(target_cols_in_abund, df_b$target)
  
  msg("Selected %d metabolites; %d present; %d missing; %d extra.",
      nrow(df_b), length(common_targets), length(missing_targets), length(nonsig_targets))
  if (length(missing_targets)) {
    msg("Missing (up to 10): %s%s",
        paste(utils::head(missing_targets, 10L), collapse = ", "),
        if (length(missing_targets) > 10) " ..." else "")
  }
  if (!length(common_targets)) stop("None of the selected metabolites are present in abund_df.")
  
  # align to common_targets 
  df_b_common <- df_b[match(common_targets, df_b$target), , drop = FALSE]
  X <- as.matrix(abund_df[, common_targets, drop = FALSE]) #is abundances
  w <- as.numeric(df_b_common$beta) #is betas
  X2 <-scale(X) #Center and scale before use
  
  if (na_rm) {
    # elementwise multiply scaled vars, then rowMeans with na.rm=TRUE
    WX <- t( t(X2) * w )
    score <- rowMeans(WX, na.rm = TRUE)
    
    
    # optional: guard against rows that are almost entirely missing (below missing threshold set in function)
    n_nonmiss <- rowSums(!is.na(X))
    too_sparse <- n_nonmiss < min_non_missing
    if (any(too_sparse)) {
      msg("Warning: %d rows had < %d non-missing metabolites; setting score to NA for those rows.",
          sum(too_sparse), min_non_missing)
      score[too_sparse] <- NA_real_
      too_sparse_excl <- abund_df[too_sparse,]
      too_sparse_excl <- c("spacer", too_sparse_excl[,names(too_sparse)==id_col | names(too_sparse)==time_col])
    }
    
  } else {
    # any NA -> NA
    score <- as.numeric(X %*% w)
  }
  
  
  # build output: include id and (optionally) time
  if (requireNamespace("rlang", quietly = TRUE)) {
    out <- tibble::tibble(
      !!rlang::sym(id_col) := abund_df[[id_col]],
      !!rlang::sym(score_name) := score
    )
    if (!is.null(time_col)) {
      out[[time_col]] <- abund_df[[time_col]]
      out <- out[, c(id_col, time_col, score_name)]
    }
  } else {
    out <- tibble::tibble(tmp_id = abund_df[[id_col]], tmp_score = score)
    names(out) <- c(id_col, score_name)
    if (!is.null(time_col)) {
      out[[time_col]] <- abund_df[[time_col]]
      out <- out[, c(id_col, time_col, score_name)]
    }
  }
  
  
  #----------------Score info -----------------------------#
  Score_info = list(Sig_targets_n = Sig_targets_n,
                    Sig_targets_names = Sig_targets_names,
                    Duplicate_targets_n = n_with_duplicates,
                    Included_targets = common_targets,
                    Missing_targets = missing_targets,
                    min_non_missing = min_non_missing)
  
  #----------------All Outputs -----------------------------#
  list(score_info = Score_info, 
       scores = out)
}