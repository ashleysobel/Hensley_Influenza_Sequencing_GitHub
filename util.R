# 6) Function to process one run ------------------------------------------------

# 'library_id' (e.g. "6HW5TB"), 'min_cov' minimum coverage threshold for masking.
process_run <- function(library_id, min_cov, fastq_dir, log_file = NULL, use_parallel = TRUE, irma_exec) {
  
  
  # 1) Define output directory for this run
  output_base <- file.path(user_directory, "IRMA_output", library_id)
  if (!dir.exists(fastq_dir)) stop("FASTQ directory not found: ", fastq_dir)
  if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)
  
  # 2) List and name FASTQ files
  fastq_paths <- base::list.files(
    path = fastq_dir,
    full.names = TRUE
  )
  
  if (length(fastq_paths) == 0) stop("No FASTQ files found in: ", fastq_dir)
  orig        <- basename(fastq_paths)                                 # e.g. "PJ6FV5_1_CHOI-005_N1.fastq"
  no_ext      <- tools::file_path_sans_ext(orig)                       # "PJ6FV5_1_CHOI-005_N1"
  sample_id   <- sub(paste0("^", library_id, "_[0-9]+_"), "", no_ext)  # "CHOI-005_N1"
  sample_names <- paste0(sample_id, "_", library_id)                   # "CHOI-005_N1_PJ6FV5"
  sample_success <- logical(length(sample_names))
  
  # 3) Run IRMA on each sample
  for (i in seq_along(fastq_paths)) {
    fq   <- fastq_paths[i]
    samp <- sample_names[i]
    outd <- file.path(output_base, samp)
    
    message(sprintf(">>> [%d/%d] IRMA on %s", i, length(fastq_paths), samp))
    message("\nProcessing run: ", fq)
    
    res <- system2(
      irma_exec,
      args   = c("FLU", fq, outd),
      stdout = TRUE,
      stderr = TRUE
    )
    cat(res, sep = "\n")
    sample_success[i] <- dir.exists(file.path(outd, "amended_consensus"))
  }
  
  # 4) Prepare consensus output directories
  seg_abbr    <- c("PB2","PB1","PA","HA","NP","NA","MP","NS")
  masked_dir   <- file.path(user_directory, "IRMA_output", "Masked_consensus")
  unmasked_dir <- file.path(user_directory, "IRMA_output", "Unmasked_consensus")
  dir.create(masked_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(unmasked_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Create HA/NA subdirectories under masked_consensus
  ha_dir <- file.path(masked_dir, "HA")
  na_dir <- file.path(masked_dir, "NA")
  dir.create(ha_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(na_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 5) Inner function to process one sample
  process_sample <- function(samp) {
    amended_dir <- file.path(output_base, samp, "amended_consensus")
    seg_files <- list.files(
      amended_dir,
      pattern    = paste0("^", samp, "_[0-9]+\\.fa$"),
      full.names = TRUE
    )
    if (length(seg_files) == 0) {
      warn_msg <- sprintf("No consensus segments for sample %s", samp)
      warning(warn_msg)
      if (!is.null(log_file)) cat(warn_msg, "\n", file = log_file, append = TRUE)
      return(NULL)
    }
    
    # Determine segment indices & names
    seg_index <- as.numeric(sub(
      paste0("^", samp, "_([0-9]+)\\.fa$"), "\\1", basename(seg_files)
    ))
    ord       <- order(seg_index)
    seg_files <- seg_files[ord]
    seg_index <- seg_index[ord]
    seg_names <- seg_abbr[seg_index]
    
    # Log any missing segments
    missing_idx <- setdiff(seq_along(seg_abbr), seg_index)
    if (length(missing_idx) > 0) {
      for (idx in missing_idx) {
        warn_msg <- sprintf("Missing consensus for segment %s in sample %s", seg_abbr[idx], samp)
        warning(warn_msg)
        if (!is.null(log_file)) cat(warn_msg, "\n", file = log_file, append = TRUE)
      }
    }
    
    # Read coverage tables
    cov_dir   <- file.path(output_base, samp, "tables")
    cov_files <- list.files(cov_dir, pattern = "-coverage\\.txt$", full.names = TRUE)
    if (length(cov_files) == 0) {
      warn_msg <- sprintf("No coverage tables for sample %s", samp)
      warning(warn_msg)
      if (!is.null(log_file)) cat(warn_msg, "\n", file = log_file, append = TRUE)
      return(NULL)
    }
    
    # Collect coverage depths per segment
    seg_cov_depths <- vector("list", length(seg_index))
    for (k in seq_along(seg_index)) {
      gene    <- seg_names[k]
      pat     <- if (gene %in% c("HA","NA")) {
        sprintf("A_%s_.*-coverage\\.txt$", gene)
      } else {
        sprintf("%s.*-coverage\\.txt$", gene)
      }
      match_f <- cov_files[grepl(pat, basename(cov_files))]
      if (length(match_f) == 0) {
        warn_msg <- sprintf("Missing coverage for segment %s in sample %s", gene, samp)
        warning(warn_msg)
        if (!is.null(log_file)) cat(warn_msg, "\n", file = log_file, append = TRUE)
        seg_cov_depths[[k]] <- numeric(0)
      } else {
        df <- read.delim(match_f[1], check.names = FALSE)
        seg_cov_depths[[k]] <- df[["Coverage Depth"]]
      }
    }
    coverage_vector <- unlist(seg_cov_depths, use.names = FALSE)
    if (length(coverage_vector) == 0) {
      warn_msg <- sprintf("No usable coverage for sample %s", samp)
      warning(warn_msg)
      if (!is.null(log_file)) cat(warn_msg, "\n", file = log_file, append = TRUE)
      return(NULL)
    }
    
    # Open full masked & unmasked FASTA files
    masked_fa   <- file.path(masked_dir,   paste0(samp, "_masked.fasta"))
    unmasked_fa <- file.path(unmasked_dir, paste0(samp, "_unmasked.fasta"))
    con_masked  <- file(masked_fa,   open = "w");    on.exit(close(con_masked), add = TRUE)
    con_unmasked<- file(unmasked_fa, open = "w");    on.exit(close(con_unmasked), add = TRUE)
    
    # Open HA & NA masked FASTA files
    con_ha <- file(file.path(ha_dir, paste0(samp, "_masked_HA.fasta")), open = "w")
    on.exit(close(con_ha), add = TRUE)
    con_na <- file(file.path(na_dir, paste0(samp, "_masked_NA.fasta")), open = "w")
    on.exit(close(con_na), add = TRUE)
    
    # Write sequences for each segment
    pos_idx <- 1
    for (k in seq_along(seg_files)) {
      lines     <- readLines(seg_files[k])
      seq_vec   <- strsplit(paste0(lines[-1], collapse = ""), "")[[1]]
      len_seg   <- length(seq_vec)
      cov_seg   <- coverage_vector[pos_idx:(pos_idx + len_seg - 1)]
      masked_seq<- seq_vec
      masked_seq[cov_seg < min_cov] <- "N"
      
      # Build header with abbreviation
      header <- paste0(">", samp, "_", seg_names[k])
      
      # Write to full masked & unmasked files
      writeLines(header,                      con_masked)
      writeLines(paste(masked_seq, collapse=""), con_masked)
      writeLines(header,                      con_unmasked)
      writeLines(paste(seq_vec,     collapse=""), con_unmasked)
      
      # If HA or NA, also write to dedicated file
      if (seg_names[k] == "HA") {
        writeLines(header,                      con_ha)
        writeLines(paste(masked_seq, collapse=""), con_ha)
      }
      if (seg_names[k] == "NA") {
        writeLines(header,                      con_na)
        writeLines(paste(masked_seq, collapse=""), con_na)
      }
      
      pos_idx <- pos_idx + len_seg
    }
    
    message(sprintf("▶ FASTAs written for sample %s", samp))
    return(samp)
  }
  
  # 6) Run sample processing (leave 1 core free)
  if (use_parallel) {
    message("using parallel")
    num_cores <- max(1, parallel::detectCores() - 1)
    mclapply(sample_names, process_sample, mc.cores = num_cores)
  } else {
    lapply(sample_names, process_sample)
  }
  
  # 7) Build and return coverage summary (unchanged)…
  summary_list <- vector("list", length(sample_names))
  for (j in seq_along(sample_names)) {
    samp <- sample_names[[j]]
    if (!sample_success[j]) {
      summary_list[[j]] <- data.frame(
        sample               = samp,
        run_id               = library_id,
        strain               = NA_character_,
        avg_coverage         = NA_real_,
        prop_low_depth       = NA_real_,
        ha_coverage_prop_low = NA_real_,
        na_coverage_prop_low = NA_real_,
        stringsAsFactors     = FALSE
      )
      next
    }
    cov_dir   <- file.path(output_base, samp, "tables")
    cov_files <- list.files(cov_dir, pattern = "-coverage\\.txt$", full.names = TRUE)
    depths    <- unlist(lapply(cov_files, function(cf) {
      df <- read.delim(cf, check.names = FALSE)
      df[["Coverage Depth"]]
    }))
    if (length(depths)>0 && is.numeric(depths)) {
      avg_cov  <- mean(depths, na.rm = TRUE)
      prop_low <- mean(depths < min_cov, na.rm = TRUE)
    } else {
      avg_cov  <- NA_real_; prop_low <- NA_real_
    }
    ha_cov_files <- list.files(cov_dir, pattern = "A_HA_.*-coverage\\.txt$", full.names = TRUE)
    na_cov_files <- list.files(cov_dir, pattern = "A_NA_.*-coverage\\.txt$", full.names = TRUE)
    ha_depths <- unlist(lapply(ha_cov_files, function(f) {
      df <- read.delim(f, check.names = FALSE); df[["Coverage Depth"]]
    }))
    ha_prop <- if (length(ha_depths)>0) mean(ha_depths >= 10, na.rm = TRUE) else NA_real_
    na_depths <- unlist(lapply(na_cov_files, function(f) {
      df <- read.delim(f, check.names = FALSE); df[["Coverage Depth"]]
    }))
    na_prop <- if (length(na_depths)>0) mean(na_depths >= 10, na.rm = TRUE) else NA_real_
    ha_files <- list.files(file.path(output_base, samp), pattern = "_HA_H[0-9]+\\.bam$", full.names = FALSE)
    na_files <- list.files(file.path(output_base, samp), pattern = "_NA_N[0-9]+\\.bam$", full.names = FALSE)
    ha_tag <- if (length(ha_files)) sub(".*_HA_(H[0-9]+)\\..*", "\\1", ha_files[1]) else NA_character_
    na_tag <- if (length(na_files)) sub(".*_NA_(N[0-9]+)\\..*", "\\1", na_files[1]) else NA_character_
    strain <- if (!is.na(ha_tag) && !is.na(na_tag)) paste0(ha_tag, na_tag) else NA_character_
    summary_list[[j]] <- data.frame(
      sample               = samp,
      run_id               = library_id,
      strain               = strain,
      avg_coverage         = avg_cov,
      prop_low_depth       = prop_low,
      ha_coverage_prop_low = 1 - ha_prop,
      na_coverage_prop_low = 1 - na_prop,
      stringsAsFactors     = FALSE
    )
  }

  # 9) Return combined summary
  final_summary <- do.call(rbind, summary_list)
  return(final_summary)
}
