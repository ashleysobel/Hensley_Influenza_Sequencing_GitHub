# 1) Prepare environment --------------------------------------------------------------
## 1.1) Clear the R environment
rm(list = ls())

## 1.2) Auto‑set working directory to this script’s folder
get_script_path <- function(){
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    normalizePath(sub("^--file=", "", file_arg))
  } else if (requireNamespace("rstudioapi", quietly = TRUE) &&
             rstudioapi::isAvailable()) {
    normalizePath(rstudioapi::getSourceEditorContext()$path)
  } else {
    stop("Unable to detect script location; please run via Rscript or RStudio.")
  }
}
script_dir <- dirname(get_script_path())
setwd(script_dir)

## 1.3) Install and load required libraries
required_pkgs <- c("parallel", "renv", "rstudioapi")
installed_pkgs <- rownames(installed.packages())
missing_pkgs   <- setdiff(required_pkgs, installed_pkgs)
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))


# 2) Define user options -----------------------------------------------------------
## 2.1) Path to the file storing your project directory
project_config <- "project_path.txt"
if (!file.exists(project_config)) {
  stop(
    "Could not find 'project_path.txt' in the project directory.\n",
    "Please create this file and enter the full path to your project root on the first line."
  )
}
user_directory <- trimws(readLines(project_config, n = 1))
if (!dir.exists(user_directory)) {
  stop("Project directory does not exist: ", user_directory)
}
setwd(user_directory)


## 2.2) Establish list of run IDs (six‑character IDs as folder names)
# Provide your RunIDs here to override auto‑detection. 
# e.g. runs <- c("DG2BND", "WPGFLD")
runs <- c()  

# If no RunIDs given, auto‑detect them from your Dino_Fastq folders:
if (length(runs) == 0) {
  fastq_base <- file.path(user_directory, "Dino_Fastq")
  # List subfolders named Dino_<RunID>
  dinos <- list.files(fastq_base, pattern = "^Dino_[A-Za-z0-9]{6}$", full.names = FALSE)
  # Strip the "Dino_" prefix
  runs <- sub("^Dino_", "", dinos)
  if (length(runs) == 0) {
    stop("No Run IDs provided and no Dino_<RunID> folders found under ", fastq_base)
  }
  message(sprintf(
    "No RunIDs specified; auto‑detected %d runs: %s",
    length(runs), paste(runs, collapse = ", ")
  ))
}




## 2.3) Coverage threshold for masking consensus positions
min_cov <- 5

## 2.4) Parallel processing flag
use_multicore <- TRUE  # TRUE = multicore; FALSE = one sample at a time


# 3) Generate package citations ----------------------------------------------------
all_citations <- sapply(required_pkgs, function(pkg) utils::toBibtex(citation(pkg)))
writeLines(
  paste(all_citations, collapse = "\n\n"),
  file.path(user_directory, "pipeline_processing_citations.bib")
)


# 4) Create required directories ---------------------------------------------------
## 4.1) Ensure base output & FASTQ folders exist
dir_list <- c(
  file.path(user_directory, "IRMA_output"),
  file.path(user_directory, "Dino_Fastq")
)
for (d in dir_list) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## 4.2) Verify FASTQ subfolder structure
fastq_base   <- file.path(user_directory, "Dino_Fastq")
missing_dirs <- character(0)
if (!dir.exists(fastq_base)) {
  dir.create(fastq_base, recursive = TRUE)
  missing_dirs <- c(missing_dirs, fastq_base)
}
for (run_id in runs) {
  subdir <- file.path(fastq_base, paste0("Dino_", run_id))
  if (!dir.exists(subdir)) {
    dir.create(subdir, recursive = TRUE)
    missing_dirs <- c(missing_dirs, subdir)
  }
}
if (length(missing_dirs) > 0) {
  message(
    "Created the following missing directories:\n",
    paste(missing_dirs, collapse = "\n"),
    "\nPlease add FASTQ files and rerun."
  )
  quit(save = "no", status = 0)
}

## 4.3) Session info directory
session_dir <- file.path(user_directory, "session_info")
if (!dir.exists(session_dir)) {
  dir.create(session_dir, recursive = TRUE)
}


# 5) Define IRMA executable path ---------------------------------------------------
config_file <- file.path(user_directory, "irma_path.txt")
if (!file.exists(config_file)) {
  stop(
    "Could not find 'irma_path.txt' in the project directory.\n",
    "Please create this file and enter the full path to your IRMA executable on the first line."
  )
}
irma_exec <- trimws(readLines(config_file, n = 1))
if (!file.exists(irma_exec)) {
  stop("IRMA executable not found at: ", irma_exec)
}
if (file.access(irma_exec, mode = 1) != 0) {
  stop(
    sprintf(
      "IRMA at '%s' is not executable.\nFix with:\n  chmod +x '%s'",
      irma_exec, irma_exec
    )
  )
}


# 6) Function to process one run ------------------------------------------------
# 'library_id' (e.g. "6HW5TB"), 'min_cov' minimum coverage threshold for masking.
process_run <- function(library_id, min_cov, fastq_dir, log_file = NULL, use_parallel = TRUE) {
  message("\nProcessing run: ", library_id)
  
  # 1) Define output directory for this run
  output_base <- file.path(user_directory, "IRMA_output", library_id)
  if (!dir.exists(fastq_dir)) stop("FASTQ directory not found: ", fastq_dir)
  if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)
  
  # 2) List and name FASTQ files
  fastq_paths <- list.files(
    fastq_dir,
    pattern    = paste0("^", library_id, "_.*\\.fastq$"),
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
  num_cores <- max(1, parallel::detectCores() - 1)
  if (use_parallel) {
    mclapply(sample_names, process_sample, mc.cores = num_cores)
  } else {
    lapply(sample_names, process_sample)
  }
  
  # 7) Build and return coverage summary (unchanged)
  summary_list <- vector("list", length(sample_names))
  for (j in seq_along(sample_names)) {
    samp <- sample_names[j]
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
  
  # 8) Ensure all entries are data.frames
  for (j in seq_along(summary_list)) {
    if (!is.data.frame(summary_list[[j]])) {
      summary_list[[j]] <- data.frame(
        sample               = sample_names[j],
        run_id               = library_id,
        strain               = NA_character_,
        avg_coverage         = NA_real_,
        prop_low_depth       = NA_real_,
        ha_coverage_prop_low = NA_real_,
        na_coverage_prop_low = NA_real_,
        stringsAsFactors     = FALSE
      )
    }
  }
  
  # 9) Return combined summary
  final_summary <- do.call(rbind, summary_list)
  return(final_summary)
}




# 7) Run pipeline across multiple runs -------------------------------------------

# 7.1) Loop over specified runs
run_summaries <- vector("list", length(runs))
r <- 1
for (r in seq_along(runs)) {
  run_id <- runs[r]
  fastq_dir <- file.path(
    user_directory,
    "Dino_Fastq",
    paste0("Dino_", run_id)
  )
  
  # Create a run-specific log file in session_dir
  log_file <- file.path(session_dir, paste0("missing_segments_", run_id, ".log"))
  file.create(log_file)  # Overwrite any old file
  
  run_summary <- tryCatch(
    process_run(run_id, min_cov, fastq_dir, log_file = log_file),
    error = function(e) {
      warning(sprintf("Processing failed for run %s: %s", run_id, e$message))
      # Return an empty data.frame to preserve length
      return(data.frame(
        sample = NA_character_,
        run_id = run_id,
        strain = NA_character_,
        avg_coverage = NA_real_,
        prop_low_depth = NA_real_,
        ha_coverage_prop_low = NA_real_,
        na_coverage_prop_low = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
  )
  
  # After run_summary is generated for this run
  # (Assume run_id, run_summary, session_dir are already defined in your loop)
  
  log_file <- file.path(session_dir, paste0("missing_segments_", run_id, ".log"))
  
  # Function to parse missing segments log and return a list: sample_name -> vector of missing segments or "All"
  parse_missing_segments_log <- function(log_file) {
    if (!file.exists(log_file)) return(list())
    lines <- readLines(log_file)
    sample_missing <- list()
    for (i in seq_along(lines)) {
      line <- lines[i]
      if (nchar(line) == 0) next
      # Find all positions where phrases start
      matches <- gregexpr("Missing consensus for segment [A-Za-z0-9]+ in sample [^ ]+|No consensus segments for sample [^ ]+", line, perl=TRUE)
      spans <- regmatches(line, matches)[[1]]
      for (msg in spans) {
        msg <- trimws(msg)
        if (grepl("^No consensus segments for sample", msg)) {
          sample <- sub("^No consensus segments for sample (.+)$", "\\1", msg)
          sample_missing[[sample]] <- "All"
        }
        if (grepl("^Missing consensus for segment", msg)) {
          matches <- regmatches(
            msg, regexec("^Missing consensus for segment ([A-Za-z0-9]+) in sample (.+)$", msg)
          )[[1]]
          if (length(matches) == 3) {
            seg <- matches[2]
            sample <- matches[3]
            if (is.null(sample_missing[[sample]])) {
              sample_missing[[sample]] <- seg
            } else if (!identical(sample_missing[[sample]], "All")) {
              sample_missing[[sample]] <- c(sample_missing[[sample]], seg)
            }
          }
        }
      }
    }
    return(sample_missing)
  }
  
  # Parse the log
  sample_missing_list <- parse_missing_segments_log(log_file)
  
  # Add columns to run_summary
  if (!is.null(run_summary) && is.data.frame(run_summary) && nrow(run_summary) > 0) {
    # Make new columns (default: no missing segments)
    run_summary$No_missing_seg <- 0L
    run_summary$missing_seg <- ""
    # For each sample, fill in the log-based values
    for (i in seq_len(nrow(run_summary))) {
      sample_name <- run_summary$sample[i]
      missing <- sample_missing_list[[sample_name]]
      if (!is.null(missing)) {
        if (identical(missing, "All")) {
          run_summary$No_missing_seg[i] <- 8L
          run_summary$missing_seg[i] <- "All"
        } else {
          run_summary$No_missing_seg[i] <- length(missing)
          run_summary$missing_seg[i] <- paste(missing, collapse = ", ")
        }
      }
      # If missing is NULL, these remain at 0 and "" (i.e. all segments present)
    }
    
    # FINAL OVERRIDE for samples where all segments are missing based on avg_coverage
    idx_missing_all <- which(is.na(run_summary$avg_coverage))
    if (length(idx_missing_all) > 0) {
      run_summary$No_missing_seg[idx_missing_all] <- 8L
      run_summary$missing_seg[idx_missing_all] <- "All"
    }
  }
  
  
  
  run_summaries[[r]] <- run_summary
  
  # Only save CSV if run_summary exists and is a non-empty data.frame
  if (!is.null(run_summary) && is.data.frame(run_summary) && nrow(run_summary) > 0 && !all(is.na(run_summary$sample))) {
    out_file <- file.path(user_directory, paste0(run_id, "_coverage_summary.csv"))
    write.csv(run_summary, file = out_file, row.names = FALSE)
    message("Saved summary for run ", run_id, " to ", out_file)
  } else {
    message(sprintf("No summary generated for run %s (see warnings above).", run_id))
    
  }
}

# 7.2) Set names so length always matches, even if some runs failed completely
names(run_summaries) <- runs

# 7.3) Generate robust summary block for the specified runs
message("\nAll runs processed. Summary:")
for (run_id in runs) {
  rs <- run_summaries[[run_id]]
  if (is.null(rs) || !is.data.frame(rs) || all(is.na(rs$sample))) {
    n_samples <- 0
    n_failures <- 0
  } else {
    n_samples  <- nrow(rs)
    n_failures <- sum(is.na(rs$strain))
  }
  message(sprintf("  %s: %d samples, %d failures", run_id, n_samples, n_failures))
}
invisible(run_summaries)



# 8) Processing summary ----------------------------------------------------------

# Provide an overall summary so that users sourcing this script see results at a glance
message("\nAll runs processed. Summary:")
names(run_summaries) <- runs
for (run_id in runs) {
  rs <- run_summaries[[run_id]]
  n_samples  <- if(is.data.frame(rs)) nrow(rs) else length(rs)
  n_failures <- if(is.data.frame(rs) && "strain" %in% names(rs)) sum(is.na(rs$strain)) else 0
  message(sprintf("  %s: %d samples, %d failures", run_id, n_samples, n_failures))
}
# Return invisibly for programmatic use
invisible(run_summaries)



# 9) Store session information ---------------------------------------------------

# 9.1) Create session directory and save sessionInfo()
session_dir <- file.path(user_directory, "session_info")
if (!dir.exists(session_dir)) {
  dir.create(session_dir, recursive = TRUE)
  message("Created session info directory: ", session_dir)
}

capture.output(
  sessionInfo(),
  file = file.path(session_dir, "sessionInfo.txt")
)
message("Wrote session information to sessionInfo.txt")

# 9.2) Write a table of all installed packages and their versions
installed <- installed.packages()[, c("Package", "Version")]
write.table(
  installed,
  file      = file.path(session_dir, "installed_packages.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("Wrote installed packages to installed_packages.tsv")

# 9.3) Write a recursive listing of all project files
all_files <- list.files(
  path      = user_directory,
  recursive = TRUE,
  all.files = TRUE,
  no..      = TRUE
)
writeLines(
  all_files,
  con = file.path(session_dir, "project_files.txt")
)
message("Wrote project file listing to project_files.txt")

# 9.4) Update renv lockfile to capture current environment
if ("renv" %in% rownames(installed.packages())) {
  renv::snapshot(prompt = FALSE)
  message("Updated renv lockfile with current package environment.")
}

# 10) Compile individual summaries & add pass/fail flags -----


# 10.1) Combine all run_summaries into one data.frame
compiled <- do.call(rbind, run_summaries)

# 10.2) Add Pass_All, Pass_HA, Pass_NA columns
compiled$Pass_All <- ifelse(
  compiled$No_missing_seg == 0 & compiled$prop_low_depth < 0.0105,
  "Y", "N"
)

# To pass HA: must not list HA among missing_seg AND ha_coverage_prop_low < 0.0105
compiled$Pass_HA <- ifelse(
  !grepl("\\bHA\\b", compiled$missing_seg) & compiled$ha_coverage_prop_low < 0.0105,
  "Y", "N"
)

# To pass NA: must not list NA among missing_seg AND na_coverage_prop_low < 0.0105
compiled$Pass_NA <- ifelse(
  !grepl("\\bNA\\b", compiled$missing_seg) & compiled$na_coverage_prop_low < 0.0105,
  "Y", "N"
)

# 10.3) Write out the compiled summary
out_compiled <- file.path(user_directory, "compiled_coverage_summary.csv")
write.csv(compiled, out_compiled, row.names = FALSE)
message("Wrote compiled coverage summary to: ", out_compiled)



