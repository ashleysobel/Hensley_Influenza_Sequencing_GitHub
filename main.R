# Hensley IRMA Pipeline — CLI Version (Step 1: flags, no hardcoded params)
# ---------------------------------------------------------------------
# This refactor replaces hardcoded parameters with command‑line flags using optparse.
# Example usage:
#   Rscript irma_pipeline.R \
#     --project /path/to/project \
#     --runs PJ6FV5,BZPB3P,XLHBTH \
#     --min_cov 5 \
#     --multicore --ncores 8 \
#     --irma /usr/local/bin/IRMA \
#     --pattern "^%RUN%.*\.(fastq|fastq.gz)$" \
#     --overwrite
# You may also provide runs via a file (one run ID per line):
#   Rscript irma_pipeline.R --project /path --runs_file runs.txt

# Load Dependencies ---------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(renv))
source("util.R")


# Parse flags ----------------------------------------------------------------
option_list <- list(
  
  # Required inputs
  make_option(c("-f", "--fastq"), type = "character", default = NULL, 
              help = "directory containing fastq input", metavar = "DIR"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Working directory", metavar = "FILE"),
  make_option(c("-r", "--runs"), type = "character", default = NULL,
              help = "Comma-separated run IDs, e.g. PJ6FV5,BZPB3P", metavar = "RUNS"),
  
  # Optional parameter
  make_option(c("-m", "--min_cov"), type = "integer", default = 5,
              help = "Minimum coverage threshold for masking. [default %default]", metavar = "INT"),
  make_option(c("--multicore"), action = "store_true", default = FALSE,
              help = "Enable multicore processing (default)."),

  make_option("--irma", type = "character", default = "IRMA",
              help = "(Only for test) Path to IRMA executable ", metavar = "FILE"),

  make_option(c("--pattern"), type = "character", default = NULL,
              help = "Regex for FASTQ names; use %RUN% placeholder for run ID. [default %default]", metavar = "REGEX"),
  make_option(c("--dry_run"), action = "store_true", default = FALSE,
              help = "Plan only; list what would run, but do not execute IRMA."),
  make_option(c("--overwrite"), action = "store_true", default = FALSE,
              help = "Re-run samples even if outputs already exist.")

)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (Sys.which("IRMA") == "") stop("IRMA not found on PATH")
system2("IRMA", "--version")

message(">>> Checking inputs...")

# Checking input flags --------------------------------------------------
if (is.null(opt$f)) {
  stop("Missing required argument: --fastq")
}
if (is.null(opt$runs)) {
  stop("Missing required argument: --runs")
}

user_directory <- normalizePath(opt$output)
message("Working directory set to ", user_directory)
if (!dir.exists(user_directory)) {
  stop("Working directory does not exist: ", user_directory)
}

runs <- character(0)
if (!is.null(opt$runs)) {
  runs <- unlist(strsplit(opt$runs, "[,\\s]+"))
  runs <- unique(trimws(runs[nzchar(runs)]))
}
if (length(runs) == 0) {
  stop("No runs provided. Use --runs PJ6FV5,BZPB3P")
}

min_cov      <- as.integer(opt$min_cov)
use_multicore<- isTRUE(opt$multicore)
requested_nc <- as.integer(opt$ncores)
fastq_pat_tmpl <- opt$pattern

# Ensure 'IRMA_output' and 'Dino_Fastq' exist under user_directory
irma_output <- file.path(user_directory, "IRMA_output")
if (!dir.exists(irma_output)) {
  dir.create(irma_output, recursive = TRUE)
  message("Created directory: ", normalizePath(irma_output))
}

# Verify FASTQ directory structure 
# Ensure 'Dino_Fastq' exists and has subfolders 'Dino_<RunID>' for each run.
fastq_base <- normalizePath(opt$fastq)
message("Using fastq under ", fastq_base)

# Check/create subfolders for each run
for (run_id in runs) {
  subdir <- file.path(fastq_base, paste0("Dino_", run_id))
  if (!dir.exists(subdir)) {
    stop("Run folder does not exist under ", fastq_base, "for ", subdir)
  }
}

# 5) Define IRMA executable path -------------------------------------------------
irma_exec <- trimws(opt$irma)


# Create a directory for session information and logs
session_dir <- file.path(user_directory, "session_info")
if (!dir.exists(session_dir)) {
  dir.create(session_dir, recursive = TRUE)
  message("Created session info directory: ", session_dir)
}

setwd(user_directory)
message(">>> All clear. Starting pipeline...")

# 7) Run pipeline across multiple runs -------------------------------------------

# 7.1) Loop over specified runs
run_summaries <- vector("list", length(runs))
r <- 1
for (r in seq_along(runs)) {
  run_id <- runs[r]
  fastq_dir <- file.path(
    fastq_base,
    paste0("Dino_", run_id)
  )

  # Create a run-specific log file in session_dir
  log_file <- file.path(session_dir, paste0("missing_segments_", run_id, ".log"))
  file.create(log_file)  # Overwrite any old file
  
  message("searching ", fastq_dir)
  
  run_summary <- tryCatch(
    process_run(run_id, min_cov, fastq_dir, log_file = log_file,  use_parallel = opt$multicore, irma_exec = irma_exec),
    error = function(e) {
      warning(sprintf("IRMA error for run %s: %s", run_id, e$message))
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
    write.csv(run_summary, file = out_file, row.names = FALSE, quote = F)
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


message("Finished")