# README.md
## Hensley_Influenza_Sequencing_GitHub
### Influenza Sequencing Pipeline

**An R-based workflow for processing influenza sequencing data with IRMA.**

> **README.md** placed at your repository root—GitHub will render it automatically.

---

## 🔧 Prerequisites

Before running the pipeline, make sure you have:

1. **R** (version ≥ 3.6)
2. **IRMA** (Iterative Refinement Meta‑Assembler)
   - Download or clone from the official source.
   - Make the script executable, for example:
     ```bash
     chmod +x /full/path/to/IRMA
     ```
3. **R packages**: `parallel`, `renv` (the script will install these if missing).

---

## 📂 Repository Layout

```
├── README.md                      # This file (Markdown)
├── project_path.txt               # Contains your project root path
├── irma_path.txt                  # Contains full path to IRMA executable
├── pipeline.R                     # Main R script
├── Dino_Fastq/                    # Input FASTQ files
│   └── Dino_<RunID>/              # e.g. Dino_XLHBTH
├── IRMA_output/                   # Pipeline outputs
│   ├── <RunID>/                   # Raw IRMA per‐sample folders
│   ├── Masked_consensus/          # Masked FASTAs
│   │   ├── HA/                    # HA-only masked FASTAs
│   │   ├── NA/                    # NA-only masked FASTAs
│   │   └── <sample>_masked.fasta  # All‐segments masked FASTAs
│   └── Unmasked_consensus/        # Unmasked FASTAs
│       └── <sample>_unmasked.fasta
├── <RunID>_coverage_summary.csv   # Coverage summary tables
├── session_info/                  # Logs and session metadata
│   ├── missing_segments_<RunID>.log
│   ├── sessionInfo.txt
│   └── installed_packages.tsv
└── pipeline_processing_citations.bib
```

---

## Setup Guide

1. **Clone the repo**
   ```bash
   git clone https://github.com/yourusername/your-repo.git
   cd your-repo
   ```

2. **Create configuration files** in the root:
   - **project_path.txt** (one line):
     ```txt
     /full/path/to/your-repo
     ```
   - **irma_path.txt** (one line):
     ```txt
     /full/path/to/IRMA
     ```

3. **Add your FASTQ files** under `Dino_Fastq/Dino_<RunID>/`:
   ```bash
   mkdir -p Dino_Fastq/Dino_XLHBTH
   cp /path/to/*.fastq Dino_Fastq/Dino_XLHBTH/
   ```

4. **Install required R packages** (if prompted):
   ```r
   install.packages(c("parallel", "renv"))
   ```

---

## Configure Your Run

Edit the top of **pipeline.R** to set:

- `runs` &mdash; vector of one or more run IDs:
  ```r
  runs <- c("XLHBTH", "ABC123")
  ```
- `min_cov` &mdash; integer coverage threshold; bases with coverage `< min_cov` are masked to `N` (default `5`).
- `use_multicore` &mdash; `TRUE` to speed up masked FASTA generation; `FALSE` to run sequentially.

---

## How to Run

From a terminal in the project root:
```bash
Rscript pipeline.R
```

Or inside an R session:
```r
setwd(trimws(readLines("project_path.txt", 1)))
source("pipeline.R")
```

The script will:

1. Verify or create required directories.
2. Load IRMA and R packages.
3. Run IRMA on each FASTQ sample.
4. Generate masked and unmasked consensus FASTAs.
5. Produce separate HA/NA masked FASTAs.
6. Summarize coverage, logging any missing segments.

---

## Results & Outputs

- **All‐segments FASTAs**: `IRMA_output/Masked_consensus/` and `Unmasked_consensus/`
- **HA-only & NA-only FASTAs**: `IRMA_output/Masked_consensus/HA/` & `…/NA/`
- **Coverage table**: `<RunID>_coverage_summary.csv`
- **Logs**: `session_info/missing_segments_<RunID>.log`
- **Session info**: `sessionInfo.txt`, `installed_packages.tsv`

---

## Troubleshooting

- **No FASTQ files found**: Ensure FASTQs are named `<RunID>_… .fastq` in the correct subfolder.
- **Permission denied**: Confirm your IRMA script has execute permission (`chmod +x`).
- **Empty consensus**: Inspect `session_info/missing_segments_<RunID>.log` for missing segments.

---

_For further help or to report issues, please open a GitHub issue in this repository._

