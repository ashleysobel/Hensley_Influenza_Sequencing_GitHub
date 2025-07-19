# README.md
## Hensley_Influenza_Sequencing_GitHub
### Influenza Sequencing Pipeline

**An R-based workflow for processing influenza sequencing data with IRMA.**

> **README.md**

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
├── Hensley_Influenza_Sequencing_GitHub.R                     # Main R script
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
1. **Clone the repository**
  
  * **Terminal**:
  
```bash
git clone https://github.com/ashleysobel/Hensley_Influenza_Sequencing_GitHub.git
cd Hensley_Influenza_Sequencing_GitHub
```
* **RStudio (GUI)**:
  
  1. Open RStudio → **File → New Project… → Version Control → Git**
  2. Repository URL:
```
git@github.com:ashleysobel/Hensley_Influenza_Sequencing_GitHub.git
```
3. Choose or create a local folder → **Create Project**
  
  2. **Configure your paths**
  Create two files in the project root (no quotes, no extra lines):
  
  * **project\_path.txt**
  *Purpose*: tells the R script where your project directory is.
*Contents* (example):
  
  ```txt
/Users/youruser/code/Hensley_Influenza_Sequencing_GitHub
```
* **irma\_path.txt**
  *Purpose*: tells the R script where to find the IRMA executable.
*Contents* (example):
```txt
/Users/youruser/code/flu-amd/IRMA
```

You can edit these in any text editor (RStudio, TextEdit, VSCode, etc.).

3. **Add your FASTQ files**

Download your FASTQs (e.g. from Plasmidsaurus). Note the six-character **Run ID** (e.g. `PJ6FV5`).

**Terminal**  
```bash
mkdir -p Dino_Fastq/Dino_PJ6FV5
mv ~/Downloads/PJ6FV5_*.fastq Dino_Fastq/Dino_PJ6FV5/

     
## Configure Your Run

Open **Hensley_IRMA_Parallel_Processing.R** and edit these at the top:
  
```r
# 2.2) Run IDs (six-character strings)
runs <- c("PJ6FV5", "BZPB3P", "XLHBTH")

# 2.3) Coverage threshold (positions with coverage < min_cov mask to "N")
min_cov <- 5

# 2.4) Parallel processing?
use_multicore <- TRUE  # TRUE = faster, FALSE = single-threaded
```

## How to Run

This repo is intended to run through RStudio. A version of this repo intended to run Docker is in process. Please check back later.

### Automatic package installation

The script installs any missing R packages (`parallel`, `renv`, …) automatically.

### Via RStudio (GUI)

1. Open the project in RStudio.
2. Open **Hensley_Influenza_Sequencing_GitHub.R**.
3. Click **Source** (or press Cmd+Shift+S).

### Via command line

From the project root:
  
```bash
Rscript Hensley_Influenza_Sequencing_GitHub.R
```

Or in an interactive R session:
  
```r
setwd(trimws(readLines("project_path.txt", 1)))
source("Hensley_Influenza_Sequencing_GitHub.R")
```


**What it does**:
  
1. Creates/verifies directories: `IRMA_output`, `Dino_Fastq`, `session_info`.
2. Installs & loads required packages.
3. Calls IRMA on each FASTQ.
4. Generates masked/unmasked consensus FASTAs (with HA/NA splits).
5. Writes `<RunID>_coverage_summary.csv` and log files.
---

## Results & Outputs

- **All‐segments FASTAs**: `IRMA_output/Masked_consensus/` and `Unmasked_consensus/`
- **HA-only & NA-only FASTAs**: `IRMA_output/Masked_consensus/HA/` & `…/NA/`
- **Coverage table**: `<RunID>_coverage_summary.csv`
- **Logs**: `session_info/missing_segments_<RunID>.log`
- **Session info**: `sessionInfo.txt`, `installed_packages.tsv`

---

## Troubleshooting

- **IRMA install errors**: On macOS arm64 (M1/M4), installing IRMA via Conda can fail due to missing `blat`. Download the CDC v1.2.0 release ([https://github.com/CDCgov/irma/releases/tag/v1.2.0](https://github.com/CDCgov/irma/releases/tag/v1.2.0)) instead, which bundles all required third-party tools.
- **Gatekeeper pop‑ups**: macOS may block IRMA’s executables (e.g. `blat`, `pigz`) with “Not Opened” warnings. Avoid this by temporarily disabling Gatekeeper:

```bash
  sudo spctl --master-disable
```

  Or in **System Settings → Privacy & Security → Allow apps downloaded from → Anywhere**. Re-enable with:

```bash
  sudo spctl --master-enable
```
- **Per‑binary exceptions**: If you prefer not to disable Gatekeeper globally, right‑click the blocked binary in Finder, select **Open**, then confirm. This creates a persistent exception for that file.
- **No FASTQ files found**: Ensure FASTQs are named `<RunID>_… .fastq` in the correct subfolder.
- **Permission denied**: Confirm your IRMA script has execute permission (`chmod +x`).
- **Empty consensus**: Inspect `session_info/missing_segments_<RunID>.log` for missing segments.


---

_For further help or to report issues, please open a GitHub issue in this repository._

