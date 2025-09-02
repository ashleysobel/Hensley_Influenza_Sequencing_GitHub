# Hensley_Influenza_Sequencing_GitHub

## Influenza Sequencing Pipeline
This is the Rscript version of the pipeline, with a docker image stabled with IRMA version 1.2.0.

## Installation
Pull the docker image from the repository
`docker pull tianyuclu/hensley_lab_flu_pipeline:1.0` \
Test if installation is successful \
`docker run --rm --entrypoint /opt/app/test.sh tianyuclu/hensley_lab_flu_pipeline:1.0`

## Usage 
Specify the input & output directory. Make sure that input directory contains sub-paths in format of Dino_$RUNID$. \
`INPUT_DIR="project/Dino_Fastq" ` \
`OUTPUT_DIR="project/output" ` \
Create your output dir if needed. `mkdir -p "$OUTPUT_DIR"` \

Run the docker container \
```
docker run --rm \
  -v "$INPUT_DIR:/in" \
  -v "$OUTPUT_DIR:/out" \
  tianyuclu/hensley_lab_flu_pipeline:1.0 \
  -o /out \
  -f /in \
  --runs "BZPB3P" \
  --min_cov 5`

```

## Parameters
To list all flags \
`docker run tianyuclu/hensley_lab_flu_pipeline:1.0 -h `
The flags for detailed configuration is still underconstruction. Please stay with input & output for the moment.
