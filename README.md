# Hensley_Influenza_Sequencing_GitHub

## Influenza Sequencing Pipeline
This is the Rscript version of the pipeline, with a docker image stabled with IRMA version 1.2.0.

## Installation
Pull the docker image from the repository
`docker pull tianyuclu/hensley_lab_flu_pipeline:1.0`
Test if installation is successful
`docker run --rm --entrypoint /opt/app/test.sh tianyuclu/hensley_lab_flu_pipeline:1.0`

## Usage 
Specify the local mounted dir
`INPUT_DIR="$PWD/Dino_Fastq"
OUTPUT_DIR="$PWD/output" `
Run the docker container 
`docker run --rm \
  -v "$PWD:/work" -w /work \
  fluirma \
  -o "/work/$(basename "$OUTPUT_DIR")" \
  -f "/work/$(basename "$INPUT_DIR")" \
  --runs "BZPB3P" \
  --min_cov 5`

