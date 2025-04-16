#!/bin/sh
module load nextflow
PROJECT_DIR=/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification
cd $PROJECT_DIR

# Define sample IDs as an array
SAMPLES=(
    "6711_2820"
    "6712_2822"
    "6713_2821"
    "6841_2964"
    "7097_3073"
    "6864_2965"
    "6945_3142"
    "6973_2987"
    "6982_3135"
    "6990_3065"
    "7040_3064"
    "7098_3140"
)

# Loop through each sample and submit job
for SAMPLE in "${SAMPLES[@]}"; do
    bsub -q oversubscribed -M 2000 -n 1 -R'select[mem>2000] rusage[mem=2000] span[hosts=1]' \
         -o "logs/${SAMPLE}_%J.o" -e "logs/${SAMPLE}_%J.err" \
         "nextflow run scripts/kraken2/scripts/main.nf -params-file analysis/${SAMPLE}/params.json -c analysis/${SAMPLE}/nextflow.config -profile farm22"
    # Optional: Add a small delay between submissions
    # sleep 1
done