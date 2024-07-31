module load nextflow
nextflow run scripts/kraken2/main.nf \
-resume \
-params-file analysis/7097_3073/params.json \
-c scripts/kraken2/nextflow.config \
-profile farm22
