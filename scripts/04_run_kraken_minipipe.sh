module load nextflow
nextflow run scripts/kraken2/main.nf \
-params-file params.json \
-c nextflow.config \
-profile farm22