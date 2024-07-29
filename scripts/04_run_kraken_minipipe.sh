module load nextflow
nextflow run scripts/kraken2/main.nf \
-params-file /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6555_2711/params.json \
-c scripts/nextflow.config \
-profile farm22

nextflow run scripts/kraken2/main.nf \
-params-file /lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6982_3135/params.json \
-c scripts/nextflow.config \
-profile farm22