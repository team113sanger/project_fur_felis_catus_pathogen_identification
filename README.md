# fur_pathogen_identification
A set of analysis scripts, metadata and results for running Germline variant calling on FUR cats.
-  This project was setup with `scripts/01_setup_project.R`. Sample lists for each cohort were prepared from:
`/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cat_cohort_files`
- Each cohort was processed using `https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/kraken2/-/tree/develop/scripts?ref_type=heads` version 0.4.1. Cohort analyses were parameterised using a
`<COHORT_DIR>/params.json` and `<COHORT_DIR>/nextflow.config`file. 
