# fur_pathogen_identification
A set of analysis scripts, metadata and results for running Germline variant calling on FUR cats.
Working Project dir: `/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification`
-  This project was setup with `scripts/01_setup_project.R`. Sample lists for each cohort were prepared from the entire set of FUR bam files and filtered post-hoc after running kraken
- Each cohort was processed using `https://gitlab.internal.sanger.ac.uk/DERMATLAS/analysis-methods/kraken2/-/tree/develop/scripts?ref_type=heads` version 0.4.1. Cohort analyses were parameterised using a
`<COHORT_DIR>/params.json` and `<COHORT_DIR>/nextflow.config`file and submitted via
`scripts/02_run_kraken_minipipe.sh`

- Visualisation, stats and filtering performed in `scripts/03_analyse_kraken2_results.R`