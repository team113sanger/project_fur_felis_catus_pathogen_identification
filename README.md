# fur_pathogen_identification

## Summary
This repository contains code that was used for and results from an analysis searching for viral- and bacterial-derived reads in samples from the the FUR Felis Catus project. 

In brief, the project was setup to run across the 13 FUR cat cohorts with (`01_setup_project.R`). This created files that were then used to parameterise a pathogen identification nextflow pipeline that was run on each cohort (`scripts/02_run_kraken_minipipe.sh`). The results were then collated and visualised in `03_analyse_kraken2_results.R`.

Cohort analyses were parameterised using a
`<COHORT_DIR>/params.json` and `<COHORT_DIR>/nextflow.config`file and submitted via `scripts/02_run_kraken_minipipe.sh`

- Visualisation, stats and filtering performed in `scripts/03_analyse_kraken2_results.R`


## Dependencies


## Organisation

├── README.md
├── analysis
│   ├── 6555_2711
│   ├── 6711_2820
│   ├── 6712_2822
│   ├── 6713_2821
│   ├── 6841_2964
│   ├── 6864_2965
│   ├── 6945_3142
│   ├── 6973_2987
│   ├── 6982_3135
│   ├── 6990_3065
│   ├── 7040_3064
│   ├── 7097_3073
│   ├── 7098_3140
│   ├── cross-cohort
│   ├── kraken2_all_sample_summary.tsv
│   └── kraken2_passed_summary.tsv
├── logs
├── renv.lock
├── scripts
│   ├── 01_setup_project.R
│   ├── 02_run_kraken_minipipe.sh
│   ├── 03_analyse_kraken2_results.R
│   ├── constants.R
│   ├── docker-compose.yml
│   ├── kraken2
│   ├── kraken2_all_sample_summary.tsv
│   ├── kraken2_passed_summary.tsv
│   └── nextflow.config
└── study_manifest.tsv


## Contact 
- **Author**:  Jamie Billington (jb63@sanger.ac.uk)