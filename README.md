# fur_pathogen_identification

## Summary
This repository contains code that was used for, and results from, an analysis of viral- and bacterial-derived reads in samples from the the FUR Felis Catus project. 

The project was setup to run across the 13 FUR cat cohorts with `01_setup_project.R`. This created files in each `analysis/{STUDY_DIR}` directory that were then used to parameterise a pathogen identification nextflow pipeline (`scripts/main.nf`). This pipeline was run on each cohort (`scripts/02_run_kraken_minipipe.sh`) and the results were then collated across cohorts and visualised (`03_analyse_kraken2_results.R`).

Within the nextflow pipeline (`scripts/main.nf`), sample bam files were filtered to exclude reads that mapped to the host genome (Felis catus v9) with samtools (v1.19) and converted to fastq format. Kraken2 (v2.1.2) was then run on the remaining paired reads using the PlusPFP reference database (https://benlangmead.github.io/aws-indexes/k2; version March 2023) and with a confidence thershold set at `0.1`. Finally, the kraken2 results were re-formatted into an mpa-style report with Krakentools (v1.2.4). Both Kraken2 and krakentools mpa-style reports are provided for each sample in `analysis/{STUDY_DIR}/kraken2` and `analysis/{STUDY_DIR}/krakentools`

In `03_analyse_kraken2_results.R` the kraken2 and kraken-mpa style results were collated and filtered to include only those samples which passed initial sequencing QC (`metadata/*samples_to_keep*` files). Next, the proportion of minimizers from a species in the Kraken reference database found in a sample was calculated for each species: Proportion = Distinct minimizers observed / Total clade–level minimizers. 

To assess the significance of a taxon discovery, the probability of finding the observed proportion of minimizers for that taxon was calculated from a normal distribution given ***n*** reads associated with that taxon and the total ***N*** reads evaluated by Kraken2 in that sample. We applied the Benjamini-Hochberg procedure to the tests to correct for multiple comparisons and results with adjusted p-value <0.05 were considered statistically significant. 

Plots reflecting the abundance of reads from bacterial and viral taxa and statistical significance of detecting as many minimizers as were observed in each sample are provided in `analysis/{STUDY_DIR}/sparki/`

## Installation and Running 
- Re-running the initial steps of this analysis requires access to sample-level bams and to a high-performance computing (HPC) cluster running Ubuntu 20.04 (Focal Fossa). 
- Run time for the (`scripts/02_run_kraken_minipipe.sh`) stage is approximately two days wall-time. 
- The nextflow pipeline is configured to run via LSF but can be adapted to run on other schedulers through changes to the `nextflow.config` files.
- Regeneration of plots and tables in this repository from kraken2 outputs can be accomplished within the bundled project development container which has docker and R (4.4.0) installed (`.devcontainer/devcontainer.json`).

## Dependencies

- Nextflow (v22.10.4)
- samtools (v1.19)
- Kraken2 (v2.1.2)
- Krakentools (v1.2.4)
- R (v4.4.0)

All R dependencies are managed with renv (v1.0.9) and are listed in `renv.lock`. To install the dependencies, run the following command in R:

```R
renv::restore()
```

## Project organisation
```
.
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
├── metadata
│   ├── 6555_2711.samples_to_keep.nucleotide_variants.txt
│   ├── 6711_2820.samples_to_keep.nucleotide_variants.txt
│   ├── 6712_2822.samples_to_keep.nucleotide_variants.txt
│   ├── 6713_2821.samples_to_keep.nucleotide_variants.txt
│   ├── 6841_2964.samples_to_keep.nucleotide_variants.txt
│   ├── 6864_2965.samples_to_keep.nucleotide_variants.txt
│   ├── 6945_3142.samples_to_keep.nucleotide_variants.txt
│   ├── 6973_2987.samples_to_keep.nucleotide_variants.txt
│   ├── 6982_3135.samples_to_keep.nucleotide_variants.txt
│   ├── 6990_3065.samples_to_keep.nucleotide_variants.txt
│   ├── 7040_3064.samples_to_keep.nucleotide_variants.txt
│   ├── 7097_3073.samples_to_keep.nucleotide_variants.txt
│   ├── 7098_3140.samples_to_keep.nucleotide_variants.txt
│   ├── kraken2_complete_non_capped_may2023_inspect.txt
│   └── study_manifest.tsv
├── renv.lock
└── scripts
    ├── 01_setup_project.R
    ├── 02_run_kraken_minipipe.sh
    ├── 03_analyse_kraken2_results.R
    ├── constants.R
    ├── main.nf
    └── nextflow.config
```

## Contact 
- **Author**:  Jamie Billington (jb63@sanger.ac.uk)

## References
- [Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)](https://www.nature.com/articles/s41596-022-00738-y)
- [Twelve years of SAMtools and BCFtools Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008,](https://doi.org/10.1093/gigascience/giab008)