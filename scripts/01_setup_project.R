library(tidyverse)
library(glue)
library(jsonlite)

project_dir <- "/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis"
studies <- read_tsv("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/study_manifest.tsv", 
                    )
studies 

output_dirs <- studies |> 
rowwise() |> 
mutate(path = glue("{project_dir}/{SeqScape_ID}_{CanApps_ID}")) |> 
pull(path)
walk(output_dirs, dir.create)
walk(paste0(output_dirs, "/metadata"), dir.create)
walk(paste0(output_dirs, "/fastqs"), dir.create)
walk(paste0(output_dirs, "/results"), dir.create)
walk(paste0(output_dirs, "/tmp"), dir.create)
walk(paste0(output_dirs, "/scripts"), dir.create)


file_ids <- paste0(output_dirs, "/params.json")
parameters <- map(output_dirs, 
    ~toJSON(
    list("fastq_files" = paste0(.x, "/fastqs/**{R1,R2}.fastq.gz"),
    "reference_db" = "/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_capped_march2023",
    "c_score" = 0.1,
    "outdir" = paste0(.x, "/results")),
    auto_unbox = TRUE, 
    pretty = TRUE))


map2(parameters, file_ids, ~write_file(x = .x, file = .y))

# read_tsv("/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_capped_march2023/inspect.txt", comment = "#",
# col_names = (c("pct", "n_reads", "n_minimisers", "rank", "taxid", "name"))) |> 
# filter(rank == "S") |>
# select(name)