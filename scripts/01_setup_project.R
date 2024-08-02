library(tidyverse)
library(glue)
library(jsonlite)
library(fs)

study_paths <- fs::dir_ls("/lustre/scratch125/casm/team113da/users/bf14/samples", 
           type = "dir") |> 
as_tibble() |> 
mutate(study_id = basename(value)) |> 
filter(grepl(study_id, pattern = "[0-9]+")) |> 
mutate(CanApps_ID = as.numeric(study_id))


project_dir <- "/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis"
studies <- read_tsv("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/study_manifest.tsv", 
                    )
studies 

study_info <- studies |> 
left_join(study_paths, by = c("CanApps_ID")) |>
rowwise() |> 
mutate(path = glue("{project_dir}/{SeqScape_ID}_{CanApps_ID}")) 

output_dirs <- study_info |> pull(path)
walk(output_dirs, dir.create)
walk(paste0(output_dirs, "/metadata"), dir.create)
walk(paste0(output_dirs, "/fastqs"), dir.create)
walk(paste0(output_dirs, "/results"), dir.create)
walk(paste0(output_dirs, "/tmp"), dir.create)
walk(paste0(output_dirs, "/scripts"), dir.create)


file_ids <- paste0(output_dirs, "/params.json")

parameters <- map2(output_dirs, as.list(study_info[["value"]]),
    ~toJSON(
    list("bamfiles" = paste0(.y, "*/*/*.bam"),
    "reference_db" = "/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_non_capped_may2023",
    "c_score" = 0.1,
    "outdir" = paste0(.x, "/results")),
    auto_unbox = TRUE, 
    pretty = TRUE))


map2(parameters, file_ids, ~write_file(x = .x, file = .y))

# read_tsv("/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_capped_march2023/inspect.txt", comment = "#",
# col_names = (c("pct", "n_reads", "n_minimisers", "rank", "taxid", "name"))) |> 
# filter(rank == "S") |>
# select(name)

map(output_dirs, ~file_copy("nextflow.config", paste0(.x, "/nextflow.config")))

