library(tidyverse)
library(glue)

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