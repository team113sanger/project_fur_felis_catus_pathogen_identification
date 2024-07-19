source("R/utilities.R")
source("R/helper.R")
source("R/plotting.R")
source("R/constants.R")

mpa_reports <- load_MPAreports("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/6982_Feline_lymphoma/pathogen_identification/krakentools", verbose = TRUE)
std_reports <- load_STDreports("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/6982_Feline_lymphoma/pathogen_identification/kraken2", verbose = TRUE)

ref_db <- loadReference("/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_capped_march2023/inspect.txt")

mpa_reports <- mpa_reports |>
               addRank(verbose = TRUE) |>
               addConciseTaxon(verbose = TRUE) |>
               transfer_ncbiID(std_reports)

std_reports <- transferDomains(std_reports, 
                               mpa_reports, 
                               verbose = TRUE)

std_reports <- std_reports |> 
               add_nReads() |>
               add_DBinfo(ref_db)

mpa_reports <- transfer_nReads(mpa_reports, 
                               std_reports) |> 
                add_DBinfo(ref_db)

plotClassificationSummary_violin(
  std_reports, 
  include_eukaryotes = FALSE, 
  return_plot = TRUE,
  outdir = "test/outputs/",
  prefix = "Feline_lymphoma"
)

plotDomainReads_barplot(
  std_reports, 
  include_eukaryotes = FALSE, 
  include_sample_names = FALSE, 
  orientation = "vertical", 
  return_plot = TRUE,
  outdir = "test/outputs/",
  prefix = "Feline_lymphoma")


std_reports <- subset_STDreport(std_reports, include_human = FALSE)
std_reports <- assess_ratioMinimisers(std_reports)
std_reports <- assess_statSig(std_reports, ref_db)
