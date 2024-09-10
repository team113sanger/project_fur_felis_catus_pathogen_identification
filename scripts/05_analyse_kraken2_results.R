################################################################################
# Setup
################################################################################

library(tidyr)
library(dplyr)
library(readr)
library(fs)
library(stringr)
library(glue)
library(purrr)

source("constants.R")


################################################################################
# Declare inputs
################################################################################


taxonomy <- c(
  "Domain", "Kingdom", "Phylum", "Class",
  "Order", "Family", "Genus", "Species"
)

project_dir <- "/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis"
ref_db_file <- "/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_non_capped_may2023/inspect.txt"



# Use a filesystem helper to gather all files
krakentools_files <- fs::dir_ls(
  path = project_dir,
  type = "file",
  glob = "*.mpa$",
  recurse = TRUE
)
krakentools_files <- krakentools_files[!grepl(krakentools_files, pattern = "work")]


kraken_files <- fs::dir_ls(project_dir,
  recurse = TRUE,
  glob = "*0.1.kraken"
)

kraken_files <- kraken_files[!grepl(kraken_files, pattern = "work")]


################################################################################
# Read in and merge
################################################################################


krakentools_df <- purrr::map_dfr(
  krakentools_files,
  ~ safe_read_tsv(.x,
    col_names = c("id", "reads"),
    id = "file_id"
  )
) |>
  separate(id, into = taxonomy, sep = "\\|") |> # Split rows by "|" into primary taxa
  mutate(across(
    taxonomy,
    ~ str_remove(.x, pattern = "[a-z]__")
  )) |> # Cleanup the names in the taxonomy columns
  mutate(
    sample_id = str_remove(basename(file_id), ".kraken.mpa"),
    .after = file_id
  ) |> # Simplify sample ids
  mutate(cohort = str_match(file_id,
    pattern = ".*/analysis/(.*)/results/krakentools/.*"
  )[, 2], .after = sample_id) |> # Simplify sample ids
  mutate(
    taxon_leaf =
      coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom, Domain),
    .before = "Domain"
  ) |> # Collect the rightmost non-NA item in each row
  mutate(
    taxon_leaf =
      str_replace_all(taxon_leaf, pattern = "_", replacement = " ")
  )




kraken2_df <- read_tsv(kraken_files,
  col_names = c(
    "pct",
    "n_fragments_clade",
    "n_fragments_taxon",
    "tax.level",
    "n_distinct_minimisers",
    "rank",
    "ncbi",
    "name"
  ),
  id = "file_id"
) |>
  filter(!grepl(rank, pattern = "[0-9]")) |> # Remove the non-primary ranks
  mutate(sample_id = str_remove(basename(file_id), "_0.1.kraken"), .after = file_id) |>
  mutate(cohort = str_match(file_id, pattern = ".*/analysis/(.*)/results/kraken2.*")[, 2], .after = sample_id)



ref_db <- read_tsv(ref_db_file,
  comment = "#",
  c(
    COLNAME_REF_DB_PCT_FRAG_CLADE,
    COLNAME_REF_DB_MINIMISERS_CLADE,
    COLNAME_REF_DB_MINIMISERS_TAXON,
    COLNAME_REF_DB_RANK,
    COLNAME_REF_DB_NCBI_ID,
    COLNAME_REF_DB_TAXON
  )
) |>
  filter(!grepl(rank, pattern = "[0-9]")) # Remove the non-primary ranks


combined_df <- left_join(kraken2_df, krakentools_df,
  by = c("name" = "taxon_leaf", "sample_id", "cohort")
) |>
  left_join(ref_db, by = c("name" = "taxon", "rank" = "rank"))


exclude_files <- fs::dir_ls(
  path = "/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/fur_cnvkit/metadata/exclude_files",
  type = "file",
  glob = "*.txt$",
  recurse = TRUE
)

excluded_samples <- read_tsv(exclude_files, col_names = c("sample_id"), id = "cohort") |> 
mutate(cohort = str_remove(basename(cohort), pattern = ".samples_to_exclude.txt"))
combined_df 

qc_passed_samples <- combined_df |> anti_join(excluded_samples, by = c("sample_id", "cohort"))


write_tsv(qc_passed_samples, "kraken2_passed_summary.tsv")
write_tsv(combined_df, "kraken2_all_sample_summary.tsv")


cohort_analysis <- function(x, cohort, project_dir, ref_db) {
  class_unclass_df <- x |>
    dplyr::filter(name %in% c("unclassified", "root")) |>
    dplyr::rename(
      "type" = "name",
      "sample" = "sample_id",
      "n_reads" = "n_fragments_clade"
    )

  read_totals <- class_unclass_df |>
    dplyr::group_by(sample) |>
    dplyr::summarise(total_read = sum(n_reads))

  png(glue::glue("{project_dir}/{cohort}/results/sparki/classification.png"))
  classification_plot <- ggplot2::ggplot(
    class_unclass_df,
    ggplot2::aes(x = type, y = log10(n_reads))
  ) +
    ggplot2::geom_violin(scale = "width", fill = "white", color = "black") +
    ggplot2::geom_line(ggplot2::aes(group = sample), alpha = 0.25) +
    ggplot2::geom_point(ggplot2::aes(color = type), alpha = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
      axis.title.x = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 12),
      axis.title.y = ggplot2::element_text(size = 14, angle = 90),
      legend.position = "none"
    ) +
    ggplot2::xlab("\nRead classification") +
    ggplot2::ylab(expression("log"[10] ~ "(# reads)")) +
    ggplot2::scale_color_manual(
      values = c("indianred2", "royalblue")
    )
  print(classification_plot)
  dev.off()

  by_domain <- x |>
    dplyr::filter(rank == "D") |>
    dplyr::select(c(n_fragments_clade, rank, Domain, sample_id))

  png(glue::glue("{project_dir}/{cohort}/results/sparki/per_domain_reads.png"))
  per_domain_reads_plot <- ggplot2::ggplot(
    by_domain,
    ggplot2::aes(x = Domain, y = log10(n_fragments_clade), fill = Domain)
  ) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::geom_jitter(size = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.title.y = ggplot2::element_text(size = 16),
      legend.position = "none"
    ) +
    ggplot2::ylab(expression("log"[10] ~ "(# classified reads)")) +
    ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2"))
  print(per_domain_reads_plot)
  dev.off()

  png(glue::glue("{project_dir}/{cohort}/results/sparki/per_domain_proportions.png"))
  x_lab <- "\nProportion of classified reads\n(all domains)"
  filename <- "nReadsDomains_barplot_with_eukaryotes"
  colours <- c("royalblue", "snow4", "indianred2", "gold")
  domain_barplot(by_domain, x_lab = x_lab, colours)
  dev.off()

  png(glue::glue("{project_dir}/{cohort}/results/sparki/per_domain_proportions_no_euks.png"))
  x_lab <- "\nProportion of classified reads\n(non-eukaryotes only)"
  filename <- "nReadsDomains_barplot_without_eukaryotes"
  colours <- c("royalblue", "indianred2", "gold")
  domain_barplot(by_domain |> dplyr::filter(Domain != "Eukaryota"), x_lab = x_lab, colours)
  dev.off()

  filtered_df <- x |> dplyr::filter(rank %in% c("S"))

  ratio_df <- filtered_df |>
    dplyr::mutate(
      ratio_taxon = n_distinct_minimisers / n_minimisers_taxon,
      ratio_clade = n_distinct_minimisers / n_minimisers_clade
    )

  total_minimisers <- ref_db |>
    dplyr::filter(rank %in% c("S")) |>
    dplyr::summarise(single_counts = sum(n_minimisers_clade)) |>
    dplyr::pull(single_counts)

  stats_df <- ratio_df |>
    dplyr::mutate(proportion = n_minimisers_clade / total_minimisers) |>
    dplyr::left_join(read_totals, by = c("sample_id" = "sample")) |>
    dplyr::mutate(
      expected_clade = proportion * total_read,
      std_clade = sqrt((proportion * total_read) * (1 - proportion))
    ) |>
    dplyr::mutate(pvalue = stats::pnorm(
      q = n_distinct_minimisers,
      mean = expected_clade,
      sd = std_clade,
      lower.tail = FALSE
    )) |>
    dplyr::group_by(rank) |>
    dplyr::mutate(padj = stats::p.adjust(pvalue, method = "BH")) |>
    dplyr::mutate(significance = ifelse(padj <= 0.05, "Significant", "Non-significant")) |>
    dplyr::ungroup()

  q_domain <- "Bacteria"
  bacteria <- stats_df |>
    dplyr::filter(rank == "S") |>
    dplyr::filter(Domain == q_domain) |>
    dplyr::group_by(name) |>
    dplyr::filter(n() > 3) |>
    dplyr::ungroup() |>
    dplyr::filter(n_fragments_clade > 5) |>
    tidyr::complete(sample_id, name,
      fill = list(
        significance = "Non-significant",
        rank = "S",
        Domain = q_domain
      )
    ) |>
    dplyr::select(
      sample_id, Domain, name, rank,
      ratio_clade, padj, significance, n_fragments_clade
    )

  print("qdom")
  q_domain <- "Viruses"
  viruses <- stats_df |>
    dplyr::filter(rank %in% c("S")) |>
    dplyr::filter(Domain == q_domain) |>
    dplyr::filter(n_fragments_clade > 5) |>
     tidyr::complete(sample_id, name, fill = list(
      significance = "Non-significant",
      rank = "S",
      Domain = q_domain
    )) |>
    dplyr::select(
      sample_id, Domain, name, rank, ratio_clade, padj,
      significance, n_fragments_clade
    )

  write_tsv(stats_df, glue("{project_dir}/{cohort}/results/sparki/significance_table.tsv"))
  bacterial_plot <- plot_minimisers(bacteria)

  png(glue::glue("{project_dir}/{cohort}/results/sparki/bacterial_minimisers.png"),
    width = 1500, height = 2600
  )
  print(bacterial_plot)
  dev.off()

  if (nrow(viruses) >1){
  viral_plot <- plot_minimisers(viruses)
  png(glue::glue("{project_dir}/{cohort}/results/sparki/viral_minimisers.png"),
    width = 1500, height = 500
  )
  print(viral_plot) # Ensure the plot is printed
  dev.off()
  }

}



# cohort_dfs <- combined_df |> split(f = combined_df[["cohort"]])

# imap(cohort_dfs, ~ cohort_analysis(x = .x, cohort = .y, project_dir, ref_db))


passed_dfs <- qc_passed_samples |> split(f = qc_passed_samples[["cohort"]])

imap(passed_dfs, ~ cohort_analysis(x = .x, cohort = .y, project_dir, ref_db))





cat_lists <- tibble(
  sample_id = basename(krakentools_files),
  cohort = str_match(dirname(krakentools_files), pattern = ".*analysis/(.*)/results")[, 2]
) |>
  mutate(sample_id = str_remove(sample_id, pattern = ".kraken.mpa"))

cohort_lists <- cat_lists |>
  group_by(cohort) |>
  group_split() |>
  set_names(unique(cat_lists$cohort))

imap(cohort_lists, ~write_tsv(.x, glue("{project_dir}/{.y}/sample_list.tsv")))
