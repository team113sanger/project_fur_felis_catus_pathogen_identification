library(tidyr)
library(dplyr)
library(readr)
library(fs)
library(stringr)
source("constants.R")
taxonomy <- c(
  "Domain", "Kingdom", "Phylum", "Class",
  "Order", "Family", "Genus", "Species"
)

# Use a filesystem helper to gather all files
krakentools_files <- fs::dir_ls("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6555_2711/results/krakentools", glob = "*.mpa")

krakentools_df <- read_tsv(krakentools_files, col_names = c("id", "reads"), id = "sample_id") |>
  separate(id, into = taxonomy, sep = "\\|") |> # Split rows by "|" into primrary taxa
  mutate(across(taxonomy, ~ str_remove(.x, pattern = "[a-z]__"))) |> # Cleanup the names in the taxonomy columns
  mutate(sample_id = str_remove(basename(sample_id), ".kraken.mpa")) |> # Simplify sample ids
  mutate(
    taxon_leaf = coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom, Domain),
    .before = "Domain"
  ) |> # Collect the rightmost non-NA item in each row
  mutate(taxon_leaf = str_replace_all(taxon_leaf, pattern = "_", replacement = " "))

kraken_files <- fs::dir_ls("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6555_2711/results/kraken2", glob = "*.kraken")

kraken2_df <- read_tsv(kraken_files,
  col_names = c("pct", "n_fragments_clade", "n_fragments_taxon", "tax.level", "n_distinct_minimisers", "rank", "ncbi", "name"),
  id = "sample_id"
) |>
  filter(!grepl(rank, pattern = "[0-9]")) |> # Remove the non-primary ranks
  mutate(sample_id = str_remove(basename(sample_id), "_0.kraken"))

ref_db <- read_tsv("/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_capped_march2023/inspect.txt",
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
  by = c("name" = "taxon_leaf", "sample_id")
) |>
  left_join(ref_db, by = c("name" = "taxon", "rank" = "rank"))

class_unclass_df <- combined_df |>
  filter(name %in% c("unclassified", "root")) |>
  rename(
    "type" = "name",
    "sample" = "sample_id",
    "n_reads" = "n_fragments_clade"
  )

read_totals <- class_unclass_df |>
  group_by(sample) |>
  summarise(total_read = sum(n_reads))

plot <- ggplot2::ggplot(
  class_unclass_df,
  ggplot2::aes(x = type, y = log10(n_reads))
) +
  ggplot2::geom_violin(scale = "width", fill = "white", color = "black") +
  ggplot2::geom_line(ggplot2::aes(group = sample), alpha = 0.25) +
  ggplot2::geom_point(ggplot2::aes(color = type), alpha = 0.5) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    # x-axis
    axis.text.x = ggplot2::element_text(size = 12, vjust = 0.5),
    axis.title.x = ggplot2::element_text(size = 14),
    # y-axis
    axis.text.y = ggplot2::element_text(size = 12),
    axis.title.y = ggplot2::element_text(size = 14, angle = 90),
    # legend
    legend.position = "none"
  ) +
  ggplot2::xlab("\nRead classification") +
  ggplot2::ylab(expression("log"[10] ~ "(# reads)")) +
  ggplot2::scale_color_manual(
    values = c("indianred2", "royalblue")
  )
png("/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6555_2711/results/sparki/classification.png")
plot
dev.off()

by_domain <- combined_df |>
  filter(rank == "D") |>
  select(n_fragments_clade, rank, Domain, sample_id)


 view <- ggplot2::ggplot(
            by_domain, 
            ggplot2::aes(x = Domain, y = log10(n_fragments_clade), fill = Domain)
        ) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(size = 0.5) +
        ggplot2::theme_bw() +
         ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(size = 15),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title.y = ggplot2::element_text(size = 16),
            # legend
            legend.position = "none"
        ) +
        ggplot2::ylab(expression("log"[10]~"(# classified reads)")) +
        ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) +
        ggplot2::facet_wrap(~Domain, scales = "free")

png("example_violin.png")
view
dev.off()

png("barplot.png")
x_lab <- "\nProportion of classified reads\n(all domains)"
filename <- "nReadsDomains_barplot_with_eukaryotes"
colours <- c("gold1", "royalblue", "snow4", "indianred2")
domain_barplot(by_domain)
dev.off()



png("barplot_minus_euk.png")
x_lab <- "\nProportion of classified reads\n(non-eukaryotes only)"
filename <- "nReadsDomains_barplot_without_eukaryotes"
colours <- c("gold1", "royalblue", "indianred2")
domain_barplot(by_domain |> filter(Domain != "Eukaryota"))
dev.off()

ratio_df <- filtered_df |>
  mutate(
    ratio_taxon = n_distinct_minimisers / n_minimisers_taxon,
    ratio_clade = n_distinct_minimisers / n_minimisers_clade
  )


stats_df <- ratio_df |>
  mutate(proportion = ratio_clade / sum(ref_db[["n_minimisers_clade"]])) |>
  left_join(read_totals, by = c("sample_id" = "sample")) |>
  mutate(
    expected_clade = proportion * total_read)
    std_clade = sqrt(proportion * total_read * (1 - total_read * proportion))
  ) |>
  mutate(pvalue = pnorm(
    q = n_distinct_minimisers,
    mean = expected_clade,
    sd = std_clade,
    lower.tail = FALSE
  )) |>
  group_by(rank) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  mutate(significance = ifelse(padj <= 0.05, "Significant", "Non-significant")) |> 
  ungroup()

stats_df |> select(name, n_distinct_minimisers, expected_clade, pvalue, padj, significance) |> print(n =500)

stats_df |>
  dplyr::filter(Domain == "Viruses") |>
  tidyr::complete(sample_id, name, fill = list(significance = "Non-significant", Domain = "Viruses")) |>
  dplyr::select(sample_id, Domain, name, rank, ratio_clade, padj, significance, n_fragments_clade)



plot <- ggplot2::ggplot(
  stats_df,
  ggplot2::aes(x = sample_id, y = name, fill = ratio_clade, color = significance, size = log10(n_fragments_clade))
) +
  ggplot2::geom_point(shape = 21, stroke = 1.25) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(rank),
    scales = "free_y",
    space = "free_y"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 15),
    legend.text.align = 0,
    legend.title = ggplot2::element_text(size = 15),
    axis.title = ggplot2::element_text(size = 15),
    plot.title = ggplot2::element_text(size = 16.5, face = "bold", hjust = 0.5),
    strip.text.y = ggplot2::element_text(size = 15, face = "bold")
  ) +
  ggplot2::scale_fill_gradientn(
    colors = colorRampPalette((RColorBrewer::brewer.pal(9, "Reds")))(100),
    name = "Ratio between unique minimisers found in sample\nand total clade-level minimisers in database",
    limits = c(0, 1)
  ) +
  ggplot2::scale_color_manual(
    values = c("snow3", "black"),
    name = "Significance",
    labels = c("Non-significant", expression("Adjusted p-value" <= "0.05"))
  ) +
  ggplot2::scale_size_continuous(name = expression("log"[10] ~ "(clade-level reads)")) +
  ggplot2::xlab("Sample") +
  ggplot2::ylab("Taxon")
png("minimisers.png", width = 1000, height = 2500)
plot
dev.off()





 plot <- ggplot2::ggplot(
            report, 
            ggplot2::aes(x = colname_taxon, y = log10(colname_n_frag_clade), fill = colname_taxon)
        ) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(size = 0.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            # x-axis
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(size = 15),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title.y = ggplot2::element_text(size = 16),
            # legend
            legend.position = "none"
        ) +
        ggplot2::ylab(expression("log"[10]~"(# classified reads)")) +
        ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) +
        ggplot2::facet_wrap(~colname_taxon, scales = "free")



source("R/utilities.R")
source("R/helper.R")
source("R/plotting.R")
source("R/constants.R")



mpa_reports <- load_MPAreports("tests/testdata/krakentools", verbose = TRUE)
std_reports <- load_STDreports("tests/testdata/kraken2", verbose = TRUE)

ref_db <- loadReference("tests/testdata/inspect.txt")

mpa_reports <- mpa_reports |>
  addRank(verbose = TRUE) |>
  addConciseTaxon(verbose = TRUE) |>
  transfer_ncbiID(report_mpa = _, std_reports)


std_reports <- transferDomains(std_reports,
  mpa_reports,
  verbose = TRUE
)

std_reports <- std_reports |>
  add_nReads() |>
  add_DBinfo(ref_db)

mpa_reports <- transfer_nReads(
  mpa_reports,
  std_reports
) |>
  add_DBinfo(ref_db)

plotClassificationSummary_violin(
  combined_df,
  return_plot = TRUE,
  outdir = "tests/testdata",
  prefix = "Feline_lymphoma"
)

plotDomainReads_barplot(
  std_reports,
  include_eukaryotes = FALSE,
  include_sample_names = FALSE,
  orientation = "vertical",
  return_plot = TRUE,
  outdir = "tests/testdata",
  prefix = "Feline_lymphoma"
)


std_reports <- subset_STDreport(std_reports, include_human = FALSE)
std_reports <- assess_ratioMinimisers(std_reports)
std_reports <- assess_statSig(std_reports, ref_db)
std_reports

plotMinimisers_dotplot(
  std_reports,
  domain_arg = "Viruses",
  return_plot = TRUE,
  fig_width = 35,
  fig_height = 35,
  outdir = "tests/testdata",
  prefix = "Feline_lymphoma"
)



plotDomainReads_violin(
  std_reports,
  include_eukaryotes = TRUE,
  return_plot = TRUE,
  outdir = "outputs",
  prefix = "SebT"
)
