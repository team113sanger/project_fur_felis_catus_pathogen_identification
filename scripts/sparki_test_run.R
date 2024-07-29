library(tidyr)
library(dplyr)
library(readr)
library(fs)
library(stringr)
library(glue)

source("constants.R")

taxonomy <- c(
  "Domain", "Kingdom", "Phylum", "Class",
  "Order", "Family", "Genus", "Species"
)

project_dir <- "/lustre/scratch125/casm/team113da/projects/FUR/FUR_analysis/FUR_analysis_cat/pathogen_identification/analysis/6712_2822"
ref_db_file <- "/lustre/scratch124/casm/team113/ref/DERMATLAS/kraken2_complete_capped_march2023/inspect.txt"

dir.create(glue("{project_dir}/results/sparki"))

# Use a filesystem helper to gather all files
krakentools_files <- fs::dir_ls(glue("{project_dir}/results/krakentools", 
                      glob = "*.mpa"))

krakentools_df <- read_tsv(krakentools_files, 
                           col_names = c("id", "reads"), id = "sample_id") |>
                  separate(id, into = taxonomy, sep = "\\|") |> # Split rows by "|" into primrary taxa
                  mutate(across(taxonomy, ~ str_remove(.x, pattern = "[a-z]__"))) |> # Cleanup the names in the taxonomy columns
                  mutate(sample_id = str_remove(basename(sample_id), ".kraken.mpa")) |> # Simplify sample ids
                  mutate(
                    taxon_leaf = 
                    coalesce(Species, Genus, Family, Order, Class, Phylum, Kingdom, Domain),
                    .before = "Domain") |> # Collect the rightmost non-NA item in each row
                  mutate(taxon_leaf = 
                  str_replace_all(taxon_leaf, pattern = "_", replacement = " "))

kraken_files <- fs::dir_ls(glue("{project_dir}/results/kraken2"), 
                    glob = "*0.1.kraken")

kraken2_df <- read_tsv(kraken_files,
                       col_names = c("pct", 
                                     "n_fragments_clade", 
                                     "n_fragments_taxon", 
                                     "tax.level", 
                                     "n_distinct_minimisers", 
                                     "rank", 
                                     "ncbi", 
                                     "name"),
                                     id = "sample_id") |>
              filter(!grepl(rank, pattern = "[0-9]")) |> # Remove the non-primary ranks
              mutate(sample_id = str_remove(basename(sample_id), "_0.1.kraken"))

ref_db <- read_tsv(ref_db_file, comment = "#",
        c(
        COLNAME_REF_DB_PCT_FRAG_CLADE,
        COLNAME_REF_DB_MINIMISERS_CLADE,
        COLNAME_REF_DB_MINIMISERS_TAXON,
        COLNAME_REF_DB_RANK,
        COLNAME_REF_DB_NCBI_ID,
        COLNAME_REF_DB_TAXON)) |>
        filter(!grepl(rank, pattern = "[0-9]")) # Remove the non-primary ranks

print("Read in")
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
  filter(rank == "R") |>
  group_by(sample) |>
  summarise(total_read = sum(n_reads))

png(glue("{project_dir}/results/sparki/classification.png"))
ggplot2::ggplot(
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
    values = c("indianred2", "royalblue"))
dev.off()
by_domain <- combined_df |>
  filter(rank == "D") |> 
  select(c(n_fragments_clade, rank, Domain, sample_id))
print("check finish")

print("check next")
png(glue("{project_dir}/results/sparki/per_domain_reads.png"))
ggplot2::ggplot(
            by_domain, 
            ggplot2::aes(x = Domain, y = log10(n_fragments_clade), fill = Domain)
        ) +
        ggplot2::geom_violin(scale = "width") +
        ggplot2::geom_jitter(size = 0.5) +
        ggplot2::theme_bw() +
         ggplot2::theme(
            # x-axis
            # axis.text.x = ggplot2::element_blank(),
            # axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            strip.text.x = ggplot2::element_text(size = 15),
            # y-axis
            axis.text.y = ggplot2::element_text(size = 15),
            axis.title.y = ggplot2::element_text(size = 16),
            # legend
            legend.position = "none"
        ) +
        ggplot2::ylab(expression("log"[10]~"(# classified reads)")) +
        ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) 
dev.off()
print("pass from here")

png(glue("{project_dir}/results/sparki/per_domain_proportions.png"))
x_lab <- "\nProportion of classified reads\n(all domains)"
filename <- "nReadsDomains_barplot_with_eukaryotes"
colours <- c("royalblue","snow4", "indianred2", "gold")
domain_barplot(by_domain)
dev.off()



png(glue("{project_dir}/results/sparki/per_domain_proportions_no_euks.png"))
x_lab <- "\nProportion of classified reads\n(non-eukaryotes only)"
filename <- "nReadsDomains_barplot_without_eukaryotes"
colours <- c("royalblue", "indianred2", "gold")
domain_barplot(by_domain |> filter(Domain != "Eukaryota"))
dev.off()


filtered_df <- combined_df |> 
               filter(rank %in% c("S","G","F"))

ratio_df <- filtered_df |>
  mutate(
    ratio_taxon = n_distinct_minimisers / n_minimisers_taxon,
    ratio_clade = n_distinct_minimisers / n_minimisers_clade
  )


print("Sample")
total_minimisers <- ref_db |> 
filter(rank %in% c("F")) |>
summarise(single_counts = sum(n_minimisers_clade)) |> 
pull(single_counts)


stats_df <- ratio_df |>
  mutate(proportion = ratio_clade / total_minimisers) |>
  left_join(read_totals, by = c("sample_id" = "sample")) |>
  mutate(
    expected_clade = proportion * total_read,
    std_clade = sqrt(proportion * total_read * (1 - total_read * proportion))) |>
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


    q_domain <- "Bacteria"
    bacteria <- stats_df |>
    dplyr::filter(rank == "S") |>
    dplyr::filter(Domain == q_domain) |>
    tidyr::complete(sample_id, name, 
                    fill = list(significance = "Non-significant", 
                                rank = "S", 
                                Domain = q_domain)) |>
    dplyr::select(sample_id, Domain, name, rank, 
                    ratio_clade, padj, significance, n_fragments_clade)

print("qdom")
q_domain <- "Viruses"
viruses <- stats_df |>
dplyr::filter(rank %in% c("S")) |>
dplyr::filter(Domain == q_domain) |>
tidyr::complete(sample_id, name, fill = list(significance = "Non-significant", 
                                            rank = "S", 
                                            Domain = q_domain)) |>
dplyr::select(sample_id, Domain, name, rank, ratio_clade, padj, 
             significance, n_fragments_clade)



    bacterial_plot <- plot_minimisers(bacteria)
    viral_plot <- plot_minimisers(viruses)

    png(glue("{project_dir}/results/sparki/bacterial_minimisers.png"), 
        width = 1500, height = 2600)
    bacterial_plot
    dev.off()

print("end")
    png(glue("{project_dir}/results/sparki/viral_minimisers.png"), 
        width = 1500, height = 1000)
    viral_plot
    dev.off()


# process_cohort(project_dir, ref_db_file)


#  plot <- ggplot2::ggplot(
#             report, 
#             ggplot2::aes(x = colname_taxon, y = log10(colname_n_frag_clade), fill = colname_taxon)
#         ) +
#         ggplot2::geom_violin(scale = "width") +
#         ggplot2::geom_jitter(size = 0.5) +
#         ggplot2::theme_bw() +
#         ggplot2::theme(
#             # x-axis
#             axis.text.x = ggplot2::element_blank(),
#             axis.title.x = ggplot2::element_blank(),
#             axis.ticks.x = ggplot2::element_blank(),
#             strip.text.x = ggplot2::element_text(size = 15),
#             # y-axis
#             axis.text.y = ggplot2::element_text(size = 15),
#             axis.title.y = ggplot2::element_text(size = 16),
#             # legend
#             legend.position = "none"
#         ) +
#         ggplot2::ylab(expression("log"[10]~"(# classified reads)")) +
#         ggplot2::scale_fill_manual(values = c("gold1", "royalblue", "snow2", "indianred2")) +
#         ggplot2::facet_wrap(~colname_taxon, scales = "free")



