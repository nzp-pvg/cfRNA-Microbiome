############################################################
## Title   : Plasma cfRNA pathogen profiling – results & figures
## Author  : Zhaoxia Wang
## Version : v0.1.0
## Date    : 2024-09-30
##
## Purpose :
##   - Summarize Kraken2 non-host read fractions (relaxed vs strict)
##     and generate Figure 2.
##   - Summarize MetaPhlAn4 species-level profiles, build
##     phylogeny-aware community trees, compute background
##     signature scores, and generate Figure 3.
##
## NOTE:
##   - All directory paths are placeholders and must be adapted
##     to your own project structure.
##   - Script assumes that upstream steps (fastp, Bowtie2,
##     Kraken2, MetaPhlAn4) have already been completed.
############################################################

rm(list = ls())

#############################
## 0. Load libraries
#############################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(forcats)

## For trees
library(igraph)
library(tidygraph)
library(ggraph)

## For Excel export (optional)
## install.packages("writexl")
library(writexl)

#############################
## 1. User-configurable paths
#############################

## Project-level root (adapt to your own machine)
PROJECT_ROOT <- "/path/to/project"  # e.g. "/data/cfRNA_project"

## 1.1 Kraken2 summary (for Figure 2)
## Assumed to be a TSV with at least:
##   cohort     : "CAD" or "TB"
##   group      : "cad", "ctl", "tb_pos", "tb_neg"  (or similar)
##   sample     : sample ID
##   mode       : "relaxed" or "strict"
##   classified_pct : % of non-host reads classified by Kraken2
KRAKEN_SUMMARY_TSV <- file.path(
  PROJECT_ROOT,
  "results",
  "Kraken2_nonhost_summary.tsv"
)

## 1.2 MetaPhlAn species-level profiles (strict mode)
## These should be post-processed tables that already contain:
##   cohort           : "CAD" or "TB"
##   group            : disease group (e.g. "cad", "ctl", "tb_pos", "tb_neg")
##   sample           : sample ID
##   clade_name       : full MetaPhlAn clade string
##   relative_abundance : species-level relative abundance (%)
##
## You can generate these tables with your own meta-analysis pipeline,
## or adapt the inputs to the format expected here.
MPA_CAD_TSV <- file.path(
  PROJECT_ROOT,
  "results",
  "metaphlan_cad_species_strict.tsv"
)
MPA_TB_TSV <- file.path(
  PROJECT_ROOT,
  "results",
  "metaphlan_tb_species_strict.tsv"
)

## 1.3 Output directory for figures and Excel tables
OUT_DIR_FIG <- file.path(PROJECT_ROOT, "figures")
OUT_DIR_TAB <- file.path(PROJECT_ROOT, "tables")

dir.create(OUT_DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_TAB, showWarnings = FALSE, recursive = TRUE)

#############################
## 2. Helper functions
#############################

## 2.1 Simple ggplot theme
theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text       = element_text(face = "bold"),
      legend.title     = element_text(face = "bold")
    )
}

## 2.2 Utility to parse MetaPhlAn clade_name into taxonomic ranks
## clade_name example:
##   "k__Bacteria|p__Actinobacteriota|c__Actinomycetia|...|s__Cutibacterium_acnes"
parse_metaphlan_clade <- function(clade_name) {
  ranks <- str_split(clade_name, "\\|", simplify = TRUE)
  ranks <- ifelse(ranks == "", NA, ranks)

  out <- tibble(
    clade_name = clade_name,
    kingdom    = ranks[, 1],
    phylum     = ranks[, 2],
    class      = ranks[, 3],
    order      = ranks[, 4],
    family     = ranks[, 5],
    genus      = ranks[, 6],
    species    = ranks[, 7]
  )
  out
}

## 2.3 Build a taxonomic tree from species-level MetaPhlAn table
## Input: data frame with columns:
##   clade_name, cohort, group, sample, relative_abundance
## Output: tbl_graph with node attributes:
##   taxon, rank, is_tip, total_abundance, log10_total, category, label
build_taxonomy_graph <- function(mpa_df,
                                 cohort_name,
                                 group_col = "group",
                                 min_total_abund = 0.01,
                                 label_key_species = NULL) {

  # Parse taxonomy
  tax_parsed <- parse_metaphlan_clade(mpa_df$clade_name) %>%
    bind_cols(mpa_df %>% select(-clade_name))

  # Restrict to species-level entries
  tax_species <- tax_parsed %>%
    filter(!is.na(species),
           str_detect(species, "s__")) %>%
    mutate(species_clean = str_replace(species, "^s__", ""))

  # Compute group-wise means and total abundance for each species
  species_summary <- tax_species %>%
    group_by(species_clean, !!sym(group_col)) %>%
    summarise(
      mean_abund = mean(relative_abundance, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    group_by(species_clean) %>%
    summarise(
      total_abundance = sum(mean_abund, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      log10_total = log10(total_abundance + 1e-6)
    )

  # Build node table for all ranks
  node_df <- tax_species %>%
    select(
      kingdom, phylum, class, order,
      family, genus, species_clean
    ) %>%
    distinct()

  # Stack into long format: one row per node
  node_long <- node_df %>%
    pivot_longer(
      cols = c(kingdom, phylum, class, order, family, genus, species_clean),
      names_to  = "rank",
      values_to = "taxon"
    ) %>%
    filter(!is.na(taxon)) %>%
    distinct()

  # Identify species tips
  node_long <- node_long %>%
    mutate(
      is_tip = rank == "species_clean"
    )

  # Create edges between adjacent ranks
  edge_list <- tax_species %>%
    select(kingdom, phylum, class, order, family, genus, species_clean) %>%
    distinct() %>%
    mutate_all(~ ifelse(is.na(.), NA, .)) %>%
    rowwise() %>%
    do({
      row <- .
      tibble(
        from = c(row$kingdom,
                 row$phylum,
                 row$class,
                 row$order,
                 row$family,
                 row$genus),
        to   = c(row$phylum,
                 row$class,
                 row$order,
                 row$family,
                 row$genus,
                 row$species_clean)
      )
    }) %>%
    ungroup() %>%
    filter(!is.na(from), !is.na(to)) %>%
    distinct()

  # Build graph
  g <- graph_from_data_frame(
    d = edge_list,
    vertices = node_long %>% distinct(taxon, rank, is_tip),
    directed = TRUE
  )

  g_tbl <- as_tbl_graph(g)

  # Add abundance and labels
  g_tbl <- g_tbl %>%
    left_join(
      species_summary,
      by = c("taxon" = "species_clean")
    ) %>%
    mutate(
      cohort  = cohort_name,
      label   = ifelse(is_tip, taxon, NA_character_)
    )

  # If key species supplied, force labels for them
  if (!is.null(label_key_species)) {
    g_tbl <- g_tbl %>%
      mutate(
        label = ifelse(taxon %in% label_key_species, taxon, label)
      )
  }

  g_tbl
}

## 2.4 Helper to classify occurrence pattern between two groups
## E.g. TB_pos vs TB_neg, CAD vs CTL
## Input: species_summary_wide (one row per species)
## with columns mean_abund_grp1, mean_abund_grp2
annotate_occurrence_pattern <- function(df,
                                        grp1_col,
                                        grp2_col,
                                        labels = c("Shared", "Group1 only", "Group2 only")) {
  df %>%
    mutate(
      mean1 = .data[[grp1_col]],
      mean2 = .data[[grp2_col]],
      category = case_when(
        mean1 > 0 & mean2 > 0 ~ labels[1],
        mean1 > 0 & mean2 == 0 ~ labels[2],
        mean1 == 0 & mean2 > 0 ~ labels[3],
        TRUE                   ~ "Other / internal"
      )
    )
}

#############################
## 3. Figure 2 – Kraken2 non-host fractions
#############################

kraken_summary <- read_tsv(
  KRAKEN_SUMMARY_TSV,
  col_types = cols()
)

## Example expected columns:
##   cohort, group, sample, mode, classified_pct
## You may adapt mapping below to your real column names.

# Ensure factors and ordering
kraken_summary <- kraken_summary %>%
  mutate(
    cohort = factor(cohort, levels = c("CAD", "TB")),
    mode   = factor(mode, levels = c("relaxed", "strict"))
  )

# Panel A: violin plots (non-host % classified)
p_fig2A <- kraken_summary %>%
  ggplot(aes(x = interaction(cohort, group, mode),
             y = classified_pct,
             fill = mode)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_jitter(
    position = position_jitter(width = 0.1),
    size = 1,
    alpha = 0.8
  ) +
  facet_wrap(~ cohort, scales = "free_x") +
  scale_fill_manual(values = c("relaxed" = "#4575b4", "strict" = "#d73027")) +
  labs(
    x = "",
    y = "Kraken2-classified non-host reads (%)",
    fill = "Host-filtering mode"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Panel B: delta (strict – relaxed)
kraken_delta <- kraken_summary %>%
  select(cohort, group, sample, mode, classified_pct) %>%
  pivot_wider(
    names_from = mode,
    values_from = classified_pct
  ) %>%
  mutate(
    delta = strict - relaxed
  )

p_fig2B <- kraken_delta %>%
  ggplot(aes(x = interaction(cohort, group),
             y = delta,
             fill = cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(
    position = position_jitter(width = 0.15),
    size = 1,
    alpha = 0.8
  ) +
  labs(
    x = "",
    y = expression(Delta~"Kraken2-classified non-host reads (%)"~"(strict - relaxed)"),
    fill = "Cohort"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

## Save Figure 2 panels
ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure2A_kraken_nonhost_violin.pdf"),
  plot     = p_fig2A,
  width    = 7,
  height   = 4
)

ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure2B_kraken_nonhost_delta.pdf"),
  plot     = p_fig2B,
  width    = 7,
  height   = 4
)

#############################
## 4. Figure 3 – MetaPhlAn trees & signatures
#############################

## 4.1 Load species-level MetaPhlAn tables
mpa_cad <- read_tsv(MPA_CAD_TSV, col_types = cols()) %>%
  mutate(
    cohort = "CAD"
  )

mpa_tb <- read_tsv(MPA_TB_TSV, col_types = cols()) %>%
  mutate(
    cohort = "TB"
  )

## Expected columns: cohort, group, sample, clade_name, relative_abundance

## 4.2 Build CAD and TB graphs (strict mode only)

# Example: define key taxa you want to force labels for
key_species_cad <- c(
  "Cutibacterium_acnes",
  "Rothia_kristinae",
  "Corynebacterium_pseudokroppenstedtii",
  "Bifidobacterium_adolescentis",
  "Bifidobacterium_longum",
  "Collinsella_aerofaciens",
  "Ligilactobacillus_ruminis"
)

key_species_tb <- c(
  "Cutibacterium_acnes",
  "Corynebacterium_tuberculostearicum",
  "Staphylococcus_epidermidis",
  "Saccharomyces_cerevisiae"
)

g_cad <- build_taxonomy_graph(
  mpa_df            = mpa_cad,
  cohort_name       = "CAD",
  group_col         = "group",
  min_total_abund   = 0.01,
  label_key_species = key_species_cad
)

g_tb <- build_taxonomy_graph(
  mpa_df            = mpa_tb,
  cohort_name       = "TB",
  group_col         = "group",
  min_total_abund   = 0.01,
  label_key_species = key_species_tb
)

## 4.3 Compute total abundance and occurrence pattern per species

## CAD: summarize mean abundances per group
cad_species_means <- mpa_cad %>%
  parse_metaphlan_clade() %>%
  filter(!is.na(species),
         str_detect(species, "s__")) %>%
  mutate(species_clean = str_replace(species, "^s__", "")) %>%
  group_by(species_clean, group) %>%
  summarise(
    mean_abund = mean(relative_abundance, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  pivot_wider(
    names_from  = group,
    values_from = mean_abund,
    values_fill = 0
  )

## Example: expecting groups "cad" and "ctl"
cad_species_means <- cad_species_means %>%
  rename(
    mean_cad = cad,
    mean_ctl = ctl
  )

cad_species_means <- annotate_occurrence_pattern(
  cad_species_means,
  grp1_col = "mean_cad",
  grp2_col = "mean_ctl",
  labels   = c("Shared (CAD & control)", "CAD only", "Control only")
)

## TB: summarize means per group (tb_pos, tb_neg)
tb_species_means <- mpa_tb %>%
  parse_metaphlan_clade() %>%
  filter(!is.na(species),
         str_detect(species, "s__")) %>%
  mutate(species_clean = str_replace(species, "^s__", "")) %>%
  group_by(species_clean, group) %>%
  summarise(
    mean_abund = mean(relative_abundance, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  pivot_wider(
    names_from  = group,
    values_from = mean_abund,
    values_fill = 0
  )

## Example: expecting groups "tb_pos" and "tb_neg"
tb_species_means <- tb_species_means %>%
  rename(
    mean_tb_pos = tb_pos,
    mean_tb_neg = tb_neg
  )

tb_species_means <- annotate_occurrence_pattern(
  tb_species_means,
  grp1_col = "mean_tb_pos",
  grp2_col = "mean_tb_neg",
  labels   = c("Shared (TB+ & TB-)", "TB+ only", "TB- only")
)

## 4.4 Attach abundance & category back to graphs

g_cad <- g_cad %>%
  left_join(
    cad_species_means %>%
      transmute(
        taxon          = species_clean,
        mean_cad,
        mean_ctl,
        category_cad   = category
      ),
    by = "taxon"
  )

g_tb <- g_tb %>%
  left_join(
    tb_species_means %>%
      transmute(
        taxon           = species_clean,
        mean_tb_pos,
        mean_tb_neg,
        category_tb     = category
      ),
    by = "taxon"
  )

## 4.5 Circular trees (Figure 3A–B)

## CAD tree
p_cad_tree <- ggraph(g_cad, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(color = "grey70", width = 0.3) +
  geom_node_point(
    aes(size = log10_total),
    color = "grey30",
    alpha = 0.9,
    data = as_tibble(g_cad) %>% filter(is_tip)
  ) +
  geom_node_text(
    aes(label = label),
    size = 2.2,
    repel = FALSE,
    data = as_tibble(g_cad) %>% filter(!is.na(label))
  ) +
  scale_size_continuous(
    name = "Total abundance\n(log10 scale)",
    range = c(1, 6)
  ) +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure3A_CAD_tree.pdf"),
  plot     = p_cad_tree,
  width    = 6,
  height   = 6
)

## TB tree
p_tb_tree <- ggraph(g_tb, layout = "dendrogram", circular = TRUE) +
  geom_edge_diagonal(color = "grey70", width = 0.3) +
  geom_node_point(
    aes(size = log10_total),
    color = "grey30",
    alpha = 0.9,
    data = as_tibble(g_tb) %>% filter(is_tip)
  ) +
  geom_node_text(
    aes(label = label),
    size = 2.2,
    repel = FALSE,
    data = as_tibble(g_tb) %>% filter(!is.na(label))
  ) +
  scale_size_continuous(
    name = "Total abundance\n(log10 scale)",
    range = c(1, 6)
  ) +
  coord_equal() +
  theme_void() +
  theme(
    legend.position = "right"
  )

ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure3B_TB_tree.pdf"),
  plot     = p_tb_tree,
  width    = 6,
  height   = 6
)

## 4.6 Background-derived signature scores (Figure 3C)

## Example: define CAD signatures
cad_sig_cad_enriched <- c(
  "Bifidobacterium_adolescentis",
  "Bifidobacterium_longum",
  "Bifidobacterium_catenulatum",
  "Bifidobacterium_bifidum",
  "Collinsella_aerofaciens",
  "Ligilactobacillus_ruminis"
)

cad_sig_ctl_enriched <- c(
  "Rothia_kristinae",
  "Corynebacterium_pseudokroppenstedtii"
)

cad_sig_single_species <- "Cutibacterium_acnes"

## Prepare per-sample species abundance matrix for CAD
cad_abund <- mpa_cad %>%
  parse_metaphlan_clade() %>%
  filter(!is.na(species),
         str_detect(species, "s__")) %>%
  mutate(species_clean = str_replace(species, "^s__", "")) %>%
  select(cohort, group, sample, species_clean, relative_abundance)

## Compute signature scores
cad_sig_scores <- cad_abund %>%
  mutate(
    species_clean = as.character(species_clean)
  ) %>%
  group_by(cohort, group, sample) %>%
  summarise(
    CAD_signature = sum(
      relative_abundance[species_clean %in% cad_sig_cad_enriched],
      na.rm = TRUE
    ),
    CTL_signature = sum(
      relative_abundance[species_clean %in% cad_sig_ctl_enriched],
      na.rm = TRUE
    ),
    C_acnes = sum(
      relative_abundance[species_clean == cad_sig_single_species],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

## TB background-4 signature (example)
tb_background4 <- c(
  "Cutibacterium_acnes",
  "Corynebacterium_tuberculostearicum",
  "Staphylococcus_epidermidis",
  "Saccharomyces_cerevisiae"
)

tb_abund <- mpa_tb %>%
  parse_metaphlan_clade() %>%
  filter(!is.na(species),
         str_detect(species, "s__")) %>%
  mutate(species_clean = str_replace(species, "^s__", "")) %>%
  select(cohort, group, sample, species_clean, relative_abundance)

tb_sig_scores <- tb_abund %>%
  group_by(cohort, group, sample) %>%
  summarise(
    TB_background4 = sum(
      relative_abundance[species_clean %in% tb_background4],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

## Combine CAD signature scores into long-format for plotting
cad_sig_long <- cad_sig_scores %>%
  pivot_longer(
    cols = c(CAD_signature, CTL_signature, C_acnes),
    names_to  = "signature",
    values_to = "value"
  )

p_cad_sig <- cad_sig_long %>%
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(
    position = position_jitter(width = 0.1),
    size = 1,
    alpha = 0.8
  ) +
  facet_wrap(~ signature, scales = "free_y") +
  labs(
    x = "",
    y = "Relative abundance (%)",
    fill = "Group"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure3C_CAD_signatures.pdf"),
  plot     = p_cad_sig,
  width    = 7,
  height   = 4.5
)

## TB signatures (background4 only; boxplot by TB status)
p_tb_sig <- tb_sig_scores %>%
  ggplot(aes(x = group, y = TB_background4, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.9) +
  geom_jitter(
    position = position_jitter(width = 0.1),
    size = 1,
    alpha = 0.8
  ) +
  labs(
    x = "",
    y = "TB_background4 signature (%)",
    fill = "Group"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(OUT_DIR_FIG, "Figure3D_TB_background4.pdf"),
  plot     = p_tb_sig,
  width    = 5,
  height   = 4
)

#############################
## 5. Export key tables for Prism / Supplement
#############################

## 5.1 Export signature tables
write_xlsx(
  list(
    CAD_signature_scores = cad_sig_scores,
    TB_signature_scores  = tb_sig_scores
  ),
  path = file.path(OUT_DIR_TAB, "signature_scores_for_Prism.xlsx")
)

## 5.2 Export species-level mean abundance tables
write_xlsx(
  list(
    CAD_species_means = cad_species_means,
    TB_species_means  = tb_species_means
  ),
  path = file.path(OUT_DIR_TAB, "metaphlan_species_means.xlsx")
)

############################################################
## End of script
############################################################
