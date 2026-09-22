### ============================================================
### RUN FULL PIPELINE
### ============================================================
# Runs every analysis script in dependency order. Set the working
# directory to the analysis root first. Each script is run in a fresh
# environment so leftover variables from one step can't leak into the next.

scripts <- c(
  "03_normals.R",
  "04_genomic_index.R",
  "05_genome_modified.R",
  "06_gene_concordance.R",
  "07_accuracy.R",
  "08_case_graphs.R",
  "09_correlations.R"   # depends on CSVs written by 04, 05 and 07
)

for (s in scripts) {
  message("Running ", s, " ...")
  source(paste0(getwd(), "/Analytic Code/", s), local = new.env())
}
