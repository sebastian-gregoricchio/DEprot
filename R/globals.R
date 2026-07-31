# Declares the names used in non-standard evaluation (dplyr / ggplot2 aes) so that
# R CMD check does not flag them as "no visible binding for global variable".
# This is a plain top-level call (not roxygen), so it does not affect NAMESPACE.

utils::globalVariables(c(
  # --- PCoA plotting internals ---
  "axis.id", "Proportion.of.Variance", "Cumulative.Proportion", "Broken.stick",
  "loading.x.scaled", "loading.y.scaled", "variable", "shape", "color",
  "observed", "reconstructed.2D"
))

