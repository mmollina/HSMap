# Register non-standard-evaluation column names used in ggplot2/dplyr pipelines so
# R CMD check does not report "no visible binding for global variable".
utils::globalVariables(c(
  "marker", "ykey", "dose_f", "value", "stage", "interval", "dam", "phase",
  "alpha", "pos", "block", "segment"
))
