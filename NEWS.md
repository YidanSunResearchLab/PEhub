# PEhub (development version)

## New features

* New `pehub_check_power()`: a pre-flight diagnostic for pooled low-depth /
  single-cell input. PEhub's background model and significance testing were
  developed and validated on deeply sequenced bulk HiChIP / Pore-C data; this
  function checks total and per-distance-bin read depth in a `loop_file_all`
  and flags samples (`"low_power"`) or specific distance ranges
  (`"long_range_unreliable"`) where hub calling may not be reliable,
  especially beyond ~50kb. It changes no downstream computation — it is a
  diagnostic, not a correction.
* `pehub_prepare_interactions()` and `pehub_run()` gain `n_cells` (optional,
  informational) and `check_power` (default `TRUE`) arguments. When enabled,
  the power check runs automatically and its result is attached to the
  returned object as `$power_check`. Fully backward compatible: existing
  calls without these arguments behave exactly as before, aside from the
  added diagnostic message.
* New vignette section, "Single-cell / pooled low-depth input", documenting
  this gap explicitly and showing how to use the new check.

## Motivation

Found while applying PEhub to snm3C-seq (single-nucleus methyl-3C) data
pooled per donor and per age bracket: the method's reconstruction logic
(recombining the single enhancer-promoter contact each cell shows into a
population-level hub) is not inherently incompatible with pooling single
cells, but the package had no documented guidance on the minimum depth this
requires, and no validation data anywhere near that regime (its own
validation sets are bulk human-heart H3K27ac HiChIP and GM12878 Pore-C).
Empirically, two real single-donor calibration runs at 658K vs. 5.35M valid
contact pairs showed an ~8x depth difference produce a ~21x difference in
usable gene count, and a background-model fit to single-donor-pooled data was
found to become numerically degenerate beyond ~50kb. The new check surfaces
this up front rather than leaving it to be discovered downstream.
