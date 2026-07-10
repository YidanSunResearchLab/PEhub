## basilisk environment for the Python leidenalg backend.
##
## Users never touch Python: basilisk provisions an isolated conda environment
## on first use. The versions below are pinned so that community detection is
## reproducible across machines.
##
## IMPORTANT: these must match the environment in which the reference results
## were produced. Check with reticulate::py_list_packages() in that environment
## and edit if they differ.

#' @importFrom basilisk BasiliskEnvironment
NULL

.pehub_env <- basilisk::BasiliskEnvironment(
  envname  = "pehub_leiden_env",
  pkgname  = "PEhub",
  packages = c(
    "python=3.13.5",
    "leidenalg=0.10.2",
    "python-igraph=0.10.8",
    "numpy=1.26.4"
  ),
  channels = c("conda-forge")
)
