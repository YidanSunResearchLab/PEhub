## Community detection via leidenalg.
##
## DESIGN:
## We call leidenalg directly through reticulate, reproducing exactly the
## behaviour of leiden::leiden() with default arguments for a weighted undirected
## igraph:
##   partition_type  = RBConfigurationVertexPartition
##   weights         = NULL  (leidenalg reads the graph's own 'weight' attribute)
##   n_iterations    = 2L
##   seed            = NULL  (non-deterministic, as in the original code)
##   max_comm_size   = 0L
##   return value    = 1-based membership vector  (part$membership + 1)
##
## Python environment:
## We try, in order:
##   1. basilisk (isolated conda env) — cleanest, but requires conda build tools.
##   2. The current reticulate Python session — works when the user has already
##      activated a conda environment that contains leidenalg (the recommended
##      install path).
## This lets the package work regardless of whether basilisk can build its env.

## Internal cache for the shared basilisk process (used if basilisk succeeds).
.pehub_cache <- new.env(parent = emptyenv())

## Try to start basilisk once; if it fails, fall back to plain reticulate.
.pehub_proc <- function() {
  if (isTRUE(.pehub_cache$basilisk_failed)) return(NULL)
  if (!is.null(.pehub_cache$proc)) return(.pehub_cache$proc)

  proc <- tryCatch(
    basilisk::basiliskStart(.pehub_env),
    error = function(e) {
      message(
        "[PEhub] basilisk could not build its Python environment:\n",
        "  ", conditionMessage(e), "\n",
        "[PEhub] Falling back to the current reticulate Python session.\n",
        "[PEhub] Make sure 'leidenalg' is installed in that session ",
        "(e.g. via the PEhub conda environment)."
      )
      NULL
    }
  )

  if (is.null(proc)) {
    .pehub_cache$basilisk_failed <- TRUE
  } else {
    .pehub_cache$proc <- proc
    reg.finalizer(
      .pehub_cache,
      function(e) if (!is.null(e$proc)) try(basilisk::basiliskStop(e$proc), silent = TRUE),
      onexit = TRUE
    )
  }
  proc
}


#' Shut down the shared PEhub Python process
#'
#' Stops the basilisk process used for Leiden community detection. Called
#' automatically at session end; can also be called manually to free resources.
#' Has no effect when running in fallback (plain reticulate) mode.
#'
#' @return Invisibly \code{NULL}.
#' @export
pehub_shutdown <- function() {
  if (!is.null(.pehub_cache$proc)) {
    try(basilisk::basiliskStop(.pehub_cache$proc), silent = TRUE)
    .pehub_cache$proc <- NULL
  }
  invisible(NULL)
}


## Internal: run leidenalg on a weighted adjacency matrix.
## Returns a 1-based integer membership vector in vertex order.
.leiden_basilisk <- function(adj, resolution_parameter = 1.0,
                             n_iterations = 2L, seed = NULL) {

  ## shared function that does the actual Python work ---
  .run <- function(ig, la, adj_mat, res, n_iter, sd) {
    gpy <- ig$Graph$Weighted_Adjacency(
      adj_mat, mode = "undirected", attr = "weight", loops = FALSE
    )
    part <- la$find_partition(
      gpy,
      la$RBConfigurationVertexPartition,
      initial_membership   = NULL,
      weights              = NULL,
      resolution_parameter = res,
      seed                 = sd,
      n_iterations         = n_iter,
      max_comm_size        = 0L
    )
    as.integer(part$membership) + 1L
  }

  proc <- .pehub_proc()

  if (!is.null(proc)) {
    ## --- basilisk path (isolated conda env) ---
    membership <- basilisk::basiliskRun(
      proc,
      fun = function(adj_mat, res, n_iter, sd) {
        ig <- reticulate::import("igraph",    convert = TRUE)
        la <- reticulate::import("leidenalg", convert = TRUE)
        .run(ig, la, adj_mat, res, n_iter, sd)
      },
      adj_mat = adj,
      res     = resolution_parameter,
      n_iter  = as.integer(n_iterations),
      sd      = if (is.null(seed)) NULL else as.integer(seed)
    )
  } else {
    ## --- fallback: plain reticulate (user's active Python) ---
    if (!reticulate::py_module_available("leidenalg")) {
      stop(
        "The Python 'leidenalg' module is not available.\n",
        "Please activate a Python environment that has it installed, e.g.:\n",
        "  mamba activate pehub\n",
        "or install it with:\n",
        "  reticulate::py_install('leidenalg', envname = 'r-reticulate')",
        call. = FALSE
      )
    }
    ig <- reticulate::import("igraph",    convert = TRUE)
    la <- reticulate::import("leidenalg", convert = TRUE)
    membership <- .run(ig, la, adj, resolution_parameter,
                       as.integer(n_iterations),
                       if (is.null(seed)) NULL else as.integer(seed))
  }

  as.integer(membership)
}


#' Run Leiden community detection for one promoter's co-membership graph
#'
#' Partitions a promoter's enhancer co-membership matrix into communities using
#' the Leiden algorithm (Python \code{leidenalg}), then ranks communities by
#' total node strength (dominance). Python is provisioned automatically via
#' basilisk when possible; otherwise the current reticulate Python session is
#' used (see the installation instructions).
#'
#' @param promoter_id Character scalar, the promoter identifier.
#' @param comembership A square weighted adjacency matrix of enhancer co-membership.
#' @param resolution Numeric, the Leiden resolution parameter (default 0.5).
#' @param seed Integer, used for \code{set.seed()} on the R side (as in the
#'   original implementation). The Python RNG is not seeded by default,
#'   matching \code{leiden::leiden()}'s behaviour.
#'
#' @return A tibble with columns \code{promoter_id}, \code{hub_id},
#'   \code{enhancers} (list column) and \code{dominance}.
#'
#' @importFrom igraph graph_from_adjacency_matrix V vcount ecount strength E
#' @importFrom tibble tibble
#' @importFrom dplyr mutate group_by summarise arrange desc "%>%"
#' @export
run_leiden_for_promoter <- function(promoter_id,
                                    comembership,
                                    resolution = 0.5,
                                    seed = 42) {

  g <- igraph::graph_from_adjacency_matrix(
    comembership, mode = "undirected", weighted = TRUE
  )

  enh_names <- igraph::V(g)$name

  if (length(enh_names) == 0 || igraph::vcount(g) < 3 || igraph::ecount(g) == 0) {
    return(tibble::tibble(
      promoter_id = promoter_id,
      hub_id      = paste0(promoter_id, "_hub1"),
      enhancers   = list(enh_names),
      dominance   = 0
    ))
  }

  set.seed(seed)
  cl <- .leiden_basilisk(as.matrix(comembership), resolution_parameter = resolution, seed = seed)

  cl_df <- tibble::tibble(
    enhancer_id = enh_names,
    cluster     = as.integer(cl)
  )

  node_strength <- igraph::strength(g, vids = igraph::V(g), weights = igraph::E(g)$weight)

  hub_strength <- cl_df %>%
    dplyr::mutate(node_strength = node_strength[enhancer_id]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(dominance = sum(node_strength, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(dominance))

  clusters <- split(cl_df$enhancer_id, cl_df$cluster)
  tibble::tibble(
    promoter_id = promoter_id,
    hub_id      = paste0(promoter_id, "_hub", seq_along(clusters)),
    enhancers   = unname(clusters),
    dominance   = hub_strength$dominance[match(as.integer(names(clusters)),
                                               hub_strength$cluster)]
  )
}
