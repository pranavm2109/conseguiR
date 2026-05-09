#' @keywords internal
.conseguiR_backend_dir <- function(create = FALSE) {
  opt_dir <- getOption("conseguiR.backend_dir")
  if (is.character(opt_dir) && length(opt_dir) == 1L && nzchar(opt_dir)) {
    if (isTRUE(create)) {
      dir.create(opt_dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(normalizePath(opt_dir, winslash = "/", mustWork = FALSE))
  }

  candidates <- c(
    if (!is.null(.conseguiR_state$pkg_root)) file.path(.conseguiR_state$pkg_root, "data", "processed"),
    file.path(getwd(), "data", "processed"),
    file.path(tempdir(), "conseguiR_backend")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  existing <- candidates[file.exists(candidates)]
  chosen <- if (length(existing) > 0L) existing[[1]] else candidates[[length(candidates)]]

  if (isTRUE(create)) {
    dir.create(chosen, recursive = TRUE, showWarnings = FALSE)
  }

  normalizePath(chosen, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.conseguiR_backend_paths <- function(backend_dir = .conseguiR_backend_dir(create = TRUE)) {
  paths <- list(
    backend_dir = backend_dir,
    gene_reg_graph_rds = file.path(backend_dir, "gene_reg_graph_no_scores.rds"),
    gene_reg_graph_nodes = file.path(backend_dir, "gene_reg_graph_no_scores_nodes.tsv.gz"),
    gene_reg_graph_edges = file.path(backend_dir, "gene_reg_graph_no_scores_edges.tsv.gz"),
    gene_gene_graph_rds = file.path(backend_dir, "gene_gene_graph.rds"),
    gene_gene_graph_nodes = file.path(backend_dir, "gene_gene_graph_nodes.tsv.gz"),
    gene_gene_graph_edges = file.path(backend_dir, "gene_gene_graph_edges.tsv.gz")
  )
  paths
}

#' @keywords internal
.conseguiR_load_backend_graph <- function(kind = c("gene_reg", "gene_gene"), backend_dir = .conseguiR_backend_dir(create = TRUE)) {
  kind <- match.arg(kind)
  paths <- .conseguiR_backend_paths(backend_dir)
  graph_path <- switch(
    kind,
    gene_reg = paths$gene_reg_graph_rds,
    gene_gene = paths$gene_gene_graph_rds
  )
  if (!file.exists(graph_path)) {
    return(NULL)
  }

  graph <- readRDS(graph_path)
  if (!inherits(graph, "igraph")) {
    stop("Expected an igraph object at backend graph path: ", graph_path)
  }
  graph
}

#' @keywords internal
.conseguiR_existing_graph_paths <- function(paths) {
  flat <- unlist(paths, use.names = FALSE)
  flat[file.exists(flat)]
}

#' @keywords internal
.conseguiR_backend_resource_work_dir <- function(create = FALSE) {
  path <- file.path(tempdir(), "conseguiR_backend_resources")
  if (isTRUE(create)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.conseguiR_backend_resource_dirs <- function(include_work_dir = TRUE) {
  option_dir <- getOption("conseguiR.backend_resource_dir")
  if (!is.null(option_dir)) {
    option_dir <- trimws(as.character(option_dir))
    option_dir <- option_dir[nzchar(option_dir)]
  }
  backend_dir <- getOption("conseguiR.backend_dir", NULL)
  if (!is.null(backend_dir)) {
    backend_dir <- trimws(as.character(backend_dir))
    backend_dir <- backend_dir[nzchar(backend_dir)]
  }

  installed <- tryCatch(
    system.file("extdata", "backend", package = "conseguiR"),
    error = function(e) ""
  )
  candidates <- c(
    if (isTRUE(include_work_dir)) .conseguiR_backend_resource_work_dir(),
    option_dir,
    backend_dir,
    installed,
    if (!is.null(.conseguiR_state$pkg_root)) {
      file.path(.conseguiR_state$pkg_root, "extdata", "backend")
    },
    if (!is.null(.conseguiR_state$pkg_root)) {
      file.path(.conseguiR_state$pkg_root, "inst", "extdata", "backend")
    },
    file.path(getwd(), "extdata", "backend"),
    file.path(getwd(), "inst", "extdata", "backend")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates[dir.exists(candidates)]
}

#' @keywords internal
.conseguiR_backend_seed_files <- function(kind = c("gene_reg", "gene_gene")) {
  kind <- match.arg(kind)
  switch(
    kind,
    gene_reg = c(
      "gene_reg_graph_no_scores_nodes.tsv.gz",
      "gene_reg_graph_no_scores_edges.tsv.gz",
      "gene_reg_graph_no_scores_nodes_compact.tsv.xz",
      "gene_reg_graph_no_scores_edges_compact.tsv.xz"
    ),
    gene_gene = c(
      "gene_gene_graph_nodes.tsv.gz",
      "gene_gene_graph_edges.tsv.gz",
      "gene_gene_graph_nodes_compact.tsv.xz",
      "gene_gene_graph_edges_compact.tsv.xz"
    )
  )
}

#' @keywords internal
.conseguiR_backend_seed_dir <- function(kind = c("gene_reg", "gene_gene")) {
  kind <- match.arg(kind)
  existing <- .conseguiR_backend_resource_dirs(include_work_dir = TRUE)
  seed_files <- .conseguiR_backend_seed_files(kind)
  has_seed <- function(dir_path) {
    any(file.exists(file.path(dir_path, seed_files)))
  }
  seeded_dirs <- existing[vapply(existing, has_seed, logical(1))]
  if (length(seeded_dirs) == 0L) {
    return(NULL)
  }
  normalizePath(seeded_dirs[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_backend_resource_candidates <- function(filename) {
  dirs <- .conseguiR_backend_resource_dirs(include_work_dir = TRUE)
  if (length(dirs) == 0L) {
    return(character())
  }

  candidates <- file.path(dirs, filename)
  unique(candidates[file.exists(candidates)])
}

#' @keywords internal
.conseguiR_backend_resource_path <- function(filename) {
  candidates <- .conseguiR_backend_resource_candidates(filename)
  if (length(candidates) == 0L) {
    return(NULL)
  }

  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_prefer_repo_relative_path <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(relpath)
  }

  NULL
}

#' @keywords internal
.conseguiR_default_gene_loc_path <- function() {
  path <- .conseguiR_backend_resource_path("NCBI38.gene.loc")
  if (!is.null(path)) {
    return(path)
  }
  path <- .conseguiR_find_data_path("data/raw/NCBI38/NCBI38.gene.loc")
  if (!is.null(path)) {
    return(path)
  }
  path <- .conseguiR_prefer_repo_relative_path("data/raw/NCBI38/NCBI38.gene.loc")
  if (!is.null(path)) {
    return(normalizePath(path, winslash = "/", mustWork = TRUE))
  }
  NULL
}

#' @keywords internal
.conseguiR_default_encode_ccre_path <- function() {
  path <- .conseguiR_prefer_repo_relative_path("data/raw/ENCODE/GRCh38-cCREs.bed")
  if (!is.null(path)) {
    return(path)
  }
  .conseguiR_find_data_path("data/raw/ENCODE/GRCh38-cCREs.bed")
}

#' @keywords internal
.conseguiR_default_encode_gene_links_path <- function() {
  path <- .conseguiR_prefer_repo_relative_path("data/raw/ENCODE/Human-Gene-Links.zip")
  if (!is.null(path)) {
    return(path)
  }
  .conseguiR_find_data_path("data/raw/ENCODE/Human-Gene-Links.zip")
}

#' @keywords internal
.conseguiR_materialize_encode_reg_loc <- function(ccre_bed_path) {
  if (is.null(ccre_bed_path) || !file.exists(ccre_bed_path)) {
    return(NULL)
  }

  resource_dir <- tryCatch(
    .conseguiR_backend_resource_work_dir(create = TRUE),
    error = function(e) NULL
  )
  if (is.null(resource_dir) || !dir.exists(resource_dir)) {
    resource_dir <- .conseguiR_backend_dir(create = TRUE)
  }
  loc_path <- tempfile(pattern = "GRCh38-cCREs_", fileext = ".loc", tmpdir = resource_dir)

  ccre_dt <- data.table::fread(ccre_bed_path, header = FALSE, showProgress = FALSE)
  if (ncol(ccre_dt) < 5L) {
    stop("ENCODE cCRE BED must have at least 5 columns to create a loc file.")
  }

  loc_dt <- data.table::data.table(
    reg_elem_id = as.character(ccre_dt[[5L]]),
    chrom = sub("^chr", "", as.character(ccre_dt[[1L]])),
    start = as.integer(ccre_dt[[2L]]) + 1L,
    end = as.integer(ccre_dt[[3L]])
  )
  loc_dt <- unique(loc_dt)
  loc_dt <- loc_dt[
    !is.na(loc_dt[["reg_elem_id"]]) & loc_dt[["reg_elem_id"]] != "" &
      !is.na(loc_dt[["chrom"]]) & loc_dt[["chrom"]] != "" &
      !is.na(loc_dt[["start"]]) & !is.na(loc_dt[["end"]])
  ]
  data.table::fwrite(loc_dt[, .(reg_elem_id, chrom, start, end)], loc_path, sep = "\t", col.names = FALSE)
  normalizePath(loc_path, winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_reg_loc_from_nodes <- function(nodes) {
  dt <- data.table::as.data.table(nodes)
  if (!"node_type" %in% names(dt) || !"node_id" %in% names(dt)) {
    stop("Node table must contain `node_type` and `node_id` to derive a regulatory loc file.")
  }

  reg_dt <- dt[node_type == "reg"]
  if (nrow(reg_dt) == 0L) {
    return(data.table::data.table(
      reg_elem_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      reg_elem_name = character()
    ))
  }

  chr_col <- if ("reg_chr" %in% names(reg_dt)) {
    data.table::fifelse(!is.na(reg_dt$reg_chr) & reg_dt$reg_chr != "", reg_dt$reg_chr, reg_dt$chr)
  } else if ("chr" %in% names(reg_dt)) {
    reg_dt$chr
  } else {
    rep(NA_character_, nrow(reg_dt))
  }

  start_col <- if ("reg_start" %in% names(reg_dt)) {
    data.table::fifelse(!is.na(reg_dt$reg_start), reg_dt$reg_start, reg_dt$start)
  } else if ("start" %in% names(reg_dt)) {
    reg_dt$start
  } else {
    rep(NA_integer_, nrow(reg_dt))
  }

  end_col <- if ("reg_end" %in% names(reg_dt)) {
    data.table::fifelse(!is.na(reg_dt$reg_end), reg_dt$reg_end, reg_dt$end)
  } else if ("end" %in% names(reg_dt)) {
    reg_dt$end
  } else {
    rep(NA_integer_, nrow(reg_dt))
  }

  reg_name_col <- if ("node_label" %in% names(reg_dt)) {
    as.character(reg_dt$node_label)
  } else {
    rep(NA_character_, nrow(reg_dt))
  }

  loc_dt <- data.table::data.table(
    reg_elem_id = as.character(reg_dt$node_id),
    chrom = sub("^chr", "", as.character(chr_col)),
    start = as.integer(start_col),
    end = as.integer(end_col),
    reg_elem_name = reg_name_col
  )
  unique(loc_dt)[
    !is.na(reg_elem_id) & reg_elem_id != "" &
      !is.na(chrom) & chrom != "" &
      !is.na(start) & !is.na(end)
  ]
}

#' @keywords internal
.conseguiR_is_valid_magma_loc <- function(loc_path, max_rows = 50L) {
  if (is.null(loc_path) || !file.exists(loc_path) || is.na(file.info(loc_path)$size) || file.info(loc_path)$size <= 0L) {
    return(FALSE)
  }

  loc_dt <- tryCatch(
    data.table::fread(loc_path, header = FALSE, nrows = max_rows, sep = "\t", showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(loc_dt) || nrow(loc_dt) == 0L) {
    return(FALSE)
  }

  n_cols <- ncol(loc_dt)
  if (!(n_cols %in% c(4L, 5L))) {
    return(FALSE)
  }

  id_col <- as.character(loc_dt[[1L]])
  chr_col <- as.character(loc_dt[[2L]])
  start_col <- suppressWarnings(as.integer(loc_dt[[3L]]))
  end_col <- suppressWarnings(as.integer(loc_dt[[4L]]))

  core_ok <- !is.na(id_col) & id_col != "" &
    !is.na(chr_col) & chr_col != "" &
    !is.na(start_col) & !is.na(end_col)

  if (!all(core_ok)) {
    return(FALSE)
  }

  if (n_cols == 5L) {
    strand_col <- trimws(as.character(loc_dt[[5L]]))
    if (!all(strand_col %in% c("+", "-"))) {
      return(FALSE)
    }
  }

  TRUE
}

#' @keywords internal
.conseguiR_write_reg_loc <- function(nodes, loc_path) {
  loc_dt <- .conseguiR_reg_loc_from_nodes(nodes)
  dir.create(dirname(loc_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(loc_dt[, .(reg_elem_id, chrom, start, end)], loc_path, sep = "\t", col.names = FALSE)
  normalizePath(loc_path, winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_backend_reg_loc_path <- function(backend_dir = .conseguiR_backend_dir(create = TRUE)) {
  loc_path <- file.path(backend_dir, "GRCh38-cCREs.loc")
  if (.conseguiR_is_valid_magma_loc(loc_path)) {
    return(normalizePath(loc_path, winslash = "/", mustWork = TRUE))
  }

  nodes_path <- file.path(backend_dir, "gene_reg_graph_no_scores_nodes.tsv.gz")
  if (file.exists(nodes_path)) {
    nodes <- data.table::as.data.table(data.table::fread(nodes_path, showProgress = FALSE))
    return(.conseguiR_write_reg_loc(nodes, loc_path))
  }

  seed_dir <- .conseguiR_backend_seed_dir()
  if (!is.null(seed_dir)) {
    compact_nodes <- file.path(seed_dir, "gene_reg_graph_no_scores_nodes_compact.tsv.xz")
    if (file.exists(compact_nodes)) {
      nodes <- .conseguiR_read_backend_table(compact_nodes)
      return(.conseguiR_write_reg_loc(nodes, loc_path))
    }
  }

  NULL
}

#' @keywords internal
.conseguiR_default_reg_loc_path <- function() {
  encode_ccre_path <- .conseguiR_default_encode_ccre_path()
  encode_loc_path <- .conseguiR_materialize_encode_reg_loc(encode_ccre_path)
  if (!is.null(encode_loc_path)) {
    return(encode_loc_path)
  }
  backend_loc_path <- .conseguiR_backend_reg_loc_path()
  if (!is.null(backend_loc_path)) {
    return(backend_loc_path)
  }
  .conseguiR_backend_resource_path("GRCh38-cCREs.loc")
}

#' @keywords internal
.conseguiR_seed_backend_graph <- function(kind = c("gene_reg", "gene_gene"), backend_dir) {
  kind <- match.arg(kind)
  seed_dir <- .conseguiR_backend_seed_dir(kind)
  if (is.null(seed_dir)) {
    return(FALSE)
  }

  compact_materialized <- switch(
    kind,
    gene_reg = .conseguiR_materialize_compact_gene_reg_seed(seed_dir, backend_dir),
    gene_gene = .conseguiR_materialize_compact_gene_gene_seed(seed_dir, backend_dir)
  )
  if (isTRUE(compact_materialized)) {
    return(TRUE)
  }

  files <- switch(
    kind,
    gene_reg = c(
      "gene_reg_graph_no_scores_nodes.tsv.gz",
      "gene_reg_graph_no_scores_edges.tsv.gz",
      "gene_reg_graph_no_scores.rds"
    ),
    gene_gene = c(
      "gene_gene_graph_nodes.tsv.gz",
      "gene_gene_graph_edges.tsv.gz",
      "gene_gene_graph.rds"
    )
  )

  src <- file.path(seed_dir, files)
  existing_src <- src[file.exists(src)]
  if (length(existing_src) < 2L) {
    return(FALSE)
  }

  dest <- file.path(backend_dir, basename(existing_src))
  dir.create(backend_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(existing_src, dest, overwrite = TRUE)
  all(ok)
}

#' @keywords internal
.conseguiR_read_backend_table <- function(path) {
  if (grepl("\\.xz$", path, ignore.case = TRUE)) {
    con <- xzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    return(
      data.table::as.data.table(
        utils::read.delim(
          con,
          sep = "\t",
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      )
    )
  }

  data.table::as.data.table(data.table::fread(path, showProgress = FALSE))
}

#' @keywords internal
.conseguiR_compact_seed_paths <- function(kind = c("gene_reg", "gene_gene"), seed_dir) {
  kind <- match.arg(kind)
  switch(
    kind,
    gene_reg = list(
      nodes = file.path(seed_dir, "gene_reg_graph_no_scores_nodes_compact.tsv.xz"),
      edges = file.path(seed_dir, "gene_reg_graph_no_scores_edges_compact.tsv.xz"),
      rds = file.path(seed_dir, "gene_reg_graph_no_scores.rds")
    ),
    gene_gene = list(
      nodes = file.path(seed_dir, "gene_gene_graph_nodes_compact.tsv.xz"),
      edges = file.path(seed_dir, "gene_gene_graph_edges_compact.tsv.xz"),
      rds = file.path(seed_dir, "gene_gene_graph.rds")
    )
  )
}

#' @keywords internal
.conseguiR_compact_node_index <- function(nodes) {
  dt <- data.table::copy(data.table::as.data.table(nodes))
  if (!"node_index" %in% names(dt)) {
    dt[, node_index := seq_len(.N)]
  }
  dt[, node_index := as.integer(node_index)]
  dt
}

#' @keywords internal
.conseguiR_write_materialized_graph <- function(graph, nodes, edges, paths) {
  dir.create(dirname(paths$nodes), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(nodes, paths$nodes, sep = "\t")
  data.table::fwrite(edges, paths$edges, sep = "\t")
  saveRDS(graph, paths$rds)
  invisible(paths)
}

#' @keywords internal
.conseguiR_materialize_compact_gene_reg_seed <- function(seed_dir, backend_dir) {
  compact <- .conseguiR_compact_seed_paths("gene_reg", seed_dir)
  if (!file.exists(compact$nodes) || !file.exists(compact$edges)) {
    return(FALSE)
  }

  nodes <- .conseguiR_compact_node_index(
    .conseguiR_read_backend_table(compact$nodes)
  )
  edges_compact <- .conseguiR_read_backend_table(compact$edges)
  required_edge_cols <- c("from_idx", "to_idx", "confidence")
  missing_edge_cols <- setdiff(required_edge_cols, names(edges_compact))
  if (length(missing_edge_cols) > 0L) {
    stop(
      "Compact gene-reg seed edges are missing required columns: ",
      paste(missing_edge_cols, collapse = ", ")
    )
  }

  node_lookup <- nodes[, .(
    node_index,
    node_id = as.character(node_id),
    node_type = as.character(node_type),
    chr = if ("chr" %in% names(nodes)) as.character(chr) else NA_character_,
    start = if ("start" %in% names(nodes)) as.integer(start) else NA_integer_,
    end = if ("end" %in% names(nodes)) as.integer(end) else NA_integer_,
    reg_chr = if ("reg_chr" %in% names(nodes)) as.character(reg_chr) else NA_character_,
    reg_start = if ("reg_start" %in% names(nodes)) as.integer(reg_start) else NA_integer_,
    reg_end = if ("reg_end" %in% names(nodes)) as.integer(reg_end) else NA_integer_
  )]

  reg_lookup <- node_lookup[node_type == "reg", .(
    from_idx = node_index,
    from = node_id,
    reg_chr = data.table::fifelse(!is.na(reg_chr) & reg_chr != "", reg_chr, chr),
    reg_start = data.table::fifelse(!is.na(reg_start), reg_start, start),
    reg_end = data.table::fifelse(!is.na(reg_end), reg_end, end)
  )]
  gene_lookup <- node_lookup[node_type == "gene", .(
    to_idx = node_index,
    to = node_id,
    gene_chr = chr,
    gene_start = start,
    gene_end = end
  )]

  edges <- merge(edges_compact, reg_lookup, by = "from_idx", all.x = TRUE, sort = FALSE)
  edges <- merge(edges, gene_lookup, by = "to_idx", all.x = TRUE, sort = FALSE)
  if (!"weight" %in% names(edges)) {
    edges[, weight := data.table::fifelse(
      is.na(confidence),
      NA_real_,
      1 / (1 + as.numeric(confidence))
    )]
  }
  if (!"link_score" %in% names(edges)) {
    edges[, link_score := as.numeric(confidence)]
  }
  if (!"link_method" %in% names(edges)) {
    edges[, link_method := NA_character_]
  }
  edges <- edges[, .(
    from,
    to,
    weight = as.numeric(weight),
    confidence = as.numeric(confidence),
    link_score = as.numeric(link_score),
    link_method = as.character(link_method),
    reg_chr,
    reg_start = as.integer(reg_start),
    reg_end = as.integer(reg_end),
    gene_chr,
    gene_start = as.integer(gene_start),
    gene_end = as.integer(gene_end)
  )]
  reg_target_labels <- edges[, .(
    label = paste(sort(unique(to)), collapse = "|")
  ), by = .(reg_elem_id = from)]

  vertices <- as.data.frame(nodes[, !"node_index"])
  if (!"name" %in% names(vertices)) {
    vertices$name <- vertices$node_id
  }
  graph <- igraph::graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = vertices,
    directed = TRUE
  )

  .conseguiR_write_materialized_graph(
    graph = graph,
    nodes = nodes[, !"node_index"],
    edges = edges,
    paths = list(
      nodes = file.path(backend_dir, "gene_reg_graph_no_scores_nodes.tsv.gz"),
      edges = file.path(backend_dir, "gene_reg_graph_no_scores_edges.tsv.gz"),
      rds = file.path(backend_dir, "gene_reg_graph_no_scores.rds")
    )
  )
  .conseguiR_write_reg_loc(
    nodes = nodes[, !"node_index"],
    loc_path = file.path(backend_dir, "GRCh38-cCREs.loc")
  )
  data.table::fwrite(
    reg_target_labels,
    file.path(backend_dir, "reg_target_labels.tsv.gz"),
    sep = "\t"
  )
  TRUE
}

#' @keywords internal
.conseguiR_materialize_compact_gene_gene_seed <- function(seed_dir, backend_dir) {
  compact <- .conseguiR_compact_seed_paths("gene_gene", seed_dir)
  if (!file.exists(compact$nodes) || !file.exists(compact$edges)) {
    return(FALSE)
  }

  nodes <- .conseguiR_compact_node_index(
    .conseguiR_read_backend_table(compact$nodes)
  )
  edges_compact <- .conseguiR_read_backend_table(compact$edges)
  required_edge_cols <- c("from_idx", "to_idx", "confidence")
  missing_edge_cols <- setdiff(required_edge_cols, names(edges_compact))
  if (length(missing_edge_cols) > 0L) {
    stop(
      "Compact gene-gene seed edges are missing required columns: ",
      paste(missing_edge_cols, collapse = ", ")
    )
  }

  lookup <- nodes[, .(
    node_index,
    node_id = as.character(node_id)
  )]
  from_lookup <- lookup[, .(from_idx = node_index, from = node_id)]
  to_lookup <- lookup[, .(to_idx = node_index, to = node_id)]

  edges <- merge(edges_compact, from_lookup, by = "from_idx", all.x = TRUE, sort = FALSE)
  edges <- merge(edges, to_lookup, by = "to_idx", all.x = TRUE, sort = FALSE)
  if (!"weight" %in% names(edges)) {
    edges[, weight := data.table::fifelse(
      is.na(confidence),
      NA_real_,
      1 / (1 + as.numeric(confidence))
    )]
  }
  if (!"n_protein_edges" %in% names(edges)) {
    edges[, n_protein_edges := NA_integer_]
  }
  edges <- edges[, .(
    from,
    to,
    confidence = as.numeric(confidence),
    weight = as.numeric(weight),
    n_protein_edges = as.integer(n_protein_edges)
  )]

  vertices <- as.data.frame(nodes[, !"node_index"])
  if (!"name" %in% names(vertices)) {
    vertices$name <- vertices$node_id
  }
  graph <- igraph::graph_from_data_frame(
    d = as.data.frame(edges),
    vertices = vertices,
    directed = FALSE
  )

  .conseguiR_write_materialized_graph(
    graph = graph,
    nodes = nodes[, !"node_index"],
    edges = edges,
    paths = list(
      nodes = file.path(backend_dir, "gene_gene_graph_nodes.tsv.gz"),
      edges = file.path(backend_dir, "gene_gene_graph_edges.tsv.gz"),
      rds = file.path(backend_dir, "gene_gene_graph.rds")
    )
  )
  TRUE
}

#' @keywords internal
.conseguiR_data_candidates <- function(relpath) {
  unique(c(
    if (!is.null(.conseguiR_state$pkg_root)) file.path(.conseguiR_state$pkg_root, relpath),
    file.path(getwd(), relpath)
  ))
}

#' @keywords internal
.conseguiR_find_data_path <- function(relpath) {
  candidates <- .conseguiR_data_candidates(relpath)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    return(NULL)
  }
  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_internal_script_env <- function(relpath) {
  env <- new.env(parent = globalenv())
  sys.source(.conseguiR_runtime_file(relpath), envir = env)
  env
}

#' @keywords internal
.conseguiR_graph_files_exist <- function(paths) {
  all(file.exists(unlist(paths, use.names = FALSE)))
}

#' @keywords internal
.conseguiR_initialize_gene_reg_backend <- function(
  paths,
  backend_dir,
  force,
  strict
) {
  target_paths <- list(
    paths$gene_reg_graph_nodes,
    paths$gene_reg_graph_edges
  )

  if (!isTRUE(force) && .conseguiR_graph_files_exist(target_paths)) {
    return("reused")
  }

  if (.conseguiR_seed_backend_graph("gene_reg", backend_dir = backend_dir)) {
    return("seeded")
  }

  msg <- paste(
    "Unable to initialize the backend gene-regulatory graph because no",
    "packaged backend seed was available."
  )
  if (isTRUE(strict)) {
    stop(msg)
  }
  "skipped"
}

#' @keywords internal
.conseguiR_initialize_gene_gene_backend <- function(
  paths,
  backend_dir,
  force,
  strict
) {
  target_paths <- list(
    paths$gene_gene_graph_nodes,
    paths$gene_gene_graph_edges
  )

  if (!isTRUE(force) && .conseguiR_graph_files_exist(target_paths)) {
    return("reused")
  }

  if (.conseguiR_seed_backend_graph("gene_gene", backend_dir = backend_dir)) {
    return("seeded")
  }

  msg <- paste(
    "Unable to initialize the backend gene-gene graph because no",
    "packaged backend seed was available."
  )
  if (isTRUE(strict)) {
    stop(msg)
  }
  "skipped"
}

#' @keywords internal
.conseguiR_backend_init_result <- function(
  backend_dir,
  paths,
  status,
  force
) {
  result <- structure(
    list(
      backend_dir = backend_dir,
      output_paths = paths,
      status = status
    ),
    class = c("conseguiR_backend_init", "list")
  )
  result
}

#' @keywords internal
.conseguiR_backend_init_context <- function(
  backend_dir,
  build_gene_reg,
  build_gene_gene,
  force,
  strict
) {
  if (!is.null(backend_dir)) {
    options(conseguiR.backend_dir = backend_dir)
  }

  resolved_backend_dir <- .conseguiR_backend_dir(create = TRUE)
  options(conseguiR.backend_dir = resolved_backend_dir)
  paths <- .conseguiR_backend_paths(resolved_backend_dir)

  list(
    backend_dir = resolved_backend_dir,
    paths = paths
  )
}

#' @keywords internal
.conseguiR_backend_init_status <- function(
  paths,
  backend_dir,
  build_gene_reg,
  build_gene_gene,
  force,
  strict
) {
  status <- list()

  if (isTRUE(build_gene_reg)) {
    status$gene_reg <- .conseguiR_initialize_gene_reg_backend(
      paths = paths,
      backend_dir = backend_dir,
      force = force,
      strict = strict
    )
  }

  if (isTRUE(build_gene_gene)) {
    status$gene_gene <- .conseguiR_initialize_gene_gene_backend(
      paths = paths,
      backend_dir = backend_dir,
      force = force,
      strict = strict
    )
  }

  status
}

.conseguiR_initialize_backend_graphs <- function(
  backend_dir = NULL,
  build_gene_reg = TRUE,
  build_gene_gene = TRUE,
  force = FALSE,
  strict = TRUE,
  quiet = FALSE
) {
  context <- .conseguiR_backend_init_context(
    backend_dir = backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict
  )
  paths <- context$paths
  backend_dir <- context$backend_dir

  status <- .conseguiR_backend_init_status(
    paths = paths,
    backend_dir = backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict
  )

  if (!isTRUE(quiet)) {
    message("Backend graph initialization complete.")
  }

  .conseguiR_backend_init_result(
    backend_dir = backend_dir,
    paths = paths,
    status = status,
    force = force
  )
}
