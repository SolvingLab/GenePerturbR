# ================================================================================
# Gene Annotation Cache
# ================================================================================
# 启动时加载一次基因注释，存入包环境，避免重复查询
# Note: .gpdb_env is defined in zzz.R, shared across all modules
# ================================================================================

#' Get Gene Annotation Mapping (with cache)
#' @keywords internal
.gpdb_get_gene_mapping <- function() {
  # 检查缓存
  if (exists("gene_id_to_name", envir = .gpdb_env)) {
    return(.gpdb_env$gene_id_to_name)
  }

  # 第一次调用：加载并缓存
  con <- .gpdb_get_connection()
  gene_map <- DBI::dbGetQuery(con, "SELECT gene_id, gene_name FROM gene_annotation")

  # 创建快速查找向量
  mapping <- setNames(gene_map$gene_name, gene_map$gene_id)

  # 缓存到包环境
  .gpdb_env$gene_id_to_name <- mapping

  message("Loaded gene annotation mapping (", length(mapping), " genes, cached for session)")

  return(mapping)
}

#' Convert gene_id to gene_name in expression data (optimized)
#' @keywords internal
.gpdb_convert_to_gene_name <- function(expr_data) {
  if (!"gene_id" %in% names(expr_data)) {
    return(expr_data) # Already converted
  }

  # 使用缓存的映射（第一次会自动加载）
  gene_mapping <- .gpdb_get_gene_mapping()

  # 向量化转换（超快）
  gene_names <- gene_mapping[expr_data$gene_id]
  gene_names[is.na(gene_names)] <- expr_data$gene_id[is.na(gene_names)]

  # 移除 gene_id
  expr_data$gene_id <- NULL

  # 处理重复：简化版本（避免data.table复杂性）
  if (any(duplicated(gene_names))) {
    # 使用基础R聚合（够快了）
    expr_data$gene_name <- gene_names

    # 按gene_name聚合，取max
    aggregated <- stats::aggregate(
      . ~ gene_name,
      data = expr_data,
      FUN = max,
      na.rm = TRUE
    )

    rownames(aggregated) <- aggregated$gene_name
    aggregated$gene_name <- NULL
    expr_data <- aggregated
  } else {
    # 无重复，直接设置
    rownames(expr_data) <- gene_names
  }

  return(expr_data)
}
