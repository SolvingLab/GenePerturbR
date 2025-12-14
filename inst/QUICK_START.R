# ================================================================================
# GenePerturbR v1.1.0 - Quick Start Guide
# ================================================================================
# 5分钟快速上手指南
# ================================================================================

# ================================================================================
# 安装和设置
# ================================================================================

# 1. 安装包
# devtools::install("/path/to/GenePerturbR")

# 2. 安装可选依赖（用于网络图）
# install.packages(c("ggraph", "igraph"))

# 3. 设置数据库路径
Sys.setenv(APIKIT_DB_PATH = "/Users/liuzaoqu/Desktop/develop/APIKIT_Dev/API_DB")

# 4. 加载包
library(GenePerturbR)

# ================================================================================
# 基础查询（3个核心问题）
# ================================================================================

# 问题1: 敲除某个基因会发生什么？
result <- gpdb_what_happens("TP53")
cat(result$summary)
head(result$top_upregulated, 10)
head(result$top_downregulated, 10)

# 问题2: 谁调控MYC？
regulators <- gpdb_find_regulators("MYC", top_n = 20, min_confidence = "high")
head(regulators$repressors)   # 抑制MYC的基因
head(regulators$activators)   # 激活MYC的基因

# 问题3: METTL3调控什么？
targets <- gpdb_find_targets("METTL3", min_effect_size = 1.5, top_n = 30)
head(targets$upregulated)     # 上调的靶基因
head(targets$downregulated)   # 下调的靶基因

# ================================================================================
# 药物靶点预测
# ================================================================================

# 定义疾病特征（来自你的DEG分析或文献）
disease_signature <- data.frame(
  gene = c("MYC", "KRAS", "CDK4", "CCND1", "E2F1"),
  logFC = c(3.2, 2.8, 2.3, 2.9, 2.1)
)

# 寻找能逆转该特征的基因（潜在药物靶点）
candidates <- gpdb_predict_targets(disease_signature, mode = "reverse", top_n = 20)
print(head(candidates, 10))

# ================================================================================
# 可视化（全部优化，纯ggplot2！）
# ================================================================================

# 获取数据集
tp53_datasets <- gpdb_list_datasets(gene = "TP53")

# 1. 火山图（白色背景标签，统一主题）
p1 <- gpdb_plot_volcano(
  tp53_datasets$dataset_id[1],
  nlabel = 15
)
print(p1)

# 2. 多数据集热图（ggplot2引擎，纯CRAN！）
p2 <- gpdb_plot_heatmap(
  head(tp53_datasets$dataset_id, 5),
  top_n = 30,
  scale = "row"
)
print(p2)

# 3. 单数据集详细热图（带样本分组）
p3 <- gpdb_plot_heatmap_single(
  tp53_datasets$dataset_id[1],
  top_up = 20,
  top_down = 20,
  scale = "row"
)
print(p3)

# 4. 数据集分布图
p4 <- gpdb_plot_comparison("TP53", stratify_by = "tissue")
print(p4)

# 5. 🆕 调控网络图（ggraph）
p5 <- gpdb_plot_network(
  "MYC",
  top_regulators = 10,
  top_targets = 10,
  layout = "fr"  # force-directed layout
)
print(p5)

# 6. 🆕 调控级联分析
cascade <- gpdb_analyze_cascade("TP53", max_depth = 3)
print(cascade$n_paths)  # 发现的调控路径数

p6 <- gpdb_plot_cascade(cascade)
print(p6)

# ================================================================================
# 高级功能
# ================================================================================

# 比较多个基因
comparison <- gpdb_compare_genes(c("TP53", "RB1", "PTEN"))
print(paste("Common targets:", comparison$n_common))

# 基因家族分析
m6a_family <- gpdb_gene_family("m6A")
print(m6a_family$stats)

# 预测基因互作
interaction <- gpdb_predict_interaction("TP53", "MDM2")
print(interaction$prediction)
print(paste("Correlation:", round(interaction$correlation, 3)))

# ================================================================================
# 数据加载（优化版）
# ================================================================================

# 加载数据（自动转换基因ID为基因符号）
data <- gpdb_load_data(tp53_datasets$dataset_id[1], normalize = TRUE)

# Expression matrix 现在：
head(rownames(data$expression))  # 基因符号（不是ENSEMBL ID）!
print(data$info)  # 详细信息，包括样本分组

# 批量加载（带进度条）
batch_data <- gpdb_load_batch(
  head(tp53_datasets$dataset_id, 10),
  type = "deg",
  show_progress = TRUE
)

# ================================================================================
# 富集分析（NEW!）
# ================================================================================

cat("=== Enrichment Analysis ===\n")

# 从targets直接做富集
targets <- gpdb_find_targets("TP53", top_n = 100)
enrich_res <- gpdb_enrich(targets, enrich.type = "GO", GO.ont = "bp")

# 可视化：左边上调，右边下调
p_enrich <- gpdb_plot_enrichment(enrich_res, show.term.num = 15)
print(p_enrich)

# 使用KEGG数据库
enrich_kegg <- gpdb_enrich(targets, enrich.type = "KEGG")
p_kegg <- gpdb_plot_enrichment(enrich_kegg)
print(p_kegg)

# ================================================================================
# 下一步
# ================================================================================

cat("\n✨ 快速入门完成！\n\n")
cat("更多示例:\n")
cat("  - 完整教程: inst/tutorials/tutorial_comprehensive.R\n")
cat("  - 富集分析: inst/tutorials/tutorial_enrichment.R\n")
cat("  - 使用指南: inst/HOW_TO_USE.R\n")
cat("  - 性能测试: inst/tutorials/benchmark_performance.R\n")
cat("\n函数帮助:\n")
cat("  ?gpdb_what_happens\n")
cat("  ?gpdb_plot_network\n")
cat("  ?gpdb_enrich\n")

