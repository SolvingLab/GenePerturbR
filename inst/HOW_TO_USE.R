# ================================================================================
# GenePerturbR v1.1.0 - How to Use Guide
# ================================================================================
# 实际使用场景和最佳实践
# ================================================================================

Sys.setenv(APIKIT_DB_PATH = "/Users/liuzaoqu/Desktop/develop/APIKIT_Dev/API_DB")
library(GenePerturbR)

# ================================================================================
# 场景 1: 探索新基因的功能
# ================================================================================
# 研究问题: METTL3 是什么？它做什么？

cat("\n【场景1】探索新基因功能\n")
cat("════════════════════════════════════════\n")

# Step 1: 查询基因效应
mettl3_effects <- gpdb_what_happens("METTL3")
cat(mettl3_effects$summary, "\n")

# Step 2: 查看数据质量
cat("\n数据质量:\n")
cat("  数据集数量:", mettl3_effects$stats$n_datasets, "\n")
cat("  高置信度关系:", mettl3_effects$stats$n_high_confidence, "\n")

# Step 3: 查看主要靶基因
cat("\nTop 上调基因:\n")
print(head(mettl3_effects$top_upregulated, 5))

# Step 4: 可视化
mettl3_ds <- gpdb_list_datasets(gene = "METTL3")
if (nrow(mettl3_ds) > 0) {
  # 火山图
  gpdb_plot_volcano(mettl3_ds$dataset_id[1], nlabel = 15)
  
  # 分布图
  gpdb_plot_comparison("METTL3", stratify_by = "tissue")
}

# ================================================================================
# 场景 2: 寻找疾病治疗靶点
# ================================================================================
# 研究问题: 如何找到逆转癌症特征的基因？

cat("\n【场景2】药物靶点预测\n")
cat("════════════════════════════════════════\n")

# Step 1: 定义疾病特征（从DEG分析或文献）
cancer_signature <- data.frame(
  gene = c("MYC", "KRAS", "BRAF", "CDK4", "CCND1", "E2F1", "MCM2", "PCNA"),
  logFC = c(3.5, 2.8, 2.5, 2.3, 2.9, 2.1, 2.4, 2.6)
)

cat("\n癌症特征基因: ", paste(cancer_signature$gene, collapse = ", "), "\n")

# Step 2: 寻找逆转该特征的基因
candidates <- gpdb_predict_targets(
  cancer_signature,
  mode = "reverse",
  top_n = 30,
  min_confidence = "medium"
)

cat("\nTop 10 候选靶点:\n")
print(head(candidates[, c("perturbed_gene", "total_score", "n_signature_matches", "match_rate")], 10))

# Step 3: 深入研究Top候选
top_gene <- candidates$perturbed_gene[1]
cat("\n深入分析Top候选:", top_gene, "\n")

# 查看该基因的效应
top_validation <- gpdb_what_happens(top_gene)
cat(top_validation$summary, "\n")

# 可视化调控网络
gpdb_plot_network(top_gene, top_regulators = 8, top_targets = 8)

# Step 4: 检查在相关组织的可用性
top_datasets <- gpdb_list_datasets(gene = top_gene, tissue = "Liver")
cat("\n在Liver组织的数据集数:", nrow(top_datasets), "\n")

# ================================================================================
# 场景 3: 构建调控网络
# ================================================================================
# 研究问题: MYC的完整调控网络是什么样的？

cat("\n【场景3】调控网络分析\n")
cat("════════════════════════════════════════\n")

# Step 1: 找谁调控MYC
myc_regulators <- gpdb_find_regulators("MYC", top_n = 15, min_confidence = "high")

cat("\nMYC 调控因子:\n")
cat("  Repressors (抑制MYC):", nrow(myc_regulators$repressors), "个\n")
cat("  Activators (激活MYC):", nrow(myc_regulators$activators), "个\n")

# Step 2: MYC调控什么
myc_targets <- gpdb_find_targets("MYC", min_effect_size = 1.0, top_n = 15)

cat("\nMYC 靶基因:\n")
cat("  Upregulated:", nrow(myc_targets$upregulated), "个\n")
cat("  Downregulated:", nrow(myc_targets$downregulated), "个\n")

# Step 3: 🆕 可视化完整网络
p_network <- gpdb_plot_network(
  "MYC",
  top_regulators = 12,
  top_targets = 12,
  layout = "fr",
  node_size = 10
)
print(p_network)

# Step 4: 🆕 追踪调控级联
cascade <- gpdb_analyze_cascade(
  "MYC",
  max_depth = 3,
  min_effect_size = 1.2
)

cat("\n级联分析:\n")
cat("  发现路径:", cascade$n_paths, "条\n")
cat("  涉及基因:", cascade$n_genes, "个\n")

# 可视化级联
if (cascade$n_paths > 0) {
  gpdb_plot_cascade(cascade)
}

# ================================================================================
# 场景 4: 比较基因家族
# ================================================================================
# 研究问题: m6A修饰酶的功能如何比较？

cat("\n【场景4】基因家族比较\n")
cat("════════════════════════════════════════\n")

# Step 1: 定义基因家族
m6a_genes <- c("METTL3", "METTL14", "WTAP", "ALKBH5", "FTO")

# Step 2: 比较它们的靶基因
comparison <- gpdb_compare_genes(m6a_genes)

cat("\n家族比较:\n")
cat("  共同靶基因:", comparison$n_common, "个\n")
cat("  METTL3 特异:", length(comparison$unique_targets$METTL3), "个\n")
cat("  ALKBH5 特异:", length(comparison$unique_targets$ALKBH5), "个\n")

# Step 3: 或使用预定义家族
family_analysis <- gpdb_gene_family("m6A")
cat("\n家族成员:", paste(family_analysis$members, collapse = ", "), "\n")

# Step 4: 可视化比较（热图）
all_datasets <- list()
for (gene in m6a_genes[1:3]) {
  ds <- gpdb_list_datasets(gene = gene)
  if (nrow(ds) > 0) {
    all_datasets[[gene]] <- head(ds$dataset_id, 2)
  }
}

all_ids <- unlist(all_datasets)
if (length(all_ids) >= 2) {
  gpdb_plot_heatmap(all_ids, top_n = 30, scale = "row")
}

# ================================================================================
# 场景 5: 组织特异性分析
# ================================================================================
# 研究问题: TP53在不同组织的效应有何不同？

cat("\n【场景5】组织特异性分析\n")
cat("════════════════════════════════════════\n")

# Step 1: 比较不同组织
tissue_comparison <- gpdb_compare_contexts(
  "TP53",
  contexts = list(
    liver = list(tissue = "Liver"),
    lung = list(tissue = "Lung"),
    brain = list(tissue = "Brain")
  )
)

cat("\n组织比较:\n")
cat("  共同靶基因:", length(tissue_comparison$common_targets), "个\n")
cat("  Liver特异:", length(tissue_comparison$unique_targets$liver), "个\n")
cat("  Lung特异:", length(tissue_comparison$unique_targets$lung), "个\n")

# Step 2: 可视化分布
gpdb_plot_comparison("TP53", stratify_by = "tissue")

# ================================================================================
# 场景 6: 加载和分析原始数据
# ================================================================================
# 研究问题: 我想用原始数据做自定义分析

cat("\n【场景6】原始数据分析\n")
cat("════════════════════════════════════════\n")

# Step 1: 加载完整数据
data <- gpdb_load_data(tp53_datasets$dataset_id[1], normalize = TRUE)

# Step 2: 检查数据结构
cat("\nExpression matrix:\n")
cat("  Rownames (前5个):", paste(head(rownames(data$expression), 5), collapse = ", "), "\n")
cat("  维度:", paste(dim(data$expression), collapse = " x "), "\n")

cat("\nMetadata:\n")
print(head(data$metadata))
cat("  分组:", data$info$sample_groups, "\n")

# Step 3: Expression matrix 可直接用于其他包
# 例如 DESeq2, edgeR, limma 等
# rownames 已经是基因符号，无需转换！

# ================================================================================
# 场景 7: 🆕 批量分析流程（优化版）
# ================================================================================

cat("\n【场景7】批量分析（使用优化功能）\n")
cat("════════════════════════════════════════\n")

# Step 1: 批量加载（7x faster！）
batch_deg <- gpdb_load_batch(
  head(tp53_datasets$dataset_id, 10),
  type = "deg",
  show_progress = TRUE  # 显示进度条
)

# Step 2: 聚合分析
all_deg_data <- do.call(rbind, lapply(names(batch_deg), function(id) {
  deg <- batch_deg[[id]]
  deg$dataset_id <- id
  deg[!is.na(deg$adj.P.Val) & deg$adj.P.Val < 0.05, ]
}))

# Step 3: 找一致性基因
gene_freq <- table(all_deg_data$gene)
consistent <- names(gene_freq[gene_freq >= 5])
cat("  一致性显著基因:", length(consistent), "个\n")

# ================================================================================
# 💡 使用技巧
# ================================================================================

cat("\n💡 使用技巧\n")
cat("════════════════════════════════════════\n")
cat("
1. 🔍 查询优化:
   - 使用 min_confidence='high' 获取可靠结果
   - 使用 min_effect_size 过滤弱效应
   - 使用 top_n 控制结果数量

2. 🎨 可视化:
   - 所有函数支持自定义 theme 参数
   - 使用 colors=NULL 获取默认美观配色
   - 热图支持 show_values=TRUE 显示数值

3. 🚀 性能:
   - 批量操作使用 gpdb_load_batch()
   - 基因注释自动缓存，第二次更快
   - show_progress=FALSE 可关闭进度条

4. 🆕 新功能:
   - gpdb_plot_network() - 调控网络
   - gpdb_analyze_cascade() - 多层级联
   - 所有热图使用ggplot2（纯CRAN）

5. 📊 数据格式:
   - Expression matrix rownames = 基因符号
   - Metadata 中 'group' 列 = treatment/control
   - 完全兼容 DESeq2, edgeR, limma
")

cat("\n✓ 使用指南完成！\n")

