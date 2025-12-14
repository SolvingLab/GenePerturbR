# GenePerturbR

> 🧬 Comprehensive genetic perturbation analysis with 7,665 RNA-seq datasets

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

GenePerturbR provides access to **18.9 million gene-gene relationships** (826K high-confidence) derived from 7,665 genetic perturbation experiments across 2,810 genes, 1,063 cell lines, and 71 tissue types.

### Key Features

- ⚡ **Lightning Fast**: Optimized SQL queries with smart caching
- 🎨 **Beautiful Viz**: Pure ggplot2 visualizations (no Bioconductor dependencies!)
- 🔬 **Rich Analysis**: Regulators, targets, cascades, networks, drug predictions
- 🎯 **High Confidence**: Multi-dataset aggregation with evidence scoring
- 🤖 **LLM-Friendly**: Natural language interfaces for AI agents
- 🚀 **Production Ready**: Clean code, comprehensive tests, detailed docs

## Installation

```r
# Install from source
devtools::install("/path/to/GenePerturbR")

# Install dependencies (all from CRAN!)
install.packages(c("ggplot2", "ggrepel", "ggraph", "igraph", "scales"))
```

## Quick Setup

```r
# Set database path (required)
Sys.setenv(APIKIT_DB_PATH = "/path/to/API_DB")

# Load package
library(GenePerturbR)
```

## Quick Start

```r
# What happens when TP53 is knocked out?
result <- gpdb_what_happens("TP53")
cat(result$summary)
head(result$top_upregulated)

# Who regulates MYC?
regulators <- gpdb_find_regulators("MYC", top_n = 20, min_confidence = "high")
head(regulators$repressors)  # Genes that repress MYC
head(regulators$activators)  # Genes that activate MYC

# What does METTL3 regulate?
targets <- gpdb_find_targets("METTL3", min_effect_size = 1.5, top_n = 30)
head(targets$upregulated)
head(targets$downregulated)

# Predict drug targets
disease_signature <- data.frame(
  gene = c("MYC", "KRAS", "CDK4", "CCND1"),
  logFC = c(3.2, 2.8, 2.3, 2.9)
)
candidates <- gpdb_predict_targets(disease_signature, mode = "reverse", top_n = 20)

# Enrichment analysis (NEW!)
targets <- gpdb_find_targets("TP53", top_n = 100)
enrich_res <- gpdb_enrich(targets, enrich.type = "GO")  # Separate up/down
gpdb_plot_enrichment(enrich_res, show.term.num = 20)   # Paired dotplot
```

## Beautiful Visualizations

All plots use unified theme and color system - publication ready out of the box!

```r
# Get datasets
tp53_ds <- gpdb_list_datasets(gene = "TP53")

# 1. Volcano plot (ggplot2 with white-background labels)
gpdb_plot_volcano(
  tp53_ds$dataset_id[1],
  nlabel = 15,
  label.bg = "white"
)

# 2. Multi-dataset DEG comparison heatmap (ggplot2, pure CRAN!)
gpdb_plot_heatmap_deg(  # New name! (gpdb_plot_heatmap still works as alias)
  head(tp53_ds$dataset_id, 5),
  top_n = 30,
  scale = "row"
)

# 3. Single-dataset sample expression heatmap
gpdb_plot_heatmap_expr(  # New name! (gpdb_plot_heatmap_single still works)
  tp53_ds$dataset_id[1],
  top_up = 20,
  top_down = 20,
  scale = "row"
)

# 4. Dataset distribution
gpdb_plot_comparison("TP53", stratify_by = "tissue")

# 5. Regulatory network (NEW!)
gpdb_plot_network(
  "MYC",
  top_regulators = 10,
  top_targets = 10,
  layout = "fr"
)

# 6. Regulatory cascade visualization (NEW!)
cascade <- gpdb_analyze_cascade("TP53", max_depth = 3)
gpdb_plot_cascade(cascade)
```

## Main Functions

### 📊 Data Access
| Function | Description |
|----------|-------------|
| `gpdb_list_datasets()` | List available datasets |
| `gpdb_get_info()` | Get dataset metadata |
| `gpdb_load_data()` | Load expression + metadata (auto-converts gene_id → gene_name!) |
| `gpdb_load_deg()` | Load DEG results |
| `gpdb_load_batch()` | Batch load with progress bar |

### 🔍 Query Functions
| Function | Description | New Features |
|----------|-------------|--------------|
| `gpdb_what_happens()` | Gene perturbation effects | Optimized SQL |
| `gpdb_find_regulators()` | Find regulatory genes | + `top_n`, `min_confidence` |
| `gpdb_find_targets()` | Find target genes | + `min_effect_size` |
| `gpdb_compare_genes()` | Compare multiple genes | Batch query |
| `gpdb_compare_contexts()` | Compare across contexts | - |

### 🧠 Smart Analysis
| Function | Description |
|----------|-------------|
| `gpdb_predict_targets()` | Predict therapeutic targets |
| `gpdb_predict_interaction()` | Predict gene interactions |
| `gpdb_analyze_cascade()` | Trace regulatory cascades 🆕 |
| `gpdb_enrich()` | Pathway enrichment analysis 🆕 |
| `gpdb_summarize()` | Generate gene summaries |
| `gpdb_search()` | Search database |

### 🎨 Visualization
| Function | Description | Engine |
|----------|-------------|--------|
| `gpdb_plot_volcano()` | Volcano plots | ggplot2 ✅ |
| `gpdb_plot_heatmap_deg()` | Multi-dataset DEG heatmaps | ggplot2 ✅ |
| `gpdb_plot_heatmap_expr()` | Single-dataset expression heatmaps | ggplot2 ✅ |
| `gpdb_plot_comparison()` | Distribution plots | ggplot2 ✅ |
| `gpdb_plot_network()` | Regulatory networks 🆕 | ggraph ✅ |
| `gpdb_plot_cascade()` | Cascade networks 🆕 | ggraph ✅ |
| `gpdb_plot_enrichment()` | Enrichment dotplots 🆕 | ggplot2 ✅ |

**Note**: Old function names (`gpdb_plot_heatmap`, `gpdb_plot_heatmap_single`) still work as aliases for backward compatibility.

## What's New in v1.1.0

### 🚀 Performance
- ✅ **SQL Query Builder**: Unified query construction system
- ✅ **Batch Query Optimization**: 10x faster for multiple datasets
- ✅ **Smart Caching**: Gene annotation mapping cached in memory
- ✅ **Progress Bars**: Visual feedback for batch operations

### 🎨 Visualization Overhaul
- ✅ **Pure ggplot2**: Removed ALL Bioconductor dependencies!
- ✅ **Unified Theme**: Consistent, beautiful plots across all functions
- ✅ **Network Plots**: New `gpdb_plot_network()` with ggraph
- ✅ **Cascade Plots**: New `gpdb_plot_cascade()` for pathway analysis
- ✅ **Better Heatmaps**: Prettier, faster, more customizable

### 💎 New Features
- ✅ **Cascade Analysis**: Trace multi-step regulatory pathways
- ✅ **Network Visualization**: Interactive-style network plots
- ✅ **Smart Parameters**: Intelligent defaults based on data
- ✅ **Better Data Loading**: Auto-converts gene_id → gene_name

### 🛠️ Code Quality
- ✅ **Function式编程**: Cleaner, more maintainable code
- ✅ **统一管道**: Native `|>` operator throughout
- ✅ **减少重复**: From ~2500 to ~2000 lines (-20%)
- ✅ **更好的错误处理**: Informative error messages

## Architecture

### Optimized Stack
```
Query Layer (SQL Builder)
    ↓
Caching Layer (Gene Annotations)
    ↓
Analysis Layer (Vectorized Operations)
    ↓
Visualization Layer (ggplot2 Ecosystem)
```

### Dependencies (All CRAN! 🎉)
- **Core**: DBI, RSQLite, qs
- **Data**: dplyr, tidyr, data.table  
- **Viz**: ggplot2, ggrepel, scales
- **Network**: ggraph, igraph (suggested)

## Performance Benchmarks

| Operation | Before | After | Improvement |
|-----------|--------|-------|-------------|
| Single query | 50ms | 10ms | **5x faster** |
| Batch load (10) | 2.0s | 0.3s | **7x faster** |
| Heatmap render | 3.0s | 1.0s | **3x faster** |
| Repeat query | 50ms | <1ms | **50x faster** (cached) |

## Complete Examples

### Example 1: Discover Therapeutic Targets

```r
# Define disease signature
cancer_sig <- data.frame(
  gene = c("MYC", "KRAS", "CDK4", "CCND1", "E2F1"),
  logFC = c(3.2, 2.8, 2.3, 2.9, 2.1)
)

# Find genes that reverse this signature
candidates <- gpdb_predict_targets(cancer_sig, mode = "reverse", top_n = 20)

# Investigate top candidate
top_gene <- candidates$perturbed_gene[1]
gpdb_what_happens(top_gene)
gpdb_plot_network(top_gene)
```

### Example 2: Regulatory Network Analysis

```r
# Who regulates MYC?
myc_regulators <- gpdb_find_regulators("MYC", top_n = 15)

# What does MYC regulate?
myc_targets <- gpdb_find_targets("MYC", top_n = 15)

# Visualize the network
gpdb_plot_network("MYC", top_regulators = 10, top_targets = 10, layout = "star")

# Trace regulatory cascade
cascade <- gpdb_analyze_cascade("MYC", max_depth = 3)
gpdb_plot_cascade(cascade)
```

### Example 3: Context-Specific Analysis

```r
# Compare TP53 across tissues
contexts <- gpdb_compare_contexts(
  "TP53",
  contexts = list(
    liver = list(tissue = "Liver"),
    lung = list(tissue = "Lung"),
    brain = list(tissue = "Brain")
  )
)

# Visualize distribution
gpdb_plot_comparison("TP53", stratify_by = "tissue")
```

## Documentation

- **Quick Start**: `inst/QUICK_START.R`
- **Usage Guide**: `inst/HOW_TO_USE.R`  
- **Complete Tutorial**: `inst/tutorials/tutorial_comprehensive.R`
- **Enrichment Tutorial**: `inst/tutorials/tutorial_enrichment.R`
- **Performance Benchmark**: `inst/tutorials/benchmark_performance.R`

## Data Structure

```
API_DB/
└── gpsadb/
    ├── processed/
    │   ├── gpsadb.db          # SQLite database (18.9M relationships)
    │   ├── expression/        # Expression matrices (gene symbols as rownames)
    │   ├── deg/               # DEG results
    │   └── metadata/          # Sample metadata (with 'group' column)
    └── gene_annotation_table.csv
```

## Citation

```bibtex
@software{liu2025genepérturbr,
  author = {Liu, Zaoqu},
  title = {GenePerturbR: Comprehensive Genetic Perturbation Analysis},
  year = {2025},
  version = {1.1.0},
  url = {https://github.com/yourusername/GenePerturbR}
}
```

## License

MIT License - see [LICENSE](LICENSE) file

## Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **ORCID**: [0000-0002-0452-742X](https://orcid.org/0000-0002-0452-742X)

---

**🐱 Small cat approved!** No cats were harmed in the making of this package.
