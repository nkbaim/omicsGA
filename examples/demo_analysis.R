#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-17 08:05:01
# @LastEditTime: 2025-01-17 08:23:29
# @LastEditors: Du Yang
# @Description: 分析流程演示脚本
#

library(omicsCompare)

# 1. 生成模拟数据 ----
set.seed(123)

# 生成表达矩阵
n_genes <- 1000  # 基因数
n_samples <- 100 # 样本数
expression_data <- matrix(
    rnorm(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples
)
rownames(expression_data) <- paste0("gene", 1:n_genes)
colnames(expression_data) <- paste0("sample", 1:n_samples)

# 生成临床数据
clinical_data <- data.frame(
    sample_id = paste0("sample", 1:n_samples),
    group = factor(rep(c("Control", "Treatment"), each = n_samples/2)),
    age = rnorm(n_samples, mean = 60, sd = 10),
    stage = factor(sample(c("I", "II", "III"), n_samples, replace = TRUE)),
    survival_time = rexp(n_samples, rate = 0.1),
    survival_status = rbinom(n_samples, 1, 0.3)
)

# 生成基因集
gene_sets <- list(
    pathway1 = sample(rownames(expression_data), 50),
    pathway2 = sample(rownames(expression_data), 50),
    pathway3 = sample(rownames(expression_data), 50)
)

# 2. 准备数据列表 ----
data_list <- list(
    expression_data = expression_data,
    clinical_data = clinical_data,
    time_col = "survival_time",
    event_col = "survival_status",
    group_col = "group",
    group_info = clinical_data$group,
    gene_sets = gene_sets
)

# 3. 准备分析配置 ----
config <- list(
    survival = list(
        covariates = c("age", "stage")
    ),
    differential = list(
        method = "limma"
    )
)

# 4. 创建分析流程 ----
# 首先创建基础流程
pipeline <- create_default_pipeline()

# 添加GSEA模块
GSEAModule <- R6::R6Class("GSEAModule",
    inherit = AnalysisModule,
    public = list(
        initialize = function() {
            super$initialize(
                name = "gsea",
                description = "GSEA富集分析模块",
                required_data = c("expression_data", "group_info", "gene_sets")
            )
        },
        
        execute = function(data_list, params) {
            super$validate_input(data_list)
            
            # 计算每个基因的差异统计量
            diff_stats <- apply(data_list$expression_data, 1, function(x) {
                t.test(x ~ data_list$group_info)$statistic
            })
            
            return(perform_gsea(
                ranked_genes = diff_stats,
                pathway_sets = data_list$gene_sets
            ))
        }
    )
)

# 注册GSEA模块
pipeline$register_module(GSEAModule$new())

# 更新配置以包含GSEA
config$gsea <- list()  # GSEA模块不需要额外参数

# 5. 执行分析 ----
results <- pipeline$execute(data_list, config, "demo_results")

# 6. 查看结果 ----
# 生存分析结果
print("生存分析结果:")
print(summary(results$survival$cox_model))

# 差异分析结果
print("差异分析结果前几行:")
print(head(results$differential))

# GSEA结果
print("GSEA分析结果:")
print(head(results$gsea))

# 7. 结果可视化示例 ----
# 生存曲线
pdf("demo_results/survival_curve.pdf")
print(results$survival$km_plot)
dev.off()

# 差异表达火山图
volcano_plot <- ggplot2::ggplot(
    as.data.frame(results$differential),
    ggplot2::aes(x = logFC, y = -log10(adj.P.Val))
) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Volcano Plot")

pdf("demo_results/volcano_plot.pdf")
print(volcano_plot)
dev.off() 