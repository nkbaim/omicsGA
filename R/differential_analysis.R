#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-17 08:53:05
# @LastEditTime: 2025-01-17 14:03:40
# @LastEditors: Du Yang
# @Description: 差异分析模块
#

#' 执行差异表达分析
#'
#' @param expression_data 表达矩阵，行为基因，列为样本
#' @param group_info 分组信息因子
#' @param method 使用的方法 ("limma" 或 "deseq2")
#'
#' @return 差异分析结果数据框
#' @export
perform_differential_analysis <- function(expression_data, 
                                       group_info, 
                                       method = "limma",
                                       paired = FALSE,
                                       paired_info = NULL) {
    if (method == "limma") {
        # 创建设计矩阵
        design <- model.matrix(~group_info)
        
        # 执行limma分析
        fit <- limma::lmFit(expression_data, design)
        fit <- limma::eBayes(fit)
        
        # 获取结果
        results <- limma::topTable(fit, coef = 2, number = Inf)
        
    } else if (method == "deseq2") {
        # 创建DESeqDataSet对象
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = round(expression_data),
            colData = data.frame(condition = group_info),
            design = ~condition
        )
        
        # 执行DESeq2分析
        dds <- DESeq2::DESeq(dds)
        results <- DESeq2::results(dds)
        results <- as.data.frame(results)
    } else if (method == "wilcox") {
        results <- data.frame(
            row.names = rownames(expression_data),
            logFC = NA,
            pvalue = NA,
            padj = NA
        )
        
        for (gene in rownames(expression_data)) {
            expr_values <- expression_data[gene, ]
            if (paired) {
                test_result <- wilcox.test(
                    expr_values[group_info == levels(group_info)[1]],
                    expr_values[group_info == levels(group_info)[2]],
                    paired = TRUE
                )
            } else {
                test_result <- wilcox.test(
                    expr_values[group_info == levels(group_info)[1]],
                    expr_values[group_info == levels(group_info)[2]]
                )
            }
            
            # 计算log2 fold change
            mean_group1 <- mean(expr_values[group_info == levels(group_info)[1]])
            mean_group2 <- mean(expr_values[group_info == levels(group_info)[2]])
            results$logFC[gene] <- log2(mean_group2 / mean_group1)
            results$pvalue[gene] <- test_result$p.value
        }
        
        # 进行多重检验校正
        results$padj <- p.adjust(results$pvalue, method = "BH")
    }
    
    return(results)
}