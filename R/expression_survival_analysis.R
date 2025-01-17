#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-17 14:04:03
# @LastEditTime: 2025-01-17 14:04:10
# @LastEditors: Du Yang
# @Description: File content
#
perform_expression_survival_analysis <- function(expression_data, 
                                              survival_data,
                                              time_col,
                                              event_col,
                                              gene_list = NULL,
                                              cutoff_method = "median") {
    # 如果没有提供gene_list，使用所有基因
    if (is.null(gene_list)) {
        gene_list <- rownames(expression_data)
    }
    
    results <- list()
    
    for (gene in gene_list) {
        # 获取基因表达值
        expr_values <- expression_data[gene, ]
        
        # 根据表达值分组
        if (cutoff_method == "median") {
            cutoff <- median(expr_values)
            groups <- factor(ifelse(expr_values > cutoff, "High", "Low"))
        }
        
        # 合并生存数据和分组信息
        analysis_data <- data.frame(
            time = survival_data[[time_col]],
            event = survival_data[[event_col]],
            group = groups
        )
        
        # 执行生存分析
        surv_formula <- as.formula(paste("Surv(time, event) ~ group"))
        surv_fit <- survfit(surv_formula, data = analysis_data)
        
        # 执行log-rank检验
        log_rank <- survdiff(surv_formula, data = analysis_data)
        pvalue <- 1 - pchisq(log_rank$chisq, df = 1)
        
        # 保存结果
        results[[gene]] <- list(
            fit = surv_fit,
            pvalue = pvalue,
            data = analysis_data
        )
    }
    
    return(results)
} 