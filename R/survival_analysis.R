#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-17 08:53:01
# @LastEditTime: 2025-01-17 11:05:58
# @LastEditors: Du Yang
# @Description: 生存分析模块
#

#' 执行生存分析
#'
#' @param survival_data 包含生存数据的数据框
#' @param time_col 时间列名
#' @param event_col 事件列名
#' @param group_col 分组列名
#' @param covariates 协变量列名向量
#'
#' @return 生存分析结果列表，包含Cox模型和KM曲线
#' @export
perform_survival_analysis <- function(survival_data,
                                    time_col,
                                    event_col,
                                    group_col,
                                    covariates = NULL) {
    # 构建生存公式
    # formula_str <- paste0("survival::Surv(", time_col, ", ", event_col, ") ~ ", group_col)
    # if (!is.null(covariates)) {
    #    formula_str <- paste0(formula_str, " + ", paste(covariates, collapse = " + "))
    # }
    survival_data <- survival_data[, c(time_col, event_col, group_col)]
    colnames(survival_data) <- c("time", "status", "group")

    # 创建公式对象
    # surv_formula <- as.formula(formula_str)

    # 执行Cox回归
    cox_model <- survival::coxph(survival::Surv(time, status) ~ group, data = survival_data)

        # 创建生存对象用于KM曲线
    surv_fit <- survival::survfit(survival::Surv(time, status) ~ group, data = survival_data)

    # 生成KM曲线
    km_plot <- survminer::ggsurvplot(
        fit = surv_fit,  # 使用预先计算的survfit对象
        data = survival_data,  # 明确传递数据
        pval = TRUE,
        risk.table = TRUE,
        conf.int = FALSE
    )

    return(list(
        cox_model = cox_model,
        km_plot = km_plot
    ))
}
