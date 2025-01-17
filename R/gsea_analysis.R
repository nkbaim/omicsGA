#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-17 08:53:08
# @LastEditTime: 2025-01-17 10:35:19
# @LastEditors: Du Yang
# @Description: 富集分析模块
#

#' 执行GSEA分析
#'
#' @param ranked_genes 基因排序列表，名称为基因名，值为排序统计量
#' @param pathway_sets 通路基因集列表
#' @param min_size 最小基因集大小
#' @param max_size 最大基因集大小
#'
#' @return GSEA分析结果数据框
#' @export
perform_gsea <- function(ranked_genes, 
                        pathway_sets, 
                        min_size = 15, 
                        max_size = 500) {
    # 执行fgsea分析
    gsea_results <- fgsea::fgsea(
        pathways = pathway_sets,
        stats = ranked_genes,
        minSize = min_size,
        maxSize = max_size
    )
    
    return(as.data.frame(gsea_results))
}

#' 执行GSVA分析
#'
#' @param expression_data 表达矩阵，行为基因，列为样本
#' @param gene_sets 基因集列表
#' @param method GSVA方法，默认为"gsva"
#'
#' @return GSVA分析结果矩阵
#' @export
perform_gsva <- function(expression_data, 
                        gene_sets, 
                        method = "gsva") {
    # 执行GSVA分析
    gsva_results <- GSVA::gsva(
        expr = as.matrix(expression_data),
        gset.idx.list = gene_sets,
        method = method,
        parallel.sz = 1
    )
    
    return(gsva_results)
}