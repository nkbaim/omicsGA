#
# @Copyright: Peking University Cancer Hospital, All Rights Reserved.
# @Author: Du Yang
# @Date: 2025-01-17 14:04:22
# @LastEditTime: 2025-01-17 14:04:28
# @LastEditors: Du Yang
# @Description: File content
#
create_output_directory <- function(analysis_type) {
    # 创建基础输出目录
    base_dir <- file.path("results", analysis_type)
    
    # 创建时间戳子目录
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_dir <- file.path(base_dir, timestamp)
    
    # 创建必要的子目录
    dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    
    return(output_dir)
} 