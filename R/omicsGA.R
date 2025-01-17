#' 分析模块基类
#' @export
AnalysisModule <- R6::R6Class("AnalysisModule",
                              public = list(
                                #' @field name 模块名称
                                name = NULL,
                                #' @field description 模块描述
                                description = NULL,
                                #' @field required_data 所需数据列表
                                required_data = NULL,

                                #' @description 初始化函数
                                #' @param name 模块名称
                                #' @param description 模块描述
                                #' @param required_data 所需数据列表
                                initialize = function(name, description, required_data) {
                                  self$name <- name
                                  self$description <- description
                                  self$required_data <- required_data
                                },

                                #' @description 验证输入数据
                                #' @param data_list 输入数据列表
                                validate_input = function(data_list) {
                                  missing_data <- setdiff(self$required_data, names(data_list))
                                  if (length(missing_data) > 0) {
                                    stop(sprintf("缺少所需数据: %s", paste(missing_data, collapse = ", ")))
                                  }
                                  return(TRUE)
                                },

                                #' @description 执行分析
                                #' @param data_list 输入数据列表
                                #' @param params 分析参数
                                execute = function(data_list, params) {
                                  stop("需要在子类中实现execute方法")
                                }
                              )
)

#' 生存分析模块
#' @export
SurvivalModule <- R6::R6Class("SurvivalModule",
                              inherit = AnalysisModule,
                              public = list(
                                initialize = function() {
                                  super$initialize(
                                    name = "survival",
                                    description = "生存分析模块",
                                    required_data = c("clinical_data", "time_col", "event_col", "group_col")
                                  )
                                },

                                execute = function(data_list, params) {
                                  super$validate_input(data_list)
                                  print(data_list$clinical_data)
                                  return(perform_survival_analysis(
                                    survival_data = data_list$clinical_data,
                                    time_col = data_list$time_col,
                                    event_col = data_list$event_col,
                                    group_col = data_list$group_col,
                                    covariates = params$covariates
                                  ))
                                }
                              )
)

#' 差异分析模块
#' @export
DifferentialModule <- R6::R6Class("DifferentialModule",
                                  inherit = AnalysisModule,
                                  public = list(
                                    initialize = function() {
                                      super$initialize(
                                        name = "differential",
                                        description = "差异表达分析模块",
                                        required_data = c("expression_data", "group_info")
                                      )
                                    },

                                    execute = function(data_list, params) {
                                      super$validate_input(data_list)
                                      return(perform_differential_analysis(
                                        expression_data = data_list$expression_data,
                                        group_info = data_list$group_info,
                                        method = params$method
                                      ))
                                    }
                                  )
)

#' 分析流程管理器
#' @export
AnalysisPipeline <- R6::R6Class("AnalysisPipeline",
                                public = list(
                                  #' @field modules 已注册的分析模块
                                  modules = list(),
                                  #' @field results 分析结果
                                  results = list(),

                                  #' @description 注册分析模块
                                  #' @param module 分析模块对象
                                  register_module = function(module) {
                                    self$modules[[module$name]] <- module
                                    invisible(self)
                                  },

                                  #' @description 执行分析流程
                                  #' @param data_list 输入数据列表
                                  #' @param config 分析配置列表
                                  #' @param output_dir 输出目录
                                  execute = function(data_list, config, output_dir = "results") {
                                    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

                                    for (module_name in names(config)) {
                                      if (!module_name %in% names(self$modules)) {
                                        warning(sprintf("未找到模块: %s，跳过", module_name))
                                        next
                                      }

                                      message(sprintf("执行模块: %s", module_name))
                                      module <- self$modules[[module_name]]

                                      tryCatch({
                                        result <- module$execute(data_list, config[[module_name]])
                                        self$results[[module_name]] <- result
                                        self$save_result(result, module_name, output_dir)
                                      }, error = function(e) {
                                        warning(sprintf("模块 %s 执行失败: %s", module_name, e$message))
                                      })
                                    }

                                    return(self$results)
                                  },

                                  #' @description 保存分析结果
                                  #' @param result 分析结果
                                  #' @param module_name 模块名称
                                  #' @param output_dir 输出目录
                                  save_result = function(result, module_name, output_dir) {
                                    result_file <- file.path(output_dir, sprintf("%s_results.rds", module_name))
                                    saveRDS(result, result_file)
                                  }
                                )
)

#' 创建默认分析流程
#' @export
create_default_pipeline <- function() {
  pipeline <- AnalysisPipeline$new()
  pipeline$register_module(SurvivalModule$new())
  pipeline$register_module(DifferentialModule$new())
  return(pipeline)
}
