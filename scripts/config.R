# ============================================================================
# config.R
# Sargentodoxa SDM 项目的全局配置文件
#
# 物种：大血藤 (Sargentodoxa cuneata)
# 主要流程：
#   1. 出现点清洗与空间稀疏化
#   2. 伪缺失点生成
#   3. 环境变量筛选
#   4. 空间交叉验证
#   5. 单模型训练、筛选与集成
#   6. 现代 / 古气候投影与可视化
#
# 使用建议：
#   - 优先修改参数值，不要随意改动下方辅助函数。
#   - 参数分区按分析流程排列，建议按区块顺序检查。
#   - 若切换服务器目录、古气候数据源或投影版本，优先检查路径配置。
#   - 每个分区内进一步分为：
#       常用：第一次跑流程时最可能需要改
#       高级：有明确建模目的时再调整
#       一般不用动：除非排错或重构流程，否则保持默认
# ============================================================================

# ==================== 1. 路径配置 ====================
# 本仓库已去服务器化：
#   - 本地输入数据（出现点 / 化石记录 / TEYF 时期元数据）随仓库提供，走相对路径；
#   - 服务器端气候栅格（WorldClim / TEYF）通过环境变量配置，默认值保留服务器原始
#     路径供参考；本地不随仓库提供。
#   - 出图脚本 scripts/visualization/ 自带 project_root 推断，不依赖本节路径；
#     本节主要供管线脚本 01-09 在服务器端重跑时使用。
# 约定：管线脚本从仓库根目录运行，例如 Rscript scripts/01_data_preparation.R。

# 推断仓库根目录（兼容 source() 与 Rscript，且不依赖 cwd）。
# config.R 位于 scripts/config.R，故根为其上一级。
project_dir <- tryCatch({
  ofiles <- vapply(sys.frames(), function(f) f$ofile, character(1))
  ofile <- ofiles[length(ofiles)]
  if (!is.na(ofile) && nzchar(ofile)) {
    normalizePath(file.path(dirname(ofile), ".."), mustWork = FALSE)
  } else {
    args <- commandArgs(trailingOnly = FALSE)
    fa <- args[grep("^--file=", args)]
    if (length(fa)) {
      normalizePath(file.path(dirname(sub("^--file=", "", fa[1])), ".."), mustWork = FALSE)
    } else {
      normalizePath(".", mustWork = FALSE)
    }
  }
}, error = function(e) normalizePath(".", mustWork = FALSE))

data_dir <- file.path(project_dir, "data")

# 常用：出现点数据（CVH / GBIF）
# 如果要切换数据来源，通常只需要改后面的 sp_file。
sp_file_cvh <- file.path(project_dir, "data/occ_data_cvh/sargentodoxa_occ_cvh_counties_filtered.csv")
sp_file_gbif <- file.path(project_dir, "data/occ_data_cvh/sargentodoxa_occ_filtered.csv")

# 一般不用动：现代气候数据目录（WorldClim 2.1，5 arc-min）
# 服务器：/home/cedrus/data/shz/data/wc2.1_5m_bio
# 本地未随仓库提供；如需重跑 01-08，设置环境变量 SARGENTO_MODERN_CLIM_DIR 指向本地 WorldClim 目录。
modern_clim_dir <- Sys.getenv(
  "SARGENTO_MODERN_CLIM_DIR",
  unset = "/home/cedrus/data/shz/data/wc2.1_5m_bio"
)
modern_clim_pattern <- "*.tif$"

# 常用：古气候数据来源
# 如果只是在两套古气候数据之间切换，通常只改 current_past_clim 即可。
# 1. TEWT 0-3600 ka（服务器：/home/cedrus/data/shz/data/tewt_0_3600ka/masked_landsea）
past_clim_dir_tewt <- Sys.getenv(
  "SARGENTO_TEWT_CLIM_DIR",
  unset = "/home/cedrus/data/shz/data/tewt_0_3600ka/masked_landsea"
)
past_clim_subdir_tewt <- TRUE  # TRUE = 每个时期一个子目录

# 2. TEYF-mask-nc-5min（服务器：/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min）
past_clim_dir_teyf <- Sys.getenv(
  "SARGENTO_TEYF_CLIM_DIR",
  unset = "/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min"
)
past_clim_subdir_teyf <- TRUE  # TRUE = 每个时期一个子目录

# 常用：当前启用的古气候数据源
# 可选："tewt" / "teyf"
current_past_clim <- "teyf"

# 根据 current_past_clim 自动选择实际使用的古气候目录
if (current_past_clim == "tewt") {
  past_clim_dir <- past_clim_dir_tewt
  past_clim_subdir <- past_clim_subdir_tewt
} else {
  past_clim_dir <- past_clim_dir_teyf
  past_clim_subdir <- past_clim_subdir_teyf
}

# 高级：陆地边界（可选，仅用于可视化）
# 没有 shapefile 时可设为 NULL
land_file <- NULL

# 常用：化石记录
# 主要用于古气候投影结果与化石证据的对照展示。
fossil_dir <- file.path(project_dir, "data/fossil_record")
# fossil_record_file 至少需要包含：paleoLong, paleoLat, age, simulation_name
fossil_record_file <- file.path(fossil_dir, "FossilRecord-paleocor-v2.1.csv")

# 高级：古气候时期顺序信息文件
# 随仓库提供（documents/teyf_simulations.csv）。
teyf_simulations_file <- file.path(project_dir, "documents/teyf_simulations.csv")
simulations_file <- teyf_simulations_file

# ==================== 2. 物种与输入数据 ====================

species_name <- "Sargentodoxa"
resp_name <- "Sargentodoxa"

# 常用：默认使用的出现点数据源
# 第一次跑流程时，这通常是最先需要确认的输入文件。
sp_file <- sp_file_cvh

# ==================== 3. 兼容性别名 ====================
# 这一段主要给旧脚本兼容使用；新脚本若统一用 modern_clim_dir / past_clim_dir，
# 通常不需要再改这里。

# 一般不用动：现代气候目录别名
wc2_clim_dir <- modern_clim_dir

# 一般不用动：古气候目录别名
tewt_dir <- past_clim_dir

# ==================== 4. 输出目录 ====================

output_base <- "./output"
dir.create(output_base, showWarnings = FALSE)

# 一般不用动：各流程步骤的输出目录
# 如果你只是正常运行流程，通常不需要改；除非想调整项目输出结构。
output_dirs <- list(
  data_prep = "./output/01_data_prep",
  pseudo_absence = "./output/02_pseudo_absence",
  variable_select = "./output/03_variable_select",
  formatted = "./output/04_formatted",
  models = "./output/05_models",
  evaluation = "./output/06_evaluation",
  ensemble = "./output/07_ensemble",
  projection = "./output/08_projection",
  paleoclimate = "./output/09_paleoclimate",
  visualization = "./output/10_visualization"
  
)

# 启动时自动创建输出目录
for (dir_name in output_dirs) {
  dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
}

# ==================== 5. 空间稀疏化 ====================

# 常用：是否在 `01_data_preparation` 中对 presence 做空间稀疏化
# 对存在明显采样聚集的数据，通常建议保持为 TRUE。
enable_spatial_thinning <- TRUE

# 常用：稀疏化方式
# 第一次使用时，`hybrid` 往往是最稳妥的选择。
# "cell"     = 每个气候栅格单元只保留 1 个点
# "distance" = 仅按最小距离阈值稀疏化
# "hybrid"   = 先按栅格去冗余，再按距离约束进一步筛选
thinning_method <- "hybrid"

# 高级：做 cell 稀疏化时，是否使用现代气候栅格作为网格参照
# 若希望 occurrence 与预测栅格严格对齐，可保持 TRUE。
thinning_use_predictor_cell <- TRUE

# 常用：距离阈值来源
# "auto"   = 根据现代气候栅格分辨率自动估算
# "manual" = 直接使用 thinning_distance_km
thinning_distance_mode <- "manual"

# 高级：自动阈值计算公式
# 现代气候栅格分辨率（km） * thinning_distance_factor
thinning_distance_factor <- 1.0

# 常用：手动指定最小保留距离（km）
# 这是控制 occurrence 稀疏程度的核心参数之一。
thinning_distance_km <- 20

# 高级：当多个点相互冲突时，按什么顺序决定“保留谁”
# "random" = 随机顺序贪心保留（推荐，可降低输入顺序偏差）
# "first"  = 按原始输入顺序保留
thinning_keep <- "random"

# 一般不用动：稀疏化随机种子（保证可重复）
thinning_seed <- 42

# 高级：诊断输出中统计的最近邻距离阈值（km）
# 仅影响诊断表和图，不改变实际稀疏化结果。
thinning_diagnostic_thresholds_km <- c(10, 20, 30)

# 高级：最近邻距离直方图的横轴上限（km）
thinning_histogram_xlim_km <- 100

# ==================== 6. 伪缺失点 ====================

# 常用：伪缺失点生成策略
# 第一次接触本流程时，通常直接使用 `kde` 即可。
# "kde"    = 基于 occurrence 密度的加权抽样
# "disk"   = 在 presence 周围给定距离范围内抽样
# "random" = 随机抽样
# "sre"    = 基于表面范围包络（Surface Range Envelope）
pa_strategy <- "kde"

# 高级：KDE 基本参数（仅在 pa_strategy = "kde" 时生效）
# 这里主要控制核密度估计的数学细节；第一次使用一般不用改。
kde_sigma <- "Scott"  # 带宽选择: "Scott", "Stoyan", "Diggle", "likelihood", "geometry"
kde_window <- "hull"  # 窗口类型: "hull" (凸包), "rect" (矩形)

# 常用：自定义 KDE 实现参数（兼容较低版本 R 环境）
# 带宽（km）：控制 occurrence 热点向外扩散的平滑尺度
# 带宽越大，伪缺失点分布越平滑、越分散；越小则越贴近 occurrence 聚集区。
kde_bandwidth_km <- 2500

# 高级：KDE 权重图的计算分辨率（km）
# 数值越大，计算越快但越粗；设为 NULL 时使用原始气候栅格分辨率
kde_compute_resolution_km <- 60

# 高级：是否在 KDE 粗网格内加入随机抖动，以减弱“块状感”
kde_jitter_within_cell <- T

# 高级：候选点过采样倍数
# 数值越大，越容易满足最小距离约束，但运行也会更慢
kde_candidate_multiplier <- 8

# 高级：最大重采样次数
kde_max_iter <- 10

# 高级：SRE 参数（仅在 pa_strategy = "sre" 时生效）
sre_quant <- 0.1      # 半分位数

# 常用：伪缺失点与 presence 的距离限制（单位：经纬度度数）
# 这两个参数直接决定伪缺失点离已知 occurrence 有多近、多远。
# 主要对 `disk` 策略生效；在当前实现里也会被 `kde` 方案用于距离约束。
pa_dist_min <- 5  # 1 度 ≈ 111 km
pa_dist_max <- 60 # KDE/disk 伪缺失点距 presence 的最大距离（度）

# 高级：伪缺失点抽样范围限制
# 用于剔除高纬区域，避免抽到明显不合理的全球极区样本
pa_lat_limit_north <- 75   # 剔除北纬 75°N 以北 (如北极点周边)
pa_lat_limit_south <- -60  # 剔除南纬 60°S 以南 (如南极洲)

# 常用：伪缺失点方案 1：1:1 比例
# presence 数量与 pseudo-absence 数量相同，适合先做基线测试。
# 注意：这里的数量不是写死的；`nb.absences = "equal"` 会根据
# 当前 presence 数量自动设置为相同数量的 pseudo-absence。
pa_scheme1 <- list(
  name = "scheme_1to1",
  nb.absences = "equal",
  nb.rep = 1,
  desc = "1:1 ratio (equal to presence)"
)

# 常用：伪缺失点方案 2：约 1:3 比例
# 当希望增加背景约束、提高模型区分能力时可尝试。
# 注意：`1000` 只是当前手动设定的示例值，不代表通用推荐值；
# 实际需要根据 presence 数量和研究设计手动调整。
pa_scheme2 <- list(
  name = "scheme_1000",
  nb.absences = 2000,
  nb.rep = 1,
  desc = "1:10 ratio (2000 points)"
)

# ==================== 7. 数据格式化 ====================

# 常用：选择用于后续建模的数据方案
# 主流程目前只保留两种方案，便于初次使用时理解和维护。
# "scheme_1to1" = 与 presence 等量
# "scheme_1000" = 固定 1000 个伪缺失点
pa_scheme_to_use <- "scheme_1000"

# ==================== 8. 变量筛选 ====================

# 常用：相关性筛选阈值（初筛）
# 值越低，筛选越严格；值越高，保留变量越多。
corr_threshold <- 0.8

# 常用：VIF 阈值（终筛）
# 用于进一步控制多重共线性，是变量筛选的核心参数之一。
vif_threshold <- 11

# 常用：是否执行 VIF 逐步剔除
use_vif_step <- TRUE

# 常用：手动预设的候选变量
# 设为 NULL 表示完全依赖“相关性 + VIF”自动筛选
# BioClim 常用变量释义：
# 如果你对物种生态学有明确判断，可以先在这里限定候选变量范围。
#   Bio1 = 年均温, Bio2 = 昼夜温差月均值, Bio3 = 等温性, Bio4 = 温度季节性
#   Bio5 = 最暖月最高温, Bio6 = 最冷月最低温, Bio7 = 温度年较差
#   Bio8 = 最湿季度均温, Bio9 = 最干季度均温, Bio10 = 最暖季度均温, Bio11 = 最冷季度均温
#   Bio12 = 年降水, Bio13 = 最湿月降水, Bio14 = 最干月降水, Bio15 = 降水季节性
#   Bio16 = 最湿季度降水, Bio17 = 最干季度降水, Bio18 = 最暖季度降水, Bio19 = 最冷季度降水
manual_var_select <- c(#'Bio1',
                       'Bio2', 'Bio5','Bio6',#'Bio7','Bio3','Bio4',
                       'Bio8','Bio9',# 'Bio10', 'Bio11',
                       #'Bio12',#
                       'Bio13','Bio14','Bio15'#'Bio16','Bio17','Bio18','Bio19'
                       )
# ==================== 9. 数据划分与交叉验证 ====================

# 常用：是否启用交叉验证
# 对正式建模，通常建议保持 TRUE。
use_cv <- TRUE

# 常用：交叉验证策略
# 本流程主要面向空间交叉验证，因此通常推荐 `user.defined`。
# "kfold"        = 常规 k 折交叉验证
# "strat"        = 分层抽样
# "random"       = 随机划分
# "user.defined" = 使用 `04_data_formatting.R` 生成的自定义空间折
cv_strategy <- "user.defined"

cv_nb_rep <- 3        # 常用：重复次数；在 user.defined 下表示生成几套平移网格
cv_k <- 5             # 常用：每套网格划分成多少折
cv_perc <- 0.8        # 一般不用动：仅非 user.defined 策略使用

# 常用：空间 CV 具体方法
# "checkerboard" = 基于空间网格的交错分配
spatial_cv_method <- "checkerboard"

# 常用：空间块大小（km）
# 通常建议明显大于稀疏化距离，否则训练 / 验证空间分离不充分
# 这是空间交叉验证中最值得优先理解和调整的参数之一。
spatial_cv_block_size_km <- 500

# 高级：空间分块时使用的投影 CRS（单位需为米）
spatial_cv_projection_crs <- "EPSG:3857"

# 高级：若某些折缺少 presence 或 pseudo-absence，自动放大 block size 重试
spatial_cv_auto_retry <- TRUE
spatial_cv_retry_factor <- 1.5
spatial_cv_max_tries <- 4

# ==================== 10. 模型算法 ====================

# 常用：单模型算法列表
# 第一次接触流程时，不必追求算法越多越好；先保证流程稳定更重要。
# 可选："ANN", "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF", "RFd"
model_list <- c( "CTA", "FDA", "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF")
# model_list <- c("MAXNET")
# 如果部分算法在当前环境不可用，可临时删减
# model_list <- c("GAM", "GBM", "GLM", "MAXNET", "RF")

# 常用：快速测试配置（只跑 2 个模型）
# 用于排查流程是否能完整跑通。
model_list_test <- c("MAXNET", "RF")

# ==================== 11. 模型训练 ====================

# 常用：参数优化策略
# "bigboss" = biomod2 推荐参数
# "user"    = 使用用户自定义参数
# "none"    = 使用包默认参数
opt_strategy <- "bigboss"

# 常用：模型评估指标
# 这些指标既影响单模型评估，也会影响后续筛选与集成。
metric_eval <- c("TSS", "AUCroc", "AUCprg", "BOYCE")

# 高级：变量重要性置换次数
# 数值越大越稳，但耗时也越高。
var_importance_rep <- 5

# 常用：并行使用的 CPU 核数
nb_cpu <- 8

# 一般不用动：全局随机种子（保证可重复）
seed_val <- 123

# ==================== 12. 集成模型 ====================

# 常用：集成算法
# 第一次跑流程时，通常只保留 EMwmean 即可。
# "EMmedian" = 中位数集成
# "EMmean"   = 平均值集成
# "EMwmean"  = 加权平均集成（当前推荐）
# "EMca"     = 委员会平均
# "EMci"     = 置信区间
# "EMcv"     = 交叉验证集成
# 若目标是得到最简洁、稳定的结果，通常只保留 EMwmean 即可
ensemble_algo <- c("EMwmean")

# 常用：模型筛选指标
# 只有同时满足这些指标阈值的模型，才会进入后续集成。
# 可选："TSS", "AUCroc", "AUCprg", "BOYCE"
metric_select <- c("TSS", "AUCroc","AUCprg", "BOYCE")

# 常用：模型筛选阈值
# 顺序必须与 metric_select 一一对应
# 当前规则：
#   TSS >= 0.8
#   AUCroc >= 0.8
#   AUCprg >= 0.8
#   BOYCE >= 0.6
metric_select_thresh <- c(0.8, 0.8, 0.8, 0.6)

# 高级：集成方式
# "all"     = 全部模型合并为一个集成模型（推荐）
# "algo"    = 按算法分别集成
# "all_run" = 按重复 run 分别集成
em.by <- "all"

# ==================== 13. 单模型筛选 ====================

# 常用：模型筛选模式
# 初学者如果还不确定阈值是否合理，先用 manual 能更直观看结果。
# "auto"   = 按阈值自动筛选
# "manual" = 手动指定保留模型
model_selection_mode <- "manual"

# 常用：手动模式下保留的模型名称（可写一个或多个）
# 例如: manual_selected_models <- c("MAXNET", "RF", "GLM")
manual_selected_models <- c("CTA", "GAM", "GLM", "MARS", "MAXNET", "RF")

# 高级：自动模式下最多保留前 N 个模型
top_n_models <- 5

# 高级：EMwmean 权重衰减方式
# "proportional" = 按指标成比例加权
# "exponential"  = 指数放大高分模型权重
emwmean_decay <- "proportional"

# ==================== 14. 投影设置 ====================

# 常用：投影版本标签
# 用于区分不同运行结果；会直接影响输出目录和文件命名
# 例如："v1_modern", "v2_teyf", "v3_tewt_petm"
projection_name <- "v10_test_sudo_number"

# 常用：需要投影的古气候时期
# 设为 NULL 时自动从数据目录扫描并投影全部时期。
# 也可以手动指定：
#   1. 字符向量，如 c("tewta_0k", "tewtb_1k")
#   2. 时间范围，如 list(time_range = c(0, 100))
#      表示只投影 0-100 ka（含边界）的时期
# 当前 09_paleoclimate_projection.R 会从时期名称中提取数字来做范围筛选。
# 对 TEWT 数据，试跑 100 ka 以内可写：
#   projection_periods <- list(time_range = c(0, 100))
projection_periods <- list(time_range = c(3602, 999999)) # 动态读取，设为 NULL 自动扫描

# 一般不用动：投影栅格输出格式
# ".tif" = GeoTIFF（推荐）
# ".grd" = Raster 原生格式
output_format <- ".tif"

# 常用：连续适宜度转二值分布图时使用的阈值规则
# "TSS" = 基于最大 TSS
# "ROC" = 基于最大 ROC
# "all" = 输出所有可用规则
metric_binary <- "TSS"

# ==================== 15. 可视化 ====================

# 常用：地图显示范围
# 只影响绘图显示，不影响建模本身。
map_extent <- c(-180, 180, -45, 35)  # 经度 min, max; 纬度 min, max

# 高级：适宜度颜色带
suitability_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(100)

# 高级：出现点样式
presence_color <- "#336600"
presence_pch <- 19

# 高级：伪缺失点样式
pseudoabsence_color <- "#FF6B6B"

# 常用：可视化时优先使用哪一种集成 metric 对应的投影图
# 可选："TSS", "AUCroc", "AUCprg", "BOYCE"
# 10_visualization.R 会优先为每个时期选择 `EMwmeanBy<metric>` 对应的 tif；
# 如果该 metric 不存在，再自动回退到其他可用结果。
visualization_projection_metric <- "TSS"

# 高级：标题字号
title_size <- 14

# ==================== 16. 地理参数 ====================

# 一般不用动：默认地理坐标系
target_crs <- "EPSG:4326"  # WGS84

# 高级：默认缓冲区大小（度）
buffer_degree <- 0.5

# ==================== 17. 其他设置 ====================

# 常用：是否保存中间结果
# 第一次跑流程时建议保持 TRUE，方便排错和回溯。
save_intermediate <- TRUE

# 常用：是否输出详细日志
verbose <- TRUE

# 高级：内存管理
# 每处理若干时期后触发一次 gc()，主要影响长时间批量投影时的稳定性。
gc_trigger <- 20  # 每处理 N 个时期后调用 gc()

# ============================================================================
# 依赖包
# ============================================================================

required_packages <- c(
  # 核心建模 / 空间数据
  "biomod2",
  "terra",
  "raster",
  "sf",
  
  # 数据整理
  "dplyr",
  "stringr",
  "tidyr",
  "readr",
  "gtools",        # mixedsort 函数
  
  # 可视化
  "ggplot2",
  "plotly",
  "tidyterra",
  "viridis",      # 颜色映射
  "corrplot",
  "magick",
  
  # 统计与辅助建模
  "usdm",          # VIF 计算
  # "spatialEco",    # pseudo.absence 版本不对
  "gamlss",
  
  # 并行
  "doParallel",
  "foreach",
  
  # 其他工具
  "reshape",
  "ggtext",
  "htmlwidgets",
  "R.utils",
  "lme4",
  "cito", #版本不对，手动安装
  "xgboost" #版本不对，手动安装
)

# 一般不用动：检查并加载依赖包
# 注意：这里会自动安装缺失包；若服务器环境受限，建议提前手动安装。
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# # 高级：低版本 R 环境下的手动安装示例
# Matrixurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-1.tar.gz"
# 
# # 从源码安装指定版本
# install.packages(Matrixurl, repos = NULL, type = "source")

# # 高级：xgboost 1.7.7.1 历史版本下载地址
# url_xgb <- "https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.7.1.tar.gz"
# 
# # 从源码安装指定版本
# install.packages(url_xgb, repos = NULL, type = "source")

# ============================================================================
# 全局辅助函数
# ============================================================================

#' 标准化加载气候数据
#' 自动读取目录中的栅格文件，按自然数字顺序排序后重命名为 Bio1, Bio2...
#' @param clim_dir 气候数据所在目录
#' @param pattern 文件名匹配模式，默认为 "*.tif$"
#' @return terra::SpatRaster 对象
load_standard_climate <- function(clim_dir, pattern = "*.tif$") {
  files <- list.files(clim_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop(paste("Error: No climate files found in", clim_dir))
  }
  
  # 读取栅格
  clim <- terra::rast(files)
  
  # 按自然顺序重新排列图层，避免出现 1, 10, 11, 2 这种字典序问题
  clim <- clim[[gtools::mixedsort(names(clim))]]
  
  # 统一图层命名
  names(clim) <- paste0("Bio", 1:terra::nlyr(clim))
  
  return(clim)
}

# ============================================================================
# 空间稀疏化辅助函数
# ============================================================================

#' 估算栅格分辨率（km）
#' 根据第一层气候栅格的经纬度分辨率，换算研究区中心附近的公里尺度
#' @param clim SpatRaster 对象
#' @param ref_points 可选，SpatVector 点数据，用于估算中心纬度
#' @return list(x_km, y_km, mean_km)
estimate_raster_resolution_km <- function(clim, ref_points = NULL) {
  raster_res <- terra::res(clim)
  
  if (!is.null(ref_points) && inherits(ref_points, "SpatVector")) {
    center_lat <- mean(terra::crds(ref_points)[, 2], na.rm = TRUE)
  } else {
    center_lat <- mean(terra::ext(clim)[3:4])
  }
  
  lat_km <- 111.32 * abs(raster_res[2])
  lon_km <- 111.32 * cos(center_lat * pi / 180) * abs(raster_res[1])
  
  list(
    x_km = lon_km,
    y_km = lat_km,
    mean_km = mean(c(lon_km, lat_km), na.rm = TRUE)
  )
}

#' 根据点集中心自动猜测 UTM EPSG
#' @param points SpatVector 点数据（经纬度）
#' @return EPSG 字符串
guess_utm_epsg <- function(points) {
  coords <- terra::crds(points)
  center_lon <- mean(coords[, 1], na.rm = TRUE)
  center_lat <- mean(coords[, 2], na.rm = TRUE)
  zone <- floor((center_lon + 180) / 6) + 1
  
  if (center_lat >= 0) {
    return(sprintf("EPSG:%d", 32600 + zone))
  }
  
  sprintf("EPSG:%d", 32700 + zone)
}

#' 按气候栅格单元进行稀疏化
#' 每个 raster cell 最多保留 1 个点
#' @param points SpatVector 点数据
#' @param clim_ref 参考气候栅格
#' @param keep 策略："random" 或 "first"
#' @param seed 随机种子
#' @return list(points, kept_index, cell_ids)
thin_points_by_cell <- function(points, clim_ref, keep = "random", seed = 123) {
  stopifnot(inherits(points, "SpatVector"))
  stopifnot(inherits(clim_ref, "SpatRaster"))
  
  coords <- terra::crds(points)
  cell_ids <- terra::cellFromXY(clim_ref[[1]], coords)
  valid_idx <- which(!is.na(cell_ids))
  
  if (length(valid_idx) == 0) {
    stop("Error: No occurrence points overlap with the predictor raster extent.")
  }
  
  idx_pool <- valid_idx
  if (keep == "random") {
    set.seed(seed)
    idx_pool <- sample(idx_pool)
  }
  
  seen_cells <- integer(0)
  keep_idx <- integer(0)
  
  for (idx in idx_pool) {
    cell_id <- cell_ids[idx]
    if (!(cell_id %in% seen_cells)) {
      keep_idx <- c(keep_idx, idx)
      seen_cells <- c(seen_cells, cell_id)
    }
  }
  
  keep_idx <- sort(keep_idx)
  
  list(
    points = points[keep_idx, ],
    kept_index = keep_idx,
    cell_ids = cell_ids[keep_idx]
  )
}

#' 按最小距离进行贪心稀疏化
#' @param points SpatVector 点数据
#' @param distance_km 最小距离阈值（km）
#' @param keep 策略："random" 或 "first"
#' @param seed 随机种子
#' @return list(points, kept_index, crs_used)
thin_points_by_distance <- function(points, distance_km, keep = "random", seed = 123) {
  stopifnot(inherits(points, "SpatVector"))
  
  if (terra::nrow(points) <= 1) {
    return(list(points = points, kept_index = seq_len(terra::nrow(points)), crs_used = terra::crs(points)))
  }
  
  utm_epsg <- guess_utm_epsg(points)
  points_proj <- terra::project(points, utm_epsg)
  coords <- terra::crds(points_proj)
  threshold_m <- distance_km * 1000
  
  order_idx <- seq_len(nrow(coords))
  if (keep == "random") {
    set.seed(seed)
    order_idx <- sample(order_idx)
  }
  
  keep_idx <- integer(0)
  kept_coords <- matrix(numeric(0), ncol = 2)
  
  for (idx in order_idx) {
    candidate <- coords[idx, , drop = FALSE]
    
    if (length(keep_idx) == 0) {
      keep_idx <- c(keep_idx, idx)
      kept_coords <- rbind(kept_coords, candidate)
      next
    }
    
    distances <- sqrt((kept_coords[, 1] - candidate[1, 1])^2 + (kept_coords[, 2] - candidate[1, 2])^2)
    
    if (all(distances >= threshold_m)) {
      keep_idx <- c(keep_idx, idx)
      kept_coords <- rbind(kept_coords, candidate)
    }
  }
  
  keep_idx <- sort(keep_idx)
  
  list(
    points = points[keep_idx, ],
    kept_index = keep_idx,
    crs_used = utm_epsg
  )
}

#' 计算每个点的最近邻距离（km）
#' @param points SpatVector 点数据
#' @return numeric vector
nearest_neighbor_distances_km <- function(points) {
  if (!inherits(points, "SpatVector") || terra::nrow(points) < 2) {
    return(NA_real_)
  }
  
  utm_epsg <- guess_utm_epsg(points)
  points_proj <- terra::project(points, utm_epsg)
  coords <- terra::crds(points_proj)
  dist_mat <- as.matrix(stats::dist(coords)) / 1000
  diag(dist_mat) <- Inf
  apply(dist_mat, 1, min, na.rm = TRUE)
}

#' 汇总最近邻距离统计量（km）
#' @param points SpatVector 点数据
#' @return named numeric vector
nearest_neighbor_summary_km <- function(points) {
  nn <- nearest_neighbor_distances_km(points)
  
  if (length(nn) == 1 && is.na(nn)) {
    return(c(min_km = NA_real_, median_km = NA_real_, mean_km = NA_real_))
  }
  
  c(
    min_km = min(nn, na.rm = TRUE),
    median_km = stats::median(nn, na.rm = TRUE),
    mean_km = mean(nn, na.rm = TRUE)
  )
}

#' 确保点数据具有 CRS
#' 若输入点缺失 CRS，则回填配置中的 target_crs
#' @param points SpatVector 点数据
#' @param default_crs 默认 CRS
#' @return 带 CRS 的 SpatVector
ensure_points_crs <- function(points, default_crs = target_crs) {
  stopifnot(inherits(points, "SpatVector"))
  
  current_crs <- terra::crs(points)
  if (is.null(current_crs) || !nzchar(current_crs)) {
    terra::crs(points) <- default_crs
  }
  
  points
}

# ============================================================================
# 古气候数据加载函数
# ============================================================================

#' 列出所有可用古气候时期
#' @param clim_dir 古气候数据根目录 (默认使用 past_clim_dir)
#' @param subdir 是否扫描子目录 (默认使用 past_clim_subdir)
#' @return 时期名称向量
list_paleo_periods <- function(clim_dir = past_clim_dir, subdir = past_clim_subdir) {
  if (subdir) {
    # 扫描时期子目录
    period_dirs <- list.dirs(clim_dir, full.names = FALSE, recursive = FALSE)
    # 如果目录名带有 "_masked" 后缀，则去掉后再返回
    periods <- gsub("_masked$", "", period_dirs)
    return(sort(unique(periods))) # 去重并排序
  } else {
    # 直接在目录内扫描文件名
    clim_files <- list.files(clim_dir, pattern = "\\.(tif|nc)$", full.names = FALSE)
    periods <- unique(gsub("_bio\\d+.*$", "", clim_files))
    return(sort(periods))
  }
}

#' 加载指定时期的古气候数据
#' @param period 时期名称 (如 "tewta_0k", "tewtb_1000k")
#' @param clim_dir 古气候数据根目录 (默认使用 past_clim_dir)
#' @param subdir 时期数据是否在子目录中 (默认使用 past_clim_subdir)
#' @param pattern 文件匹配模式 (默认 "*.nc$|*.tif$")
#' @return SpatRaster 对象
#' 会自动兼容 `period/` 和 `period_masked/` 两种目录结构
load_paleo_climate <- function(period, 
                               clim_dir = past_clim_dir, 
                               subdir = past_clim_subdir,
                               pattern = "\\.(nc|tif)$") {
  
  cat("  Loading paleoclimate for period:", period, "\n")
  
  if (subdir) {
    # 自动尝试两种常见目录命名方式
    path_plain <- file.path(clim_dir, period)
    path_masked <- file.path(clim_dir, paste0(period, "_masked"))
    
    if (dir.exists(path_plain)) {
      period_dir <- path_plain
    } else if (dir.exists(path_masked)) {
      period_dir <- path_masked
    } else {
      stop("Error: Neither '", period, "' nor '", period, "_masked' directory exists in ", clim_dir)
    }
    
    clim_files <- list.files(period_dir, pattern = pattern, full.names = TRUE)
    
    if (length(clim_files) == 0) {
      stop("Error: No climate files found in ", period_dir)
    }
    
    cat("    Found", length(clim_files), "files in", basename(period_dir), "\n")
    clim <- terra::rast(clim_files)
    
  } else {
    # 文件直读模式
    clim_files <- list.files(clim_dir, 
                             pattern = paste0("^", period, "_.*", pattern), 
                             full.names = TRUE)
    if (length(clim_files) == 0) stop("Error: No files found for ", period)
    clim <- terra::rast(clim_files)
  }
  
  # 变量名标准化：兼容 bio1 / bio01 / bio_1 / bio_01 等格式
  layer_names <- names(clim)
  standard_names <- gsub(".*bio_?0*(\\d+).*", "Bio\\1", layer_names, ignore.case = TRUE)
  names(clim) <- standard_names
  
  cat("    Standardized layers:", paste(names(clim), collapse = ", "), "\n")
  
  return(clim)
}

#' 获取投影输出目录
#' 根据 projection_name 自动拼接并创建目录
#' @param type 投影类型 ("modern" 或 "paleo")
#' @param period 古气候时期名称 (仅 paleo 需要)
#' @return 输出目录路径
get_projection_output_dir <- function(type = c("modern", "paleo"), period = NULL) {
  type <- match.arg(type)
  
  if (type == "modern") {
    out_dir <- file.path(output_dirs$projection, paste0("modern_", projection_name))
  } else {
    if (is.null(period)) {
      stop("Error: period is required for paleo projection")
    }
    out_dir <- file.path(output_dirs$projection, paste0("paleo_", period, "_", projection_name))
  }
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  return(out_dir)
}

# ============================================================================
# 启动时打印当前配置摘要
# ============================================================================

cat("\n")
cat("=================================================================\n")
cat("  Sargentodoxa SDM - Configuration Loaded\n")
cat("=================================================================\n")
cat("  Running mode: Local (de-servered repo; climate rasters via env vars)\n")
cat("  Project:", project_dir, "\n")
cat("  Species:", species_name, "\n")
cat("  Data file:", sp_file, "\n")
cat("  Modern climate:", modern_clim_dir, "\n")
cat("  Past climate:", past_clim_dir, "\n")
cat("  Past climate source:", current_past_clim, "\n")
cat("  Past climate subdirs:", past_clim_subdir, "\n")
cat("\n")
cat("  Pseudo-absence strategy:", pa_strategy, "\n")
cat("  Pseudo-absence scheme:", pa_scheme_to_use, "\n")
cat("  Variable selection: VIF threshold =", vif_threshold, "\n")
cat("  Cross-validation:", cv_strategy, " (k =", cv_k, ", rep =", cv_nb_rep, ")\n")
cat("\n")
cat("  Models:", paste(model_list, collapse = ", "), "\n")
cat("  Model selection mode:", model_selection_mode, "\n")
if (model_selection_mode == "manual") {
  cat("  Manual selected models:", paste(manual_selected_models, collapse = ", "), "\n")
} else {
  cat("  Auto selection: Top", top_n_models, "models (TSS >=", metric_select_thresh[1], ", AUC >=", metric_select_thresh[2], ")\n")
}
cat("  Ensemble:", paste(ensemble_algo, collapse = ", "), "\n")
cat("\n")
cat("  Output base:", output_base, "\n")
cat("  Projection name:", projection_name, "\n")
cat("=================================================================\n")
cat("\n")
