# Reproducibility

整理日期：2026-07-28

## 总览

| 阶段 | 脚本 | 本地可复现？ | 依赖 |
|---|---|---|---|
| 01 数据准备 / 空间稀疏化 | `01_data_preparation.R` | 否（需完整 WorldClim） | WorldClim 2.1 5min 全部 Bio 栅格（服务器） |
| 02 伪缺失点 | `02_pseudo_absence.R` | 否（依赖 01） | 01 输出 + WorldClim |
| 03 变量筛选 | `03_variable_selection.R` | 否（依赖 01） | 01 输出 + WorldClim |
| 04 数据格式化 / 空间 CV | `04_data_formatting.R` | 否（依赖 01–03） | 前置输出 |
| 05 模型训练 | `05_model_training.R` | 否（依赖 04） | 04 格式化数据 |
| 06 模型评估 | `06_model_evaluation.R` | 否（依赖 05） | 05 模型 |
| 07 集成建模 | `07_ensemble_modeling.R` | 否（依赖 05–06） | 05/06 输出 |
| 08 现代投影 | `08_projection.R` | 否（需 WorldClim） | 07 集成模型 + WorldClim |
| 09 古气候投影 | `09_paleoclimate_projection.R` | 否（需 TEYF） | 07 集成模型 + TEYF 栅格（服务器） |
| 出图 M1–M4, R1–R8 | `scripts/visualization/*.R` | **是** | 本仓库 frozen `output/` + `data/` |

## 本地复现出图

所有 `scripts/visualization/` 脚本从仓库根目录运行，读取随仓库提供的
`output/`（frozen 模型/投影结果）与 `data/`（出现点、化石记录），不依赖任何
服务器数据。详见 `README.md` 的 Quickstart。

> 注：`output/04_formatted/biomod2_formatted_data.RData`（72 MB）与
> `output/05_models/model_results.RData`（73 MB）未随仓库提供——出图非必需，
> 已有提取好的 CSV（变量重要性、评估指标等）。如需重跑 05–07，需先在服务器完成 01–04。

## 服务器端重跑管线

如需从原始气候数据重跑 01–09，需要服务器上的气候栅格。在服务器上克隆本仓库后，
通过环境变量指定气候目录（默认值即服务器原始路径，故在原服务器上可直接运行）：

```bash
export SARGENTO_MODERN_CLIM_DIR="/home/cedrus/data/shz/data/wc2.1_5m_bio"
export SARGENTO_TEYF_CLIM_DIR="/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min"
# 可选：TEWT
export SARGENTO_TEWT_CLIM_DIR="/home/cedrus/data/shz/data/tewt_0_3600ka/masked_landsea"
```

然后从仓库根目录依次运行 `Rscript scripts/01_…_.R` … `09_…_.R`。
路径配置见 `scripts/config.R` 第 1 节。

## R7（重要变量古气候背景）特殊说明

`scripts/visualization/04_server_important_variable_paleoclimate_background.R`
需要原始 TEYF Bio 变量 nc 文件（服务器
`/home/cedrus/data/shz/data/from linlin/teyf-mask-nc-5min/`），设计为服务器端运行。
其输出（R7 统计表与图）已随仓库提供于 `output/figures_paper/`，本地重绘脚本
`06_redraw_r7_…`、`07_redraw_r7_contrast_…` 可直接基于这些输出重绘热图，无需服务器。

## 随机种子

- 空间稀疏化 / CV：`thinning_seed = 42`、`seed_val = 123`（见 `config.R`）
- 模型训练、集成、投影均受上述种子约束，保证可重复。

## 未随仓库提供的原工作目录内容

以下原工作目录内容未纳入本仓库（旧版 / 未用 / 过大）：

- `output/09_paleoclimate/`（旧版古气候投影，被 `09_paleoclimate_teyf/` 取代）
- `output/10_visualization/`、`output/10_visualization_teyf/`（旧可视化）
- `output/05_models_random/`（随机对照运行）
- `output/04_formatted/biomod2_formatted_data.RData`、`output/05_models/model_results.RData`（大文件，出图非必需）
- R6 化石地理位移图脚本与输出（与正文图功能重复，未纳入文章图版）
- 根目录历史脚本、`.RData` / `.Rhistory`
