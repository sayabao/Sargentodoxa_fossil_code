# Article figure ↔ code mapping

整理日期：2026-07-28（2026-07-28 按图清单精简后更新）

本文档映射 v5.3 手稿每张「代码生成」的图到本仓库的脚本与输出。
**本仓库只产出组件图；最终文章图由作者手动 PS 拼合**（组合图版 plates/ 已移除）。
非代码图（分布图、形态 / 标本照片）见末尾。

## 正文图（组件图 → 文章图，手动 PS 拼合）

| 文章图 | 内容 | 代码组件 | 脚本 | 组件输出（output/figures_paper/…） |
|---|---|---|---|---|
| Fig 3 | 变量重要性 + 重要变量古气候背景分异 | M4 + R7 | `visualization/01_…`（M4）、`04_server_…`+`07_…`（R7） | `model_diagnostics/M4_variable_importance.*`、`trends/R7_important_variable_background_contrast_heatmap.*` |
| Fig 4 | 生态位模拟：古气候适宜图 + 面积趋势 + 纬度梯度 | R2b + R3 + R3b | `visualization/01_…`（R3, R3b）、`02_…`（R2b） | `maps/R2b_all_period_suitability_panels.*`、`trends/R3_suitability_area_trend_from_earliest_fossil.*`、`trends/R3b_high_suitability_latitudinal_area_heatmap.*` |
| Fig 5 | 化石点与投影对比 + 逐点适宜度 | R4 + R5 | `visualization/03_…` | `fossil_comparison/R4_fossil_projection_panels.*`、`fossil_comparison/R5_fossil_point_suitability.*` |

## 附录图

| 文章图 | 内容 | 代码组件 | 脚本 | 组件输出 |
|---|---|---|---|---|
| Fig S1 | 出现点 / 伪缺失点采样图 | M1 | `visualization/00_…` | `model_diagnostics/M1_occurrence_pseudoabsence_map.*` |
| Fig S2 | 变量相关矩阵 + VIF | M2 | `visualization/01_…` | `model_diagnostics/M2_variable_correlation_vif.*` |
| Fig S3 | 模型评价指标 | M3 | `visualization/01_…` | `model_diagnostics/M3_model_evaluation_metrics.*` |
| Fig S4 | Present-day suitability (A) + Presence vs Background (B) | 08_projection 内置图 | `08_projection.R` | `output/08_projection/projection_modern_map.pdf`：**A = 第 1 页**，**B = 第 6 页**（"Presence vs Background: Ensemble (EMwmean)" 左图）；作者手动加绿色虚线 |
| Fig S11 | 化石点适宜度显著性检验 | R8 | `visualization/09_…` | `fossil_comparison/R8_fossil_significance_figures10.*` |

## 非代码图（不在本仓库）

| 文章图 | 内容 | 来源 |
|---|---|---|
| Fig 1 | 化石记录与现代分布 | 手工制图（PSD） |
| Fig 2 | 新化石与修订标本照片 | 标本照片 |
| Fig S5–S9 | 各国化石标本照片（guipingensis/cuneata、France、Poland、Germany、Japan） | 标本照片 |
| **Fig S10** | **现生 *S. cuneata* 形态特征**（聚合果、浆果、种子、SEM） | 标本照片（**非**现代适宜区图；此前误标为 R1，已更正） |

> 注：R1（现代适宜区图）未进入文章图版--现代适宜度已在 Fig 4 的 R2b 0 Ma 面板展示，故 R1 已从仓库移除。

## 阈值与口径（图版一致性）

- 适宜度分级：0–0.1 / 0.1–0.3 / 0.3–0.7 / 0.7–1（least / slightly / moderately / highly suitable）
- 高适宜区阈值：0.7（R3b、R7 均用此定义）
- biomod2 投影原始为 0–1000 存储尺度，出图前转换为 0–1
- 化石记录 16 条，与手稿 Table 1 一致；v2.1 文件已采用手稿终版分类学标签（*tajimiensis*、*guipingensis*、*cretaceum* comb. nov. 正式作者），仅标签更新不影响模型

## 已移除的图版（2026-07-28 精简）

作者实际使用手动 PS 拼图，故以下已从仓库移除：
- 全部组合图版 `plates/`（Figure_1…Figure_4、Appendix_Figure_S1）及生成脚本 `08_combine_paper_figure_plates.R`
- R1 现代适宜区图、R2 关键时期图（被 R2b 取代）
- R3b 探索版（latitudinal_gradient、curves、continuous_shift、suitable_03_1_shift、hemisphere_shift 等）及 R7 探索版（paleoclimate_background_heatmap、zscore、contrast_trends）
- R6 化石地理位移图（与正文图重复）及脚本 `05_…`、`06_…`

## 与原工作目录图版计划的差异

原工作目录的 `PAPER_FIGURE_PLAN.md`（2026-05-21）规划的是 4 个合并图版
（Fig1=M4+R7, Fig2=R3+R3b, Fig3=R2b, Fig4=R4+R5, S1=M1+M2+M3）。v5.3 手稿最终采用
5 主图 + 11 附图、多数单独成图，且作者手动 PS 拼合。本仓库以 **v5.3 实际图集**为准，
只保留各文章图对应的组件图。
