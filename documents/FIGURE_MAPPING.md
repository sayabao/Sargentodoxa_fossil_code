# Article figure ↔ code mapping

整理日期：2026-07-28

本文档映射 v5.3 手稿每张「代码生成」的图到本仓库的脚本与输出。
非代码图（分布图、形态 / 标本照片）不在此列，见末尾。

## 正文图

| 文章图 | 内容 | 代码组件 | 脚本 | 输出 |
|---|---|---|---|---|
| Fig 3 | 变量重要性 + 重要变量古气候背景分异 | M4 + R7 | `visualization/01_…`（M4）、`04_server_…`+`06_`+`07_`（R7）、`08_combine_…` | `output/figures_paper/plates/Figure_1_M4_R7_variable_importance_and_background.*` |
| Fig 4 | 生态位模拟：古气候适宜图 + 面积趋势 + 纬度梯度 | R2b + R3 + R3b | `visualization/01_…`（R3, R3b）、`02_…`（R2b）、`08_combine_…` | `…/plates/Figure_2_R3_R3b_temporal_latitudinal_dynamics.*`、`…/plates/Figure_3_R2b_all_period_suitability_panels.*` |
| Fig 5 | 化石点与投影对比 + 逐点适宜度 | R4 + R5 | `visualization/03_…`、`08_combine_…` | `…/plates/Figure_4_R4_R5_fossil_model_consistency.*` |

## 附录图

| 文章图 | 内容 | 代码组件 | 脚本 | 输出 |
|---|---|---|---|---|
| Fig S1 | 出现点 / 伪缺失点采样图 | M1 | `visualization/00_…` | `…/model_diagnostics/M1_occurrence_pseudoabsence_map.*` |
| Fig S2 | 变量相关矩阵 + VIF | M2 | `visualization/01_…` | `…/model_diagnostics/M2_variable_correlation_vif.*` |
| Fig S3 | 模型评价指标 | M3 | `visualization/01_…` | `…/model_diagnostics/M3_model_evaluation_metrics.*` |
| Fig S4 | Presence vs Background | 08_projection 内置图 | `08_projection.R`（`plot_presence_background_comparison`） | `output/08_projection/projection_modern_map.pdf`（多页 PDF 中的 Presence vs Background 页） |
| Fig S10 | 现代适宜区图 | R1 | `visualization/02_…` | `…/maps/R1_modern_suitability_map.*` |
| Fig S11 | 化石点适宜度显著性检验 | R8 | `visualization/09_…` | `…/fossil_comparison/R8_fossil_significance_figures10.*` |

## 非代码图（不在本仓库）

| 文章图 | 内容 | 来源 |
|---|---|---|
| Fig 1 | 化石记录与现代分布 | 手工制图（PSD） |
| Fig 2 | 新化石与修订标本照片 | 标本照片 |
| Fig S5–S9 | 各国化石标本照片（guipingensis/cuneata、France、Poland、Germany、Japan） | 标本照片 |

## 阈值与口径（图版一致性）

- 适宜度分级：0–0.1 / 0.1–0.3 / 0.3–0.7 / 0.7–1（least / slightly / moderately / highly suitable）
- 高适宜区阈值：0.7（R3b、R7 均用此定义）
- biomod2 投影原始为 0–1000 存储尺度，出图前转换为 0–1
- 化石记录 16 条，与手稿 Table 1 一致；Tajimi 记录在分析文件中标为 *cuneata*，手稿已修订为 *tajimiensis* sp. nov.（仅标签，不影响模型）

## 与原工作目录图版计划的差异

原工作目录的 `PAPER_FIGURE_PLAN.md`（2026-05-21）规划的是 4 个合并图版
（Fig1=M4+R7, Fig2=R3+R3b, Fig3=R2b, Fig4=R4+R5, S1=M1+M2+M3）。v5.3 手稿最终采用
5 主图 + 11 附图、多数单独成图。本仓库以 **v5.3 实际图集**为准（见上表）。
