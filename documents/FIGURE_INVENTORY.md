# Figure inventory - 请标注保留哪些

> **2026-07-28 已处理**：已按你的标注精简脚本与输出。保留的组件图：M1/M2/M3/M4、R2b、R3、R3b_high_suitability_latitudinal_area_heatmap、R7_background_contrast_heatmap、R4/R5/R8；Fig S4 取自 `output/08_projection/projection_modern_map.pdf` 第 1 页(A)+第 6 页(B)。删除：R1、R2、R3 另一版、R3b 其余变体、R7 其余变体、R6、全部 plates/，及脚本 05/06/08。Fig S10 经核对为现生形态照片（非代码图），R1 确认未使用。最终映射见 `FIGURE_MAPPING.md`。

本表列出本仓库 `output/figures_paper/` 当前生成的所有图（每个组件图同时有 `.png` 和 `.pdf`）。
「v5.3 建议」列是我根据手稿 v5.3 图集的推断（见 `FIGURE_MAPPING.md`）。
**请在「你的决定」列填 `保留` / `删除`，或直接编辑本表**；确认后我会精简脚本（删除不要的出图代码与输出）。

> 注：你说「仓库生成的图和实际用的图不一样」，所以以你的标注为准。组合图版（plates/）与组件图（model_diagnostics/maps/trends/fossil_comparison/）可能重复，请说明你实际用的是哪一组。

## 组件图

| 文件（output/figures_paper/…） | 内容 | v5.3 建议 | 你的决定 |
|---|---|---|---|
| model_diagnostics/M1_occurrence_pseudoabsence_map | 出现点/伪缺失点采样图 | → Fig S1 |保留 |
| model_diagnostics/M2_variable_correlation_vif | 变量相关矩阵 + VIF | → Fig S2 |保留 |
| model_diagnostics/M3_model_evaluation_metrics | 模型评价指标 | → Fig S3 |保留 |
| model_diagnostics/M4_variable_importance | 变量重要性 | → Fig 3 组件(A) |保留 |
| maps/R1_modern_suitability_map | 现代适宜区图 | → Fig S10 |此图未使用，Fig S10理论上是modern 的结构图，需要你重新核对 |
| maps/R2_paleoclimate_key_period_panels | 关键时期适宜图（7 时期） | 未用（被 R2b 取代） |未使用 |
| maps/R2b_all_period_suitability_panels | 全时期适宜图（69–0 Ma，15 面板） | → Fig 4 组件(A) |保留 |
| trends/R3_suitability_area_trend_from_earliest_fossil | 适宜面积时间趋势 | → Fig 4 组件(B) |保留 |
| trends/R3_suitability_trend_from_earliest_fossil | 适宜趋势（另一版） | 探索版，疑似重复 |未使用 |
| trends/R3b_high_suitability_latitudinal_area_heatmap | 高适宜区纬度梯度热图 | → Fig 4 组件(C) |保留 |
| trends/R3b_latitudinal_suitability_area_gradient | 纬度梯度（同上另一版） | 探索版 |未使用 |
| trends/R3b_high_suitability_latitudinal_curves | 高适宜纬度曲线 | 探索版 |未使用 |
| trends/R3b_high_suitability_latitudinal_hemisphere_shift | 高适宜南北半球位移 | 探索版 |未使用 |
| trends/R3b_continuous_suitability_weighted_latitudinal_shift | 连续适宜度加权纬度位移 | 探索版 |未使用 |
| trends/R3b_continuous_suitability_weighted_latitudinal_hemisphere_shift | 同上南北半球版 | 探索版 |未使用 |
| trends/R3b_suitable_03_1_latitudinal_shift | 0.3–1 纬度位移 | 探索版 |未使用 |
| trends/R3b_suitable_03_1_latitudinal_hemisphere_shift | 0.3–1 南北半球位移 | 探索版 |未使用 |
| trends/R7_important_variable_background_contrast_heatmap | 重要变量高适宜 vs 背景对比热图 | → Fig 3 组件(B) |保留 |
| trends/R7_important_variable_background_contrast_heatmap_trends | 对比热图+趋势 | 探索版 |未使用 |
| trends/R7_important_variable_paleoclimate_background_heatmap | 重要变量古气候背景热图 | 备选 |未使用 |
| trends/R7_important_variable_paleoclimate_background_zscore | 重要变量 z-score 折线 | 探索版 |未使用 |
| fossil_comparison/R4_fossil_projection_panels | 化石点叠加投影图 | → Fig 5 组件(A) |保留 |
| fossil_comparison/R5_fossil_point_suitability | 化石点适宜度提取 | → Fig 5 组件(B) |保留 |
| fossil_comparison/R6_fossil_geographic_shift | 化石现代/古坐标位移 | 未用（与正文图重复） |未使用 |
| fossil_comparison/R8_fossil_significance_figures10 | 化石点适宜度显著性检验 | → Fig S11 |保留 |

## 组合图版（plates/）

| 文件（output/figures_paper/plates/…） | 内容 | v5.3 建议 | 你的决定 |
|---|---|---|---|
| Figure_1_M4_R7_variable_importance_and_background | M4 + R7 合图 | → Fig 3 |所有组合图版都不保留，我实际使用的是手动ps拼的图 |
| Figure_2_R3_R3b_temporal_latitudinal_dynamics | R3 + R3b 合图 | → Fig 4 的一部分 |同上 |
| Figure_3_R2b_all_period_suitability_panels | R2b 全时期面板 | → Fig 4 的一部分 |同上 |
| Figure_4_R4_R5_fossil_model_consistency | R4 + R5 合图 | → Fig 5 |同上 |
| Appendix_Figure_S1_M1_M2_M3_model_diagnostics | M1 + M2 + M3 合图 | → Fig S1？还是拆开 |同上 |

## 不在 figures_paper 的代码图

| 文章图 | 来源 | 说明 |
|---|---|---|
| Fig S4（Presence vs Background） | `output/08_projection/projection_modern_map.pdf` 的其中一页 | 由 `08_projection.R` 生成；如需独立文件，我可单独导出 |

决定：Fig. S4 的A 是output/08_projection/projection_modern_map.pdf 的第一页，B 是同文档第六页的Presence vs Background: Ensemble (EMwmean)这个图（也就是左边的那个图），然后手动加上了绿色虚线。

## 非代码图（不在本仓库，仅记录）

Fig 1（分布图，PSD）、Fig 2（形态照片）、Fig S5–S9（各国标本照片）-- 这些不来自本代码仓库。

---

**填表说明**：你只需在「你的决定」列标注。常见选项：
- `保留` -- 保留该图与生成它的脚本
- `删除` -- 删除该图输出，并移除/注释生成它的脚本段
- `合并到 XXX` -- 该图并入某张文章图（告诉我怎么并）
- 也可以直接说"组件图全删，只保留 plates/" 或反之。
