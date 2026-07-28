# Variable Selection Methods / 环境变量筛选方法

## 中文

为降低气候变量间多重共线性对模型的影响，本研究首先对 19 个 BioClim 变量进行 Pearson 相关性分析，并结合变量间相关性、生态学意义和信息冗余程度进行初步筛选。初步选取 8 个候选变量，分别为 Bio2、Bio3、Bio4、Bio8、Bio9、Bio15、Bio18 和 Bio19。

随后采用方差膨胀因子（Variance Inflation Factor, VIF）进行逐步筛选。初始阈值设为 10，在该阈值下 Bio3 被剔除，Bio19 也会被排除。然而，Bio19 表征最冷季度降水量，对解释寒冷季节水分限制及物种潜在分布边界具有明确生态学意义，且其 VIF 值接近 10。考虑到 Bio19 的生态学意义，且最终变量间未出现高相关性，本研究将 VIF 阈值适度放宽至 11，以保留 Bio19。

最终用于建模的环境变量为 7 个：Bio2、Bio4、Bio8、Bio9、Bio15、Bio18 和 Bio19。

## English

To reduce the influence of multicollinearity among climatic predictors on model performance, we first performed Pearson correlation analysis on the 19 BioClim variables and initially selected variables by considering pairwise correlation, ecological relevance, and information redundancy. Eight candidate variables were initially selected: Bio2, Bio3, Bio4, Bio8, Bio9, Bio15, Bio18, and Bio19.

Variance inflation factor (VIF) analysis was then conducted using a stepwise filtering procedure. The initial threshold was set to 10, under which Bio3 was removed and Bio19 would also be excluded. However, Bio19 represents precipitation of the coldest quarter, which is ecologically meaningful for explaining cold-season water limitation and potential distribution limits. Considering the ecological relevance of Bio19 and the absence of high correlation among the final predictors, the VIF threshold was moderately relaxed to 11 to retain Bio19.

The final set of environmental variables used for modelling included seven variables: Bio2, Bio4, Bio8, Bio9, Bio15, Bio18, and Bio19.
