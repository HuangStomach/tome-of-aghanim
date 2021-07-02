
# VAEN 复现

> 更多专注于 rank+sigmoid 模式

## 数据文件说明

### CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct

不同癌症细胞类型中不同基因的表达（RPKM）

> RPKM是Reads Per Kilobase per Million mapped reads的缩写，代表每百万reads中来自于某基因每千碱基长度的reads数。RPKM是将map到基因的read数除以map到基因组上的所有read数(以million为单位)与RNA的长度(以KB为单位)。
>
> RNA-seq是二代测序技术中用来表示基因表达量或丰度的方法。在衡量基因表达量时，若是单纯以map到的read数来计算基因的表达量，在统计上是不合理的。因为在随机抽样的情况下，序列较长的基因被抽到的机率本来就会比序列短的基因较高，如此一来，序列长的基因永远会被认为表达量较高，而错估基因真正的表现量，所以Ali Mortazavi等人在2008年提出以RPKM在估计基因的表现量。

|Name|Description|22RV1_PROSTATE (ACH-000956)| ... |
| ---- | ---- | ---- | ---- |
| 基因ID | 基因名称 | 该癌细胞下的表达 | ... |

### HiSeqV2

每个HiSeqV2文件所在的目录代表一种癌症类型，描述不同样本对不同基因的表达

> 14-15位为2位数字，01-09表示肿瘤样本，10-16表示正常对照样本
> [sample-type-codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)

|sample|TCGA-OR-A5LC-01| ... |
| ---- | ---- | ---- |
| 基因ID | 基因名称 | ... |

### CCLE_NP24.2009_Drug_data_2015.02.24.csv

24种抗癌症药物对不同癌症的表现

### v17.3_fitted_dose_response.txt

Genomics of Drug Sensitivity in Cancer (GDSC) 肿瘤药物敏感性基因组学数据库 提供的癌症药物反应数据

### *.model.list.RData

使用EN训练后的结果，包含训练模型以及许多训练结果和数据。

``` R
res_list <- list()
res_list[["model_summary"]] <- model_summary # 包含超参数与损失的数据集合
res_list[["weight.mat"]] <- weight.mat # 权重中的矩阵
res_list[["model"]] <- fit # 训练模型 可用来直接预测
# res_list[["original_Y"]] <- Y
# res_list[["self_pred"]] <- adj_expr_pred
res_list[['ys']] <- ["original_Y", "self_pred"] # 包含真值与预测值的list
```

## 具体步骤

* 先对数据进行预处理，只采用CCLE和TCGA数据均有交集，样本数量足够大且活跃的数据，保存为 `V15.CCLE\TCGA.4VAE.*.tsv`。
* 使用VAE（？）对数据进行去噪和降维，保存为 `CCLE\TCGA_latent.tsv`。
* 再构建elastic net模型，使用CCLE和GDSC进行训练，保存为 `CCLE\GDSC.model.list.Rdata`。
> 其中使用`glmnet`包，在`cv.glmnet`函数中设置`alpha`为0-1之间的数即可使用elastic net进行回归
* 之后选择其中表现最好的模型对tcga数据进行预测分析，保存为 `VAEN_CCLE\GDSC.*.pred_TCGA.txt`。
* 对ccle数据进行预测分析，保存为 `VAEN_CCLE\GDSC.A.pred_CCLE.full.txt`。
* 将训练模型过程中对训练数据的预测保存为 `VAEN_CCLE\GDSC.A.pred_CCLE\GDSC.txt`。
* 针对训练过的所有药物，将每个药物表现最好的模型存入 `dr.CCLE\GDSC.A.models.RData`
* 再将CCLE全部数据和固体数据（？）进行整合，存入 `VAEN_CCLE.MIX.pred_CCLE\TCGA.*.txt`。
