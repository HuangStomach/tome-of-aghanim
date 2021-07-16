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
> （？）中问题，是否无须VAE，只使用神经网络对数据进行无监督学习降维？
> 说是tsne降维看不出来，但是是使用组织进行降维而不是使用癌症类型，TCGA癌症类型降维效果明显
* 再构建elastic net模型，使用CCLE和GDSC进行训练，保存为 `CCLE\GDSC.model.list.Rdata`。
> 其中使用`glmnet`包，在`cv.glmnet`函数中设置`alpha`为0-1之间的数即可使用elastic net进行回归
> 还有疑问，这里进行训练的时候没有区分癌症，而是对某种药物的药物反应数据全部进行训练，但是预测的时候又按照癌症进行了区分
* 之后选择其中表现最好的模型对tcga数据进行预测分析，保存为 `VAEN_CCLE\GDSC.*.pred_TCGA.txt`。
* 对ccle数据进行预测分析，保存为 `VAEN_CCLE\GDSC.A.pred_CCLE.full.txt`。
* 将训练模型过程中对训练数据的自我预测保存为 `VAEN_CCLE\GDSC.A.pred_CCLE\GDSC.txt`。
> 这里有疑问，实际上训练过程使用的就是CCLE_latent，也就是pred_CCLE是训练latent时的自我预测，然后pred_CCLE.full是用这个模型单独对训练数据又预测一遍，需要再研究。
* 针对训练过的所有药物，将每个药物表现最好的模型存入 `dr.CCLE\GDSC.A.models.RData`
* 再将CCLE全部数据和固体数据（？）进行整合，存入 `VAEN_CCLE.MIX.pred_CCLE\TCGA.*.txt`。

## 改进计划

* 不使用vae，直接使用auto-encoder进行降维度。
* 用tcga数据训练分类，之后匹配ccle数据再进行标注。(单纯分类，或者domain-adversarial training)

## TCGA癌症中英文对照

|Cohort|英文名称|中文名称|
| ---- | ---- | ---- |
|ACC|Adrenocortical carcinoma|肾上腺皮质癌|
|BLCA|Bladder Urothelial Carcinoma|膀胱尿路上皮癌|
|BRCA|Breast invasive carcinoma|乳腺浸润癌|
|CESC|Cervical squamous cell carcinoma and endocervical adenocarcinoma|宫颈鳞癌和腺癌|
|CHOL|Cholangiocarcinoma|胆管癌|
|COAD|Colon adenocarcinoma|结肠癌|
|COADREAD|Colon adenocarcinoma/Rectum adenocarcinoma Esophageal carcinoma|结直肠癌|
|DLBC|Lymphoid Neoplasm Diffuse Large B-cell Lymphoma|弥漫性大B细胞淋巴瘤|
|ESCA|Esophageal carcinoma|食管癌|
|FPPP|FFPE Pilot Phase II|FFPE试点二期|
|GBM|Glioblastoma multiforme|多形成性胶质细胞瘤|
|GBMLGG|Glioma|胶质瘤|
|HNSC|Head and Neck squamous cell carcinoma|头颈鳞状细胞癌|
|KICH|Kidney Chromophobe|肾嫌色细胞癌|
|KIPAN|Pan-kidney cohort (KICH+KIRC+KIRP)|混合肾癌|
|KIRC|Kidney renal clear cell carcinoma|肾透明细胞癌|
|KIRP|Kidney renal papillary cell carcinoma|肾乳头状细胞癌|
|LAML|Acute Myeloid Leukemia|急性髓细胞样白血病|
|LGG|Brain Lower Grade Glioma|脑低级别胶质瘤|
|LIHC|Liver hepatocellular carcinoma|肝细胞肝癌|
|LUAD|Lung adenocarcinoma|肺腺癌|
|LUSC|Lung squamous cell carcinoma|肺鳞癌|
|MESO|Mesothelioma|间皮瘤|
|OV|Ovarian serous cystadenocarcinoma|卵巢浆液性囊腺癌|
|PAAD|Pancreatic adenocarcinoma|胰腺癌|
|PCPG|Pheochromocytoma and Paraganglioma|嗜铬细胞瘤和副神经节瘤|
|PRAD|Prostate adenocarcinoma|前列腺癌|
|READ|Rectum adenocarcinoma|直肠腺癌|
|SARC|Sarcoma|肉瘤|
|SKCM|Skin Cutaneous Melanoma|皮肤黑色素瘤|
|STAD|Stomach adenocarcinoma|胃癌|
|STES|Stomach and Esophageal carcinoma|胃和食管癌|
|TGCT|Testicular Germ Cell Tumors|睾丸癌|
|THCA|Thyroid carcinoma|甲状腺癌|
|THYM|Thymoma|胸腺癌|
|UCEC|Uterine Corpus Endometrial Carcinoma|子宫内膜癌|
|UCS|Uterine Carcinosarcoma|子宫肉瘤|
|UVM|Uveal Melanoma|葡萄膜黑色素瘤|