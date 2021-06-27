
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

|sample|TCGA-OR-A5LC-01| ... |
| ---- | ---- | ---- |
| 基因ID | 基因名称 | ... |

``` R
# Z01
setwd("/path/to/work/Z01.ReLU/")

system("cp /work/Z01/CCLE.4VAE.Z01.tsv /work/Z01.ReLU/")
system("cp /work/Z01/TCGA.4VAE.Z01.tsv /work/Z01.ReLU/")
dir.create("/work/Z01.ReLU/result")
dir.create("/work/Z01.ReLU/result.EN/")
dir.create("/work/Z01.ReLU/result.EN/dr.CCLE")
dir.create("/work/Z01.ReLU/result.EN/dr.GDSC")

for(b in 1:100) {
	cmd = paste("python3 run.VAE_100_0005_100_100.ReLU.py ", b, 
	" /work/Z01.ReLU/ CCLE.4VAE.Z01.tsv TCGA.4VAE.Z01.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}

# ZS
setwd("/path/to/work/ZS.ReLU/")

system("cp /work/ZS/CCLE.4VAE.ZS.tsv /work/ZS.ReLU/")
system("cp /work/ZS/TCGA.4VAE.ZS.tsv /work/ZS.ReLU/")
dir.create("/work/ZS.ReLU/result")
dir.create("/work/ZS.ReLU/result.EN/")
dir.create("/work/ZS.ReLU/result.EN/dr.CCLE")
dir.create("/work/ZS.ReLU/result.EN/dr.GDSC")

for(b in 1:100){
	cmd = paste("python3 run.VAE_100_0005_100_100.ReLU.py ", b, 
	" /work/ZS.ReLU/ CCLE.4VAE.ZS.tsv TCGA.4VAE.ZS.tsv CCLE.latent.tsv CCLE.weight.tsv TCGA.latent.tsv encoder.hdf5 decoder.hdf5", sep="" )
	system(cmd)
}
```