### 下载训练数据

``` bash
# Data/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct
curl -x host:port -O 'https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct'
# Data/TCGA/ACC/HiSeqV2
curl -x host:port -O 'https://bioinfo.uth.edu/VAEN/DATA/TCGA/ACC/HiSeqV2?csrt=15326058242117628311'
```

### 数据预处理 以及脚本所需数据

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