import pandas
import numpy as np
np.set_printoptions(threshold=np.inf)

ccle = pandas.read_table("./CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", header = 0, skiprows = 2).to_numpy() # 读取对应的肿瘤细胞系
tcga_rpkm = pandas.read_table("./TCGA/ACC/HiSeqV2", header = 0).to_numpy() # 读取肿瘤基因组图谱中对应肿瘤的测序结果

genes = np.intersect1d(ccle[:, 1], tcga_rpkm[:, 0]) # 获取有交集的基因 len(genes) 18217
new_ccle = ccle[np.in1d(ccle[:, 1], genes)] # 过滤包含在交集基因中的数据
# 排除表达度低的基因, avg.RPKM > 1 in CCLE shape (12437, 1156)
ccle_gene_mat = np.array(list(filter(lambda cell_line: np.mean(cell_line) > 1, new_ccle[:, 2:]))) 

# 选择有>=20个细胞系的组织