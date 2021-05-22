import os
import pandas
import numpy as np
#np.set_printoptions(threshold=np.inf)

ccle_table = pandas.read_table("./CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", header = 0, skiprows = 2) # 读取对应的肿瘤细胞系
ccle = ccle_table.to_numpy()
tcga_rpkm = pandas.read_table("./TCGA/ACC/HiSeqV2", header = 0).to_numpy() # 读取肿瘤基因组图谱中对应肿瘤的测序结果

genes = np.intersect1d(ccle[:, 1], tcga_rpkm[:, 0]) # 获取有交集的基因 len(genes) 18217
new_ccle = ccle[np.in1d(ccle[:, 1], genes)] # 过滤包含在交集基因中的数据
# 排除表达度低的基因, avg.RPKM > 1 in CCLE shape (12437, 1156)
ccle_gene_mat = np.array(list(filter(lambda cell_line: np.mean(cell_line[2:]) > 1, new_ccle)))
ccle_gene_mat_names = ccle_gene_mat[:, :2] # 分离符合条件的基因名称
ccle_gene_mat = ccle_gene_mat[:, 2:].astype(np.float64) # 分离符合条件的基因数据

# 选择有 >=20 个细胞系的组织
columns = []
columns_count = {}
for col in ccle_table.columns[2:]:
    v = '_'.join(col.split('(')[0].split('_')[1:])
    columns_count.setdefault(v, 0)
    columns_count[v] += 1

for i, col in enumerate(ccle_table.columns[2:]):
    v = '_'.join(col.split('(')[0].split('_')[1:])
    if columns_count[v] >= 20: columns.append(i)

# 对ccle数据进行log2转换
ccle_gene_mat = ccle_gene_mat[:, columns] # shape (12437, 1100)
log2_ccle_gene_mat = np.log2(ccle_gene_mat + 1)

# 选择方差变化较大的前一半样本
variance = log2_ccle_gene_mat.var(1)
mid = np.median(variance)
mad5000_ccle_gene_mat = []
mad_genes = []
for i, cell in enumerate(log2_ccle_gene_mat):
    if cell.var() > mid: 
        # 同时筛选数据和基因名称
        mad5000_ccle_gene_mat.append(cell)
        mad_genes.append(ccle_gene_mat_names[i])

mad5000_ccle_gene_mat = np.array(mad5000_ccle_gene_mat) # shape (6219, 1100)
mad_genes = np.array(mad_genes)

# 按照TCGA目录下的癌症类型进行数据整理
cancer_types = os.listdir('./TCGA')
cancer_types.sort()

rpkm_mat = np.array([])
cur_genes = []
for i, cancer in enumerate(cancer_types):
    tcga_rpkm = pandas.read_table("./TCGA/{}/HiSeqV2".format(cancer), header = 0).to_numpy()
    # 排除行和为0和错误的数据的基因
    tcga_rpkm = np.array([row for row in tcga_rpkm if sum(row[1:]) != np.NaN or sum(row[1:]) != 0])
    tcga_rpkm = tcga_rpkm[np.in1d(tcga_rpkm[:, 0], mad_genes[:, 1])] # 筛选基因有交集的数据

    if i == 0: cur_genes = tcga_rpkm[:, 0]
    else: cur_genes = np.intersect1d(cur_genes, tcga_rpkm[:, 0]) # 不停获取基因的交集

    t_tcga_rpkm = np.transpose(tcga_rpkm)
    if len(rpkm_mat) == 0: rpkm_mat = t_tcga_rpkm[:, np.in1d(tcga_rpkm[:, 0], cur_genes)]
    else:
        rpkm_mat = np.row_stack(
            (rpkm_mat[:, np.in1d(rpkm_mat[0, :], cur_genes)], t_tcga_rpkm[1:, np.in1d(tcga_rpkm[:, 0], cur_genes)])
        )
    print(cancer, "\t", len(tcga_rpkm[0]), len(rpkm_mat[0]), len(rpkm_mat))

#ACC     79 6203 79
#BLCA    426 6203 505
#BRCA    1218 6203 1723
#CESC    308 6203 2031
#CHOL    45 6203 2076
#COAD    329 6203 2405
#DLBC    48 6203 2453
#ESCA    196 6203 2649
#GBM     172 6203 2821
#HNSC    566 6203 3387
#KICH    91 6203 3478
#KIRC    606 6203 4084
#KIRP    323 6203 4407
#LAML    173 6203 4580
#LGG     530 6203 5110
#LIHC    423 6203 5533
#LUAD    576 6203 6109
#LUSC    553 6203 6662
#MESO    87 6203 6749
#OV      308 6203 7057
#PAAD    183 6203 7240
#PCPG    187 6203 7427
#PRAD    550 6203 7977
#READ    105 6203 8082
#SARC    265 6203 8347
#SKCM    474 6203 8821
#STAD    450 6203 9271
#TGCT    156 6203 9427
#THCA    572 6203 9999
#THYM    122 6203 10121
#UCEC    201 6203 10322
#UCS     57 6203 10379
#UVM     80 6203 10459
