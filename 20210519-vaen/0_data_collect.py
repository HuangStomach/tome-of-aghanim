import sys
import time
from urllib import request

def Schedule(blocknum, blocksize, totalsize):
    speed = (blocknum * blocksize) / (time.time() - start_time)
    speed_str = " 速度: %s" % format_size(speed)
    recv_size = blocknum * blocksize
     
    # 设置下载进度条
    f = sys.stdout
    pervent = min(1.0, recv_size / totalsize)
    percent_str = "%.2f%%" % (pervent * 100)
    n = round(pervent * 50)
    s = ('#' * n).ljust(50, '-')
    f.write(percent_str.ljust(8, ' ') + '[' + s + ']' + speed_str)
    f.flush()
    f.write('\r')
 
# 字节bytes转化K\M\G
def format_size(bytes):
    try:
        bytes = float(bytes)
        kb = bytes / 1024
    except:
        print("传入的字节格式不对")
        return "Error"
    
    if kb >= 1024:
        M = kb / 1024
        if M >= 1024:
            G = M / 1024
            return "%.3fG" % (G)
        else:
            return "%.3fM" % (M)
    else:
        return "%.3fK" % (kb)
 
start_time = time.time()
filename = './CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct'
url = 'https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct'
request.urlretrieve(url, filename, Schedule)

start_time = time.time()
filename = './CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv'
url = 'https://bioinfo.uth.edu/VAEN/DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv'
request.urlretrieve(url, filename, Schedule)

tcga = ['BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

for t in tcga:
    url = 'https://bioinfo.uth.edu/VAEN/DATA/TCGA/{}/HiSeqV2'.format(t)
    filename = './TCGA/{}/HiSeqV2'.format(t)
    start_time = time.time()
    print("开始下载 {}".format(url))
    request.urlretrieve(url, filename, Schedule)
    print('\n下载完毕')