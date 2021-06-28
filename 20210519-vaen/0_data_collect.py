import sys
import os
import time
from urllib import request

start_time = time.time()

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

def tcga():
    tcga = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']

    for t in tcga:
        url = 'https://bioinfo.uth.edu/VAEN/DATA/TCGA/{}/HiSeqV2'.format(t)
        filename = './TCGA/{}/HiSeqV2'.format(t)
        global start_time
        start_time = time.time()
        print("开始下载 {}".format(url))
        request.urlretrieve(url, filename, Schedule)
        print('\n下载完毕')

files = [
    ('./Data/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct', 'https://data.broadinstitute.org/ccle/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct'),
    ('./Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv', 'https://bioinfo.uth.edu/VAEN/DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv'),
    ('./Data/CCLE/DepMap-2018q3-celllines.csv', 'https://bioinfo.uth.edu/VAEN/DATA/CCLE/DepMap-2018q3-celllines.csv'),
    ('./Data/GDSC/v17.3_fitted_dose_response.txt', 'https://bioinfo.uth.edu/VAEN/DATA/GDSC/v17.3_fitted_dose_response.txt'),
    ('./Data/GDSC/Screened_Compounds.txt', 'https://bioinfo.uth.edu/VAEN/DATA/GDSC/Screened_Compounds.txt'),
    ('./Data/Match/drugs_match_2.txt', 'https://bioinfo.uth.edu/VAEN/DATA/drugs.match-2.txt'),
    ('./Data/Match/drugs_match.txt', 'https://bioinfo.uth.edu/VAEN/DATA/drugs.match.txt'),
    ('./Data/TCGA/*/HiSeqV2', 'dir')
]

str = ""
while True:
    print("需要准备的数据文件：")
    for i, item in enumerate(files):
        filename, url = item
        sign = "✅" if os.path.exists(filename) else "{}]".format(i)
        print("{} {}".format(sign, filename))

    str = input("请选择需要下载的数据文件, a为全部下载, q则退出: ");

    if str.lower() == "q":
        break
    elif str.lower() == 'a':
        tcga()
        for i, item in enumerate(files):
            filename, url = item
            start_time = time.time()
            request.urlretrieve(url, filename, Schedule)
    elif str.isdigit():
        index = int(str)
        if index < 0 or index >= len(files): continue
        filename, url = files[index]
        start_time = time.time()
        request.urlretrieve(url, filename, Schedule)
    else:
        continue
