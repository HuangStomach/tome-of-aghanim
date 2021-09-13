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

files = [
    ('./Data/adj_E.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/e_admat_dgc.txt'),
    ('./Data/adj_IC.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ic_admat_dgc.txt'),
    ('./Data/adj_GPCR.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/gpcr_admat_dgc.txt'),
    ('./Data/adj_NR.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/nr_admat_dgc.txt'),
    ('./Data/drug_E.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/e_simmat_dc.txt'),
    ('./Data/drug_IC.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ic_simmat_dc.txt'),
    ('./Data/drug_GPCR.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/gpcr_simmat_dc.txt'),
    ('./Data/drug_NR.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/nr_simmat_dc.txt'),
    ('./Data/target_E.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/e_simmat_dg.txt'),
    ('./Data/target_IC.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/ic_simmat_dg.txt'),
    ('./Data/target_GPCR.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/gpcr_simmat_dg.txt'),
    ('./Data/target_NR.txt', 'http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/nr_simmat_dg.txt'),
]

str_in = ""
while True:
    print("需要准备的数据文件：")
    for i, item in enumerate(files):
        filename, url = item
        sign = "✅" if os.path.exists(filename) else "{}]".format(i)
        print("{} {}".format(sign, filename))

    str_in = input("请选择需要下载的数据文件, a为全部下载, q则退出: ");

    if str_in.lower() == "q":
        break
    elif str_in.lower() == 'a':
        for i, item in enumerate(files):
            filename, url = item
            start_time = time.time()
            request.urlretrieve(url, filename, Schedule)
            print('')
    elif str_in.isdigit():
        index = int(str_in)
        if index < 0 or index >= len(files): continue
        filename, url = files[index]
        if isinstance(url, str) == False :
            url()
            continue
        start_time = time.time()
        request.urlretrieve(url, filename, Schedule)
    else:
        continue