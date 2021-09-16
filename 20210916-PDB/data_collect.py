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
    ['DB00622', 'P28223', ['6A93', '6A94', '6WGT', '6WH4', '6WHA']],
    ['DB01118', 'P07550', ['1GQ4', '2R4R', '2R4S', '2RH1', '3D4S', '3KJ6', '3NY8', '3NY9', '3NYA', '3P0G', '3PDS', '3SN6', '4GBR', '4LDE', '4LDL', '4LDO', '4QKX', '5D5A', '5D5B', '5D6L', '5JQH', '5X7D', '6E67', '6MXT', '6N48', '6OBA', '6PRZ', '6PS0', '6PS1', '6PS2', '6PS3', '6PS4', '6PS5', '6PS6']],
    ['DB00938', 'P08588', []],
    ['DB00734', 'P28222', ['4IAQ', '4IAR', '5V54', '7C61']],
    ['DB00834', 'P10275', ['1E3G', '1GS4', '1T5Z', '1T63', '1T65', '1XJ7', '1XOW', '1XQ3', '1Z95', '2AM9', '2AMA', '2AMB', '2AO6', '2AX6', '2AX7', '2AX8', '2AX9', '2AXA', '2HVC', '2OZ7', '2PIO', '2PIP', '2PIQ', '2PIR', '2PIT', '2PIU', '2PIV', '2PIW', '2PIX', '2PKL', '2PNU', '2Q7I', '2Q7J', '2Q7K', '2Q7L', '2YHD', '2YLO', '2YLP', '2YLQ', '2Z4J', '3B5R', '3B65', '3B66', '3B67', '3B68', '3BTR', '3L3X', '3L3Z', '3RLJ', '3RLL', '3V49', '3V4A', '3ZQT', '4HLW', '4K7A', '4OEA', '4OED', '4OEY', '4OEZ', '4OFR', '4OFU', '4OGH', '4OH5', '4OH6', '4OHA', '4OIL', '4OIU', '4OJ9', '4OJB', '4OK1', '4OKB', '4OKT', '4OKW', '4OKX', '4OLM', '4QL8', '5CJ6', '5JJM', '5T8E', '5T8J', '5V8Q', '5VO4']],
    ['DB00320', 'P08908', []],
    ['DB00308', 'Q14524', ['4DCK', '4DJC', '4JQ0', '4OVN', '5DBR', '6MUD']],
    ['DB00640', 'Q07343', ['1F0J', '1RO6', '1RO9', '1ROR', '1TB5', '1XLX', '1XLZ', '1XM4', '1XM6', '1XMU', '1XMY', '1XN0', '1XOS', '1XOT', '1Y2H', '1Y2J', '2CHM', '2QYL', '3D3P', '3FRG', '3G45', '3GWT', '3HC8', '3HDZ', '3HMV', '3KKT', '3LY2', '3O0J', '3O56', '3O57', '3W5E', '3WD9', '4KP6', '4MYQ', '4NW7', '4WZI', '4X0F', '5K6J', '5LAQ', '5OHJ', '6BOJ']],
    ['DB00193', 'P11229', ['5CXV', '6WJC']],
    ['DB00843', 'P14416', ['5AER', '6CM4', '6LUQ']],
    ['DB00482', 'Q9UPY5', []],
    ['DB01069', 'P08908', []],
    ['DB01151', 'P25100', []],
    ['DB00907', 'P20309', []],
    ['DB00933', 'P21728', []],
    ['DB00661', 'O95180', []],
    ['DB06216', 'Q9H3N8', []],
    ['DB00321', 'O43525', ['5J03']],
    ['DB00799', 'P48443', ['2GL8', '7A79']],
    ['DB00334', 'P41595', ['4IB4', '4NC3', '5TUD', '5TVN', '6DRX', '6DRY', '6DRZ', '6DS0']],
    ['DB01183', 'P62158', []],
    ['DB00441', 'P00374', ['1BOZ', '1DHF', '1DLR', '1DLS', '1DRF', '1HFP', '1HFQ', '1HFR', '1KMS', '1KMV', '1MVS', '1MVT', '1OHJ', '1OHK', '1PD8', '1PD9', '1PDB', '1S3U', '1S3V', '1S3W', '1U71', '1U72', '2C2S', '2C2T', '2DHF', '2W3A', '2W3B', '2W3M', '3EIG', '3F8Y', '3F8Z', '3F91', '3FS6', '3GHC', '3GHV', '3GHW', '3GI2', '3GYF', '3L3R', '3N0H', '3NTZ', '3NU0', '3NXO', '3NXR', '3NXT', '3NXV', '3NXX', '3NXY', '3NZD', '3OAF', '3S3V', '3S7A', '4DDR', '4G95', '4KAK', '4KBN', '4KD7', '4KEB', '4KFJ', '4M6J', '4M6K', '4M6L', '4QHV', '4QJC', '5HPB', '5HQY', '5HQZ', '5HSR', '5HSU', '5HT4', '5HT5', '5HUI', '5HVB', '5HVE', '6A7C', '6A7E', '6DAV', '6DE4', '6VCJ']],
    ['DB01611', 'Q9NR97', ['3W3G', '3W3J', '3W3K', '3W3L', '3W3M', '3W3N', '3WN4', '4QBZ', '4QC0', '4R07', '4R08', '4R09', '4R0A', '4R6A', '5AWA', '5AWB', '5AWC', '5AWD', '5AZ5', '5HDH', '5WYX', '5WYZ', '5Z14', '5Z15', '6KYA', '6TY5', '6V9U', '6WML', '6ZJZ', '7CRF']],
    ['DB00659', 'O60391', []],
]

headers = {
    ':authority': 'go.drugbank.com',
    ':method': 'GET',
    ':path': '/structures/small_molecule_drugs/DB00622.pdb',
    ':scheme': 'https',
    'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
    'accept-encoding': 'gzip, deflate, br',
    'accept-language': 'zh-CN,zh;q=0.9,ja;q=0.8,en;q=0.7,zh-TW;q=0.6',
    'cache-control': 'max-age=0',
    'cookie': '_ga=GA1.2.995934760.1631795845; _gid=GA1.2.1389548491.1631795845; _clck=1qbhj4s|1|eus|0; __hstc=49600953.d8d208984575045ae6165aea9e1bbe7c.1631795845885.1631795845885.1631795845885.1; hubspotutk=d8d208984575045ae6165aea9e1bbe7c; __hssrc=1; _gcl_au=1.1.1923849822.1631795848; _clsk=el7eqs|1631796420037|3|1|d.clarity.ms/collect; _omx_drug_bank_session=L0ZtdkxmS3VuYWtXWTRzVEVBQjR2dlpubko4eTRnWkVtSHg3VlQwbFozeitrYmRwQWdNWVh0SEl4Y0pVTnNabUh3UFN4cmVYUVVIVHgxbnp5MUJtKytvQ2dSdDB1TFVudURjYUp6YlNCZmE4Ym1kTlNJWXN5RHA3ZzZXblNoUyt0dFVTaFZVOG1mQWVjcVhJSlhvSWRnPT0tLVVneDFEdWVTcmZGcjhoWXhqWFIrYWc9PQ%3D%3D--fadadf2dbc5a631609c7bd57afadb1fd40757a05',
    'dnt': '1',
    'if-none-match': 'W/"38e2b3a13489db7ea7a80031a2c69828"',
    'sec-ch-ua': '"Google Chrome";v="93", " Not;A Brand";v="99", "Chromium";v="93"',
    'sec-ch-ua-mobile': '?0',
    'sec-ch-ua-platform': '"macOS"',
    'sec-fetch-dest': 'document',
    'sec-fetch-mode': 'navigate',
    'sec-fetch-site': 'none',
    'sec-fetch-user': '?1',
    'upgrade-insecure-requests': '1',
    'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/93.0.4577.63 Safari/537.36',
}

def download(url, path):
    try:
        req = request.Request(url, headers=headers)
        data = request.urlopen(req).read()
        with open(path, 'wb') as f:
            f.write(data)
            f.close()
    except Exception as e:
        print(str(e))

drug_url = 'https://go.drugbank.com/structures/small_molecule_drugs/{}.pdb'
protein_url = 'https://files.rcsb.org/download/{}.pdb'
for file in files:
    [drug, protein, identifier] = file
    start_time = time.time()
    print(drug_url.format(drug))
    download(drug_url.format(drug), './Data/{}.pdb'.format(drug))
    for item in identifier:
        start_time = time.time()
        print(protein_url.format(item))
        is_exists = os.path.exists('./Data/{}'.format(protein))
        if not is_exists:
            os.makedirs('./Data/{}'.format(protein)) 
        request.urlretrieve(protein_url.format(item), './Data/{}/{}.pdb'.format(protein, item), Schedule)
        time.sleep(2)
    print('')