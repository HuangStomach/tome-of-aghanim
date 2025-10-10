import time
from ._cli import Cli

@Cli('创建所有数据实体')
def create_all():
    print('creating...')
    time.sleep(3)
    print('done.')

@Cli('拼接字符串')
def concat(a, b):
    print(a + b)
