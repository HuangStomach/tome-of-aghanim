## 启动命令

```bash
docker run --gpus -itd -v <ORIGINAL PATH>:/data <IMAGE ID>
```

## 进入容器

```bash
docker exec -it <CONTAINER ID> bash
```

## 数据准备

宿主机路径`<ORIGINAL PATH>`映射到容器内`/data`目录下，输入文件需命名为`input.fasta`，在运行模型之后，输出会生成在同一目录下命名为`output.fasta`

```python
input_data_dir = "/data/input.fasta"
output_data_dir = "/data/output.fasta"
```

## 运行模型

```bash
$ python 2_submit2_self_addid.py
```

