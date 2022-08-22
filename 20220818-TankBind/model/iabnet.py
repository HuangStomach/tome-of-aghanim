import torch
from torch_geometric.utils import to_dense_batch
from torch import nn
from torch.nn import Linear
import torch.nn as nn

from gvp.gvp_embedding import GVP_embedding
from gnn import GNN
from gin import GIN
from triangle import TriangleProteinToCompound, TriangleSelfAttentionRowWise, Transition

class IaBNet_with_affinity(torch.nn.Module):
    def __init__(self, 
        hidden_channels=128, embedding_channels=128, c=128, mode=0, 
        protein_embed_mode=1, compound_embed_mode=1, n_trigonometry_module_stack=5, 
        protein_bin_max=30, readout_mode=2):

        super().__init__()
        self.layernorm = torch.nn.LayerNorm(embedding_channels)
        self.protein_bin_max = protein_bin_max
        self.mode = mode
        self.protein_embed_mode = protein_embed_mode
        self.compound_embed_mode = compound_embed_mode
        self.n_trigonometry_module_stack = n_trigonometry_module_stack
        self.readout_mode = readout_mode

        if protein_embed_mode == 0:
            self.conv_protein = GNN(hidden_channels, embedding_channels) # GraphSAGE 的特征聚合方法
        if protein_embed_mode == 1:
            self.conv_protein = GVP_embedding(
                (6, 3), (embedding_channels, 16), (32, 1), (32, 1), seq_in=True
            )

        if compound_embed_mode == 0:
            self.conv_compound = GNN(hidden_channels, embedding_channels) # GraphSAGE 的特征聚合方法
        elif compound_embed_mode == 1:
            self.conv_compound = GIN(
                input_dim = 56, hidden_dims = [128,56,embedding_channels], 
                edge_input_dim = 19, concat_hidden = False
            ) # 药物特征使用的torch_drug走GIN

        if mode == 0:
            self.protein_pair_embedding = Linear(16, c)
            self.compound_pair_embedding = Linear(16, c)
            self.protein_to_compound_list = []
            self.protein_to_compound_list = nn.ModuleList([
                TriangleProteinToCompound(embedding_channels, c) for _ in range(n_trigonometry_module_stack)
            ])
            self.triangle_self_attention_list = nn.ModuleList([
                TriangleSelfAttentionRowWise(embedding_channels) for _ in range(n_trigonometry_module_stack)
            ])
            self.tranistion = Transition(embedding_channels, n=4)

        self.linear = Linear(embedding_channels, 1)
        self.linear_energy = Linear(embedding_channels, 1)
        if readout_mode == 2: self.gate_linear = Linear(embedding_channels, 1)
        # self.gate_linear = Linear(embedding_channels, 1)
        self.bias = torch.nn.Parameter(torch.ones(1))
        self.leaky = torch.nn.LeakyReLU()
        self.dropout = nn.Dropout2d(p=0.25)

    def forward(self, data):
        # 首先利用GNN来聚合蛋白和药物的特征
        # 数据应该是一个异构网络
        if self.protein_embed_mode == 0:
            x = data['protein'].x.float()
            edge_index = data[("protein", "p2p", "protein")].edge_index
            protein_batch = data['protein'].batch
            protein_out = self.conv_protein(x, edge_index)
        elif self.protein_embed_mode == 1:
            nodes = (
                data['protein']['node_s'], 
                data['protein']['node_v']
            )
            edges = (
                data[("protein", "p2p", "protein")]["edge_s"], 
                data[("protein", "p2p", "protein")]["edge_v"]
            )
            protein_batch = data['protein'].batch
            protein_out = self.conv_protein(nodes, data[("protein", "p2p", "protein")]["edge_index"], edges, data.seq) # 使用GVP提取蛋白特征，利用旋转不变性提取标量特征

        compound_x = data['compound'].x.float()
        compound_batch = data['compound'].batch
        if self.compound_embed_mode == 0:
            compound_edge_index = data[("compound", "c2c", "compound")].edge_index.T
            compound_out = self.conv_compound(compound_x, compound_edge_index)
        elif self.compound_embed_mode == 1:
            compound_edge_feature = data[("compound", "c2c", "compound")].edge_attr
            edge_weight = data[("compound", "c2c", "compound")].edge_weight
            compound_out = self.conv_compound(
                compound_edge_index, edge_weight ,compound_edge_feature, compound_x.shape[0], compound_x
            )['node_feature'] # 药物直接使用gin来提取特征

        # protein_batch version could further process b matrix. better than for loop.
        # protein_out_batched of shape b, n, c
        # 图节点做batch的时候，会稀疏，这个方法把节点重新聚合成batch大小的矩阵
        # 并一起返回一个mask来指明哪些节点是fake的，方便后续运算
        protein_out_batched, protein_out_mask = to_dense_batch(protein_out, protein_batch)
        compound_out_batched, compound_out_mask = to_dense_batch(compound_out, compound_batch)

        node_xyz = data.node_xyz # 不知道具体是什么 猜测是蛋白质的坐标
        p_coords_batched, _ = to_dense_batch(node_xyz, protein_batch)
        # c_coords_batched, c_coords_mask = to_dense_batch(coords, compound_batch)

        protein_pair = self.get_pair_dis_one_hot(
            p_coords_batched, bin_size=2, bin_min=-1, bin_max=self.protein_bin_max
        ) # 蛋白质之间距离的one-hot编码

        # compound_pair = get_pair_dis_one_hot(c_coords_batched, bin_size=1, bin_min=-0.5, bin_max=15)
        compound_pair_batched, _ = to_dense_batch(data.compound_pair, data.compound_pair_batch)
        batch_n = compound_pair_batched.shape[0]
        max_compound_size_square = compound_pair_batched.shape[1]
        max_compound_size = int(max_compound_size_square**0.5)
        assert (max_compound_size**2 - max_compound_size_square)**2 < 1e-4
        compound_pair = torch.zeros((batch_n, max_compound_size, max_compound_size, 16)).to(data.compound_pair.device)

        for i in range(batch_n):
            one = compound_pair_batched[i]
            compound_size_square = (data.compound_pair_batch == i).sum()
            compound_size = int(compound_size_square**0.5)
            compound_pair[i,:compound_size, :compound_size] = one[:compound_size_square].reshape(
                (compound_size, compound_size, -1)
            )
        
        protein_pair = self.protein_pair_embedding(protein_pair.float())
        compound_pair = self.compound_pair_embedding(compound_pair.float()) # 根据batch获得药物对儿
        # 对以上的蛋白和药物数据进行简单的线性变换embedding
        # b = torch.einsum("bik,bjk->bij", protein_out_batched, compound_out_batched).flatten()

        # 对药物和蛋白的特征进行标准化
        protein_out_batched = self.layernorm(protein_out_batched)
        compound_out_batched = self.layernorm(compound_out_batched)
        # z of shape, b, protein_length, compound_length, channels.
        # 对药物和蛋白的特征进行矩阵乘
        z = torch.einsum("bik,bjk->bijk", protein_out_batched, compound_out_batched)
        z_mask = torch.einsum("bi,bj->bij", protein_out_mask, compound_out_mask)
        # z = z * z_mask.unsqueeze(-1)
        # print(protein_pair.shape, compound_pair.shape, b.shape)

        if self.mode == 0:
            for _ in range(1):
                for i_module in range(self.n_trigonometry_module_stack):
                    z = z + self.dropout(
                        self.protein_to_compound_list[i_module](
                            z, protein_pair, compound_pair, z_mask.unsqueeze(-1)
                        )
                    )
                    z = z + self.dropout(
                        self.triangle_self_attention_list[i_module](z, z_mask)
                    )
                    z = self.tranistion(z) # 普通的线性变换
        # batch_dim = z.shape[0]

        b = self.linear(z).squeeze(-1)
        y_pred = b[z_mask]
        y_pred = y_pred.sigmoid() * 10   # normalize to 0 to 10.
        if self.readout_mode == 0:
            pair_energy = self.linear_energy(z).squeeze(-1) * z_mask
            affinity_pred = self.leaky(self.bias + ((pair_energy).sum(axis=(-1, -2))))
        if self.readout_mode == 1:
            # valid_interaction_z = (z * z_mask.unsqueeze(-1)).mean(axis=(1, 2))
            valid_interaction_z = (z * z_mask.unsqueeze(-1)).sum(axis=(1, 2)) / z_mask.sum(axis=(1, 2)).unsqueeze(-1)
            affinity_pred = self.linear_energy(valid_interaction_z).squeeze(-1)
            # print("z shape", z.shape, "z_mask shape", z_mask.shape, "valid_interaction_z shape", valid_interaction_z.shape, "affinity_pred shape", affinity_pred.shape)
        if self.readout_mode == 2:
            pair_energy = (self.gate_linear(z).sigmoid() * self.linear_energy(z)).squeeze(-1) * z_mask
            affinity_pred = self.leaky(self.bias + ((pair_energy).sum(axis=(-1, -2))))
        return y_pred, affinity_pred # 后面就是简单的线性变换来得到输出

    def get_pair_dis_one_hot(d, bin_size=2, bin_min=-1, bin_max=30):
        # without compute_mode='donot_use_mm_for_euclid_dist' could lead to wrong result.
        # cdist
        # 计算两组行向量集合中每一对之间的p-norm距离。 这里是默认的欧式距离
        pair_dis = torch.cdist(d, d, compute_mode='donot_use_mm_for_euclid_dist')
        pair_dis[pair_dis>bin_max] = bin_max # 最大距离为30
        pair_dis_bin_index = torch.div(pair_dis - bin_min, bin_size, rounding_mode='floor').long() # 类似于标准化
        pair_dis_one_hot = torch.nn.functional.one_hot(pair_dis_bin_index, num_classes=16)
        # 反回表示距离的one-hot向量
        return pair_dis_one_hot
