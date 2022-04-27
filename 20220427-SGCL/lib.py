import torch

# tuple of node-tensors
def node_tensors(matrix, threshold = 0):
    src = []
    dst = []
    (r, c) = matrix.shape
    for i in range(r):
        for j in range(c):
            if matrix[i][j] > threshold:
                src.append(i)
                dst.append(j)
    return torch.tensor(src), torch.tensor(dst)
