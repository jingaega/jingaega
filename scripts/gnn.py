"""GCN with hybrid SAGPool/statistical-prior pooling (rule 10 reference) and a matched MLP.

Graphs are built on the fold's training data only. Fixed training budget: 200 epochs, no
early stopping, no checkpoint selection (rule 4).
"""
import functools
import os

import networkx as nx
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from common import PROC

EPOCHS, LR, WD, BATCH, HIDDEN = 200, 1e-3, 5e-4, 16, 64
N_EDGES_COEXPR = 2500          # top |r| edges among the 500 genes -> mean degree 10
LAMBDA, POOL_RATIO = 0.5, 0.5
torch.set_num_threads(1)

_STRING = None


# ---------------------------------------------------------------- graphs (edge lists over gene indices)
def g_coexpr(f, rng):
    r = np.corrcoef(f.Xtr.T)
    iu = np.triu_indices_from(r, 1)
    w = np.abs(r[iu])
    top = np.argsort(-w, kind="mergesort")[:N_EDGES_COEXPR]
    return np.c_[iu[0][top], iu[1][top]]


def g_string(f, rng):
    global _STRING
    if _STRING is None:
        _STRING = pd.read_csv(os.path.join(PROC, "string_v12_700.tsv"), sep="\t")
    pos = {g: i for i, g in enumerate(f.genes)}
    s = _STRING[_STRING.a.isin(pos) & _STRING.b.isin(pos)]
    return np.c_[s.a.map(pos).values, s.b.map(pos).values].astype(int)


def g_none(f, rng):
    return np.zeros((0, 2), dtype=int)


def randomise_uniform(edges, n, rng):
    """Same edge count, uniform random topology (C2)."""
    m, seen = len(edges), set()
    while len(seen) < m:
        a, b = rng.integers(0, n, 2)
        if a != b:
            seen.add((min(a, b), max(a, b)))
    return np.array(sorted(seen), dtype=int).reshape(-1, 2)


def randomise_degree(edges, n, rng):
    """Degree-preserving double-edge swaps (C2b, configuration-model null as in Brouard 2024)."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(map(tuple, edges))
    m = G.number_of_edges()
    if m >= 2:
        nx.double_edge_swap(G, nswap=10 * m, max_tries=1000 * m, seed=int(rng.integers(1 << 31)))
    return np.array(list(G.edges()), dtype=int).reshape(-1, 2)


def norm_adj(edges, n):
    """Sparse D^-1/2 (A + I) D^-1/2."""
    e = np.r_[edges, edges[:, ::-1], np.c_[np.arange(n), np.arange(n)]] if len(edges) else np.c_[np.arange(n), np.arange(n)]
    deg = np.bincount(e[:, 0], minlength=n).astype(np.float32)
    w = 1.0 / np.sqrt(deg[e[:, 0]] * deg[e[:, 1]])
    return torch.sparse_coo_tensor(torch.tensor(e.T), torch.tensor(w, dtype=torch.float32), (n, n)).coalesce()


# ---------------------------------------------------------------- model
class GCNPool(nn.Module):
    def __init__(self, A, prior, n_classes):
        super().__init__()
        self.A = A
        self.prior = prior                                     # (N,) in [0,1], from limma p (training)
        self.W = nn.ModuleList([nn.Linear(1, HIDDEN), nn.Linear(HIDDEN, HIDDEN), nn.Linear(HIDDEN, HIDDEN)])
        self.att = nn.Linear(HIDDEN, 1)
        self.out = nn.Linear(2 * HIDDEN, n_classes)

    def prop(self, H):                                         # H: (B, N, F) -> A @ H per sample
        B, N, F = H.shape
        return torch.sparse.mm(self.A, H.permute(1, 0, 2).reshape(N, B * F)).reshape(N, B, F).permute(1, 0, 2)

    def forward(self, x):                                      # x: (B, N)
        H = x.unsqueeze(-1)
        for lin in self.W:
            H = torch.relu(self.prop(lin(H)))
        z = torch.sigmoid(self.prop(self.att(H))).squeeze(-1)  # SAGPool score (GAMB-GNN eq. 13)
        zmin, zmax = z.min(1, keepdim=True)[0], z.max(1, keepdim=True)[0]
        z = (z - zmin) / (zmax - zmin + 1e-8)                  # eq. 15
        q = LAMBDA * self.prior.unsqueeze(0) + (1 - LAMBDA) * z  # eq. 17
        k = max(1, int(POOL_RATIO * x.shape[1]))
        idx = q.topk(k, dim=1).indices
        Hs = torch.gather(H, 1, idx.unsqueeze(-1).expand(-1, -1, H.shape[-1]))
        Hs = Hs * torch.gather(q, 1, idx).unsqueeze(-1)        # h_v <- Q_v h_v
        return self.out(torch.cat([Hs.mean(1), Hs.max(1)[0]], 1))


class MLP(nn.Module):
    def __init__(self, n_in, n_classes):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(n_in, HIDDEN), nn.ReLU(), nn.Dropout(0.3),
                                 nn.Linear(HIDDEN, HIDDEN), nn.ReLU(), nn.Dropout(0.3), nn.Linear(HIDDEN, n_classes))

    def forward(self, x):
        return self.net(x)


def _train(model, f, seed):
    g = torch.Generator().manual_seed(seed)
    Xtr = torch.tensor(f.Xtr, dtype=torch.float32)
    ytr = torch.tensor(f.ytr)
    nc = f.n_classes
    cw = torch.tensor(len(f.ytr) / (nc * np.bincount(f.ytr, minlength=nc).clip(1)), dtype=torch.float32)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=WD)
    lossf = nn.CrossEntropyLoss(weight=cw)
    model.train()
    for _ in range(EPOCHS):
        perm = torch.randperm(len(ytr), generator=g)
        for i in range(0, len(perm), BATCH):
            b = perm[i:i + BATCH]
            opt.zero_grad()
            lossf(model(Xtr[b]), ytr[b]).backward()
            opt.step()
    model.eval()
    with torch.no_grad():
        return torch.softmax(model(torch.tensor(f.Xva, dtype=torch.float32)), 1).numpy()


def _gcn_fn(graph_fn, randomiser, f, seed):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    n = f.Xtr.shape[1]
    edges = graph_fn(f, rng)
    if randomiser is not None:
        edges = randomiser(edges, n, rng)
    f._n_edges = len(edges)
    p = -np.log10(np.clip(f.pvals, 1e-300, 1))
    prior = torch.tensor((p - p.min()) / (p.max() - p.min() + 1e-12), dtype=torch.float32)
    return _train(GCNPool(norm_adj(edges, n), prior, f.n_classes), f, seed)


def make_gcn(graph_fn, randomiser=None):
    """Picklable model function (module-level partial) for multiprocessing."""
    return functools.partial(_gcn_fn, graph_fn, randomiser)


def m_mlp(f, seed):
    torch.manual_seed(seed)
    nc = f.n_classes
    return _train(MLP(f.Xtr.shape[1], nc), f, seed)


def edge_count(f):
    return dict(n_edges=int(getattr(f, "_n_edges", -1)))
