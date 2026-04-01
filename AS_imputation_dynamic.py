import os
import argparse
import numpy as np
import pandas as pd
from tqdm import trange, tqdm
from concurrent.futures import ThreadPoolExecutor

# ======================
#source activate scasl
#python AS_imputation_dynamic.py  workdir --k 3 --threads 100
# 参数解析
# ======================
parser = argparse.ArgumentParser(
    description="Dynamic KNN imputation for AS matrix (multithread)"
)
parser.add_argument(
    "workdir",
    help="Directory containing AS.csv and gene_matrix.csv"
)
parser.add_argument(
    "--k",
    type=int,
    default=3,
    help="Number of nearest neighbors (default: 3)"
)
parser.add_argument(
    "--threads",
    type=int,
    default=100,
    help="Number of threads (default: 100)"
)

args = parser.parse_args()
workdir = args.workdir
k = args.k
max_workers = args.threads

# ======================
# 路径定义
# ======================
as_file = os.path.join(workdir, "AS_abs.csv")
dm_file = os.path.join(workdir, "Gene_matrix.csv")
out_file = os.path.join(
    workdir,
    f"AS_imputed_Dynamic_k{k}_multithread.csv"
)

# ======================
# 读入数据
# ======================
df = pd.read_csv(as_file, index_col=0)
dm = pd.read_csv(dm_file, index_col=0)

# ======================
# 函数定义
# ======================
def fill_na(df):
    return df.fillna(0.5)

def knn_impute_dynamic_multithread(df_fillna, dm, df, k=3, max_workers=100):
    array_dm = dm.values.T   # cells × features
    m = array_dm.shape[0]
    n = df.shape[0]

    # 1. 距离矩阵
    d = np.zeros((m, m))
    for i in trange(m - 1, desc="Computing distance matrix"):
        for j in range(i, m):
            d[i, j] = d[j, i] = np.sqrt(
                np.sum((array_dm[i] - array_dm[j]) ** 2)
            )

    # 2. 邻居排序
    all_neighbors = np.argsort(d, axis=1)
    dist_std = d.std() + 1e-8

    df_values = df.values
    df_fillna_values = df_fillna.values

    def impute_for_cell(i):
        col_res = df_values[:, i].copy()

        for row_idx in range(n):
            if pd.isna(df_values[row_idx, i]):
                neighbor_ids = []
                neighbor_wts = []

                for neigh in all_neighbors[i]:
                    if neigh == i:
                        continue
                    val = df_values[row_idx, neigh]
                    if not pd.isna(val):
                        dist = d[i, neigh]
                        neighbor_ids.append(neigh)
                        neighbor_wts.append(np.exp(-2 * dist / dist_std))
                        if len(neighbor_ids) >= k:
                            break

                if neighbor_ids:
                    neighbor_vals = df_values[row_idx, neighbor_ids]
                    neighbor_wts = np.array(neighbor_wts)
                    neighbor_wts /= neighbor_wts.sum()
                    col_res[row_idx] = np.average(
                        neighbor_vals, weights=neighbor_wts
                    )
                else:
                    col_res[row_idx] = df_fillna_values[row_idx, i]

        return col_res

    # 3. 多线程
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(
            tqdm(
                executor.map(impute_for_cell, range(m)),
                total=m,
                desc="KNN multithread imputation"
            )
        )

    df_knn = pd.DataFrame(
        np.column_stack(results),
        index=df.index,
        columns=df.columns
    )

    return df_knn

# ======================
# 运行
# ======================
df_fillna = fill_na(df)
df_knn_imputed = knn_impute_dynamic_multithread(
    df_fillna,
    dm,
    df,
    k=k,
    max_workers=max_workers
)

df_knn_imputed.to_csv(out_file)

print(f"Imputation finished. Output saved to:\n{out_file}")
