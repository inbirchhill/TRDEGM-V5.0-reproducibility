"""[분석 스크립트] mid150_matched_worker.py(또는 그 220/260대 버전)가
만든 원자료(*_ckpt.npz + *_queue.json)에서 "홉1>홉3" 판정과 p값을
계산하는 정확한 로직. §3.296/302/331에서 반복 사용된 것과 동일.

사용법: 파일 상단의 CKPT/QUEUE/GRAPH 세 경로만 해당 규모의 파일로
바꿔서 실행하면 됨(174/220/260 등 공통으로 재사용 가능).
"""
import numpy as np
import networkx as nx
from scipy import stats
import json

# ---- 이 세 줄만 대상 규모에 맞게 교체 ----
GRAPH_FILE = 'hub_graph_mid150.npz'      # 예: hub_graph_mid220.npz
CKPT_FILE = 'mid150_matched_ckpt.npz'    # 예: mid220_matched_ckpt.npz
QUEUE_FILE = 'mid150_matched_queue.json' # 예: mid220_matched_queue.json
# -----------------------------------------

d1 = np.load(GRAPH_FILE)
A = d1['A']
n = len(A)
G = nx.from_numpy_array(A)

ckpt = np.load(CKPT_FILE)
with open(QUEUE_FILE) as f:
    qdata = json.load(f)
sources = qdata['sources']

all_h1, all_h3 = [], []
per_source_detail = []

for source in sources:
    # 1) 이 소스로부터 각 노드까지의 홉거리 계산
    hops = dict(nx.single_source_shortest_path_length(G, source, cutoff=10))
    hop1_nodes = [v for v, h in hops.items() if h == 1]
    hop3_nodes = [v for v, h in hops.items() if h == 3]
    if not hop1_nodes or not hop3_nodes:
        print(f"source={source}: 홉1 또는 홉3 노드 없음, 스킵")
        continue

    # 2) 원자료 불러오기: shape=(시행수, 노드수), 각 원소는 그 시행에서
    #    관찰창 동안의 |상대속도| 최댓값(run_condition의 출력)
    rk = ckpt[f's{source}_kick']    # 킥 조건, 예: (6, n)
    rn = ckpt[f's{source}_nokick']  # 무킥 조건(자연변동 기준선), 예: (6, n)

    # 3) 핵심: "성공 판정 임계값"은 노드별로 다르게 잡는다.
    #    무킥 조건에서 그 노드가 자연적으로 보이는 반응의 95백분위수를
    #    "우연으로도 나올 수 있는 최대치"로 보고, 킥 조건의 반응이
    #    이를 넘으면 "진짜 반응"으로 판정한다(노드마다 자연변동 크기가
    #    다르므로 전역 단일 임계값이 아니라 노드별 임계값을 씀).
    thresh = np.percentile(rn, 95, axis=0)  # shape=(n,), 노드별 임계값

    # 4) 각 노드 t에서 "성공확률" = 킥조건 시행들 중 임계값을 넘은 비율
    #    (예: 6회 시행 중 4회가 임계값 초과 -> 그 노드의 성공확률 4/6)
    def success_rate(node_idx):
        return (rk[:, node_idx] > thresh[node_idx]).mean()

    s1 = [success_rate(t) for t in hop1_nodes]  # 홉1 노드들 각각의 성공확률
    s3 = [success_rate(t) for t in hop3_nodes]  # 홉3 노드들 각각의 성공확률

    # 5) 이 소스에서의 대표값 = 홉거리별 성공확률의 평균
    h1_mean = np.mean(s1)
    h3_mean = np.mean(s3)
    all_h1.append(h1_mean)
    all_h3.append(h3_mean)
    per_source_detail.append((source, h1_mean, h3_mean, h1_mean > h3_mean))

# 6) 소스 12개(또는 해당 개수)를 짝지어 Wilcoxon 부호순위검정
#    (귀무가설: 홉1과 홉3 사이에 차이 없음, 대립가설: 홉1 > 홉3, 한쪽검정)
all_h1 = np.array(all_h1)
all_h3 = np.array(all_h3)
n_pos = (all_h1 > all_h3).sum()
n_total = len(all_h1)
stat, p = stats.wilcoxon(all_h1, all_h3, alternative='greater')

print("=== 소스별 상세 ===")
for source, h1, h3, pos in per_source_detail:
    print(f"  source={source}: 홉1={h1:.4f} 홉3={h3:.4f} {'(양성)' if pos else '(음성)'}")

print(f"\n=== 종합 ===")
print(f"홉1 평균={all_h1.mean():.4f}, 홉3 평균={all_h3.mean():.4f}")
print(f"홉1>홉3: {n_pos}/{n_total}")
print(f"Wilcoxon 부호순위검정(한쪽, greater): p={p:.5f}")
