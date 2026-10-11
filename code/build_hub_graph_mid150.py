"""hub_graph_mid150.npz 재생성 / Regenerate hub_graph_mid150.npz

[KO] 원본 hub_graph_mid150.npz(174허브)가 유실되어, network_base.py로 다시 만든다.
     구성: 성장껍질 depth=1, buffer_radius=5 (허브 300개) -> 원점으로부터 거리제곱
     r^2 <= 8.75 인 허브만 남긴 유도부분그래프(174허브, networkx 노드 순서 유지).
     고유진동수 omega = 2*pi (N6), 2*pi/sqrt(2) (N12)  (swing_eom.py와 같은 규약).
     검증: (1) mid150_matched_worker.py의 소스 선택(rng seed 7)이 원본 큐의 소스
     12개와 정확히 일치, (2) 무킥 시행 하나를 다시 돌린 결과가 원본 체크포인트와
     차이 0, (3) analyze_hop_result.py가 p=0.00806(홉1 0.1907, 홉3 0.1586, 10/12)을
     그대로 재현.

[EN] The original hub_graph_mid150.npz (174 hubs) was lost; this script rebuilds it
     with network_base.py. Construction: growing shell depth=1, buffer_radius=5
     (300 hubs) -> induced subgraph on hubs with squared distance r^2 <= 8.75 from
     the origin (174 hubs, networkx node order preserved). Natural frequencies
     omega = 2*pi (N6), 2*pi/sqrt(2) (N12), the same convention as swing_eom.py.
     Checks: (1) the source selection in mid150_matched_worker.py (rng seed 7)
     reproduces the 12 sources in the original queue exactly; (2) rerunning one
     no-kick trial matches the original checkpoint with zero difference;
     (3) analyze_hop_result.py reproduces p=0.00806 (hop1 0.1907, hop3 0.1586, 10/12).

실행 / Run:  python3 build_hub_graph_mid150.py   (network_base.py 필요 / requires network_base.py)
"""
import json
import numpy as np
import networkx as nx
from network_base import (build_growing_shell, build_multicell_graph,
                          classify_hub_types, build_hub_graph)

all_pts, vor, cluster_idx = build_growing_shell(buffer_radius=5, depth=1)
G, _, _ = build_multicell_graph(all_pts, vor, cluster_idx)
hub_pos, hub_N, _ = classify_hub_types(G)
H = build_hub_graph(G, hub_pos, hub_N)

nodes = [n for n in H.nodes() if np.linalg.norm(hub_pos[n]) ** 2 <= 8.75 + 1e-6]
A = nx.to_numpy_array(H.subgraph(nodes), nodelist=nodes)
types = np.array([hub_N[n] for n in nodes])
omega = 2 * np.pi * np.where(types == 6, 1.0, 1.0 / np.sqrt(2))
np.savez('hub_graph_mid150.npz', A=A, omega=omega)
print(f"저장/saved: hub_graph_mid150.npz  (hubs={len(nodes)}, edges={int(A.sum() // 2)}, "
      f"N6={int((types == 6).sum())}, N12={int((types == 12).sum())})")

# --- 소스 선택 지문 확인 / source-selection fingerprint check
deg = A.sum(axis=1).astype(int)
reg = sorted(i for i in range(len(deg)) if deg[i] == 20)
rng = np.random.default_rng(7)
picked = [int(x) for x in rng.choice(reg, size=min(6, len(reg)), replace=False)]
irr = [i for i in range(len(deg)) if deg[i] != 20]
picked += [int(x) for x in rng.choice(irr, size=6, replace=False)]
try:
    expected = json.load(open('mid150_matched_queue.json'))['sources']
    print("원본 큐 소스와 일치 / matches original queue sources:", picked == expected)
except FileNotFoundError:
    print("mid150_matched_queue.json 없음 / not found; picked sources:", picked)
