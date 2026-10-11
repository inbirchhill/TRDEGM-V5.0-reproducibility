"""
§3.133/134/135 번호충돌 규명: depth=0 (공식판 §3.134)과 depth=1 (초안 §3.135)을
'같은 절차'로 나란히 재실행 — K1 grid, settle_T, seed 전부 통일.

절차(통일):
- K1 grid: 0.05~0.55, 0.15~0.35 구간은 0.01 간격으로 촘촘히
- settle_T = 100 (dt=0.01, 10000 steps), transient 전부 포함(연속 적분, 별도 분리 없음 — 원 실험과 동일)
- v_std: 마지막 window(마지막 2000 steps)에서 각 노드 각속도의 "시간에 대한 표준편차"의 노드 평균
  (완전잠금이면 각속도가 시간에 대해 상수이므로 std->0)
- 질서변수 r: 마지막 시점의 Kuramoto order parameter |mean_i exp(i*theta_i)|
  (완전잠금 상태에서는 상대위상이 고정되므로 시간에 무관한 값)
- 시드의존성 대조: 문턱 아래/위 후보점에서 seed 0~5 반복
"""
import numpy as np, time, json
from swing_eom import build_system, rk4_step

def run_one(n_hubs, Omega, edges, gamma, K, dt=0.01, n_steps=10000, seed=0, window=2000):
    rng = np.random.default_rng(seed)
    state = np.concatenate([rng.uniform(0, 2*np.pi, n_hubs), np.zeros(n_hubs)])
    v_hist = np.zeros((window, n_hubs))
    for step in range(n_steps):
        state = rk4_step(state, dt, n_hubs, Omega, edges, gamma, K)
        if step >= n_steps - window:
            v_hist[step - (n_steps - window)] = state[n_hubs:]
    theta_final = state[:n_hubs]
    r = np.abs(np.mean(np.exp(1j*theta_final)))
    v_std = np.std(v_hist, axis=0).mean()
    return r, v_std

def scan_depth(depth, K_grid, seed=0):
    nodes, idx, n_hubs, Omega, edges, H, cluster_idx = build_system(depth=depth, buffer_radius=6)
    results = []
    for K in K_grid:
        t0=time.time()
        r, v_std = run_one(n_hubs, Omega, edges, gamma=1.0, K=K, seed=seed)
        results.append(dict(K=float(K), r=float(r), v_std=float(v_std)))
        print(f"  depth={depth} K={K:.3f} r={r:.6f} v_std={v_std:.3e}  ({time.time()-t0:.1f}s)")
    return n_hubs, results

if __name__ == "__main__":
    coarse = list(np.round(np.arange(0.05, 0.15, 0.02), 3))
    fine   = list(np.round(np.arange(0.15, 0.36, 0.01), 3))
    coarse2= list(np.round(np.arange(0.36, 0.56, 0.05), 3))
    K_grid = sorted(set(coarse+fine+coarse2))

    out = {}
    for depth in [0, 1]:
        print(f"=== depth={depth} scan (seed=0) ===")
        n_hubs, results = scan_depth(depth, K_grid, seed=0)
        out[f"depth{depth}"] = dict(n_hubs=n_hubs, scan=results)

    with open("k1_scan_results.json", "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("저장 완료: k1_scan_results.json")
