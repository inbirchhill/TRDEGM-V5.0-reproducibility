"""
항목1: 결합운동방정식(2차 쿠라모토/'swing equation') 정의 + Lyapunov 지수로 카오스 판정

배경: 기존 v(전파속도) 자유변수 스캔 접근은 '목표값 맞춰 사후조정' 함정에 빠질 위험이
확인되어 확정 못한 채 보류됨(0716 세션2~7). 이번에는 그 접근을 완전히 버리고, 전력망
동기화 연구/조셉슨접합배열 연구에서 표준적으로 쓰이는 2차(관성형) 쿠라모토 모델(swing
equation)을 채택한다. 이 프레임의 장점:
  1. 라그랑주역학(+레일리 소산함수)에서 표준적으로 유도되는 형태 (전력망 스윙방정식 문헌 확인).
  2. 각 노드의 '자유회전 진동수'가 기존에 이미 검증된 값(§3.3의 nu_p=1, nu_a=1/sqrt2)으로
     정확히 고정되므로 새 자유 진동수 파라미터 도입이 필요 없음.
  3. 리아푸노프 지수(Benettin 알고리즘)로 '진짜 카오스' 여부를 정량적으로 판정 가능
     (v스캔의 R^2 피팅과 달리, 목표값에 맞추는 절차가 아니라 동역학계의 객관적 성질을 측정).

운동방정식 (허브 i, 위상 theta_i):
  theta_i'' = omega_i - gamma*theta_i' - K * sum_{j~i} sin(theta_i - theta_j)
  (라그랑주 L = sum 0.5*I*theta_i'^2 + sum omega_i*theta_i - sum_edges K*(1-cos(theta_i-theta_j)),
   레일리 소산 R = sum 0.5*gamma*theta_i'^2, I=1(단위관성, 별도 도입된 스케일 없음))

자유회전(결합 K=0) 정상상태 각속도 = omega_i/gamma 이므로, omega_i = gamma * Omega_i
로 설정하면 결합이 없을 때 각 허브가 정확히 기존 검증된 목표 각속도로 회전한다:
  Omega_i = 2*pi*nu_p = 2*pi*1        (N=6 허브, 3중점형/짧은채널)
  Omega_i = 2*pi*nu_a = 2*pi/sqrt(2)  (N=12 허브, 4중점형/긴채널)

남는 자유 파라미터는 오직 결합강도 K 하나 뿐(gamma=1 고정, omega_i=Omega_i가 되도록 스케일).
"""
import numpy as np
import networkx as nx
from network_base import (build_growing_shell, build_multicell_graph,
                            classify_hub_types, build_hub_graph)

def build_system(depth=0, buffer_radius=6):
    all_pts, vor, cluster_idx = build_growing_shell(buffer_radius=buffer_radius, depth=depth)
    G, face_seen, vertex_coord = build_multicell_graph(all_pts, vor, cluster_idx)
    hub_pos, hub_N, n_mixed = classify_hub_types(G)
    H = build_hub_graph(G, hub_pos, hub_N)
    nodes = list(H.nodes())
    idx = {n:i for i,n in enumerate(nodes)}
    n_hubs = len(nodes)
    Omega = np.array([2*np.pi*(1.0 if H.nodes[n]['N']==6 else 1.0/np.sqrt(2)) for n in nodes])
    # 인접 리스트 (배열 형태로, numba 없이도 빠르게)
    edges = np.array([[idx[a], idx[b]] for a,b in H.edges()])
    return nodes, idx, n_hubs, Omega, edges, H, cluster_idx

def rhs(state, n_hubs, Omega, edges, gamma, K):
    theta = state[:n_hubs]
    v = state[n_hubs:]
    dtheta = v
    coupling = np.zeros(n_hubs)
    dth = theta[edges[:,0]] - theta[edges[:,1]]
    s = np.sin(dth)
    np.add.at(coupling, edges[:,0], -K*s)
    np.add.at(coupling, edges[:,1],  K*s)
    dv = gamma*Omega - gamma*v + coupling
    return np.concatenate([dtheta, dv])

def rk4_step(state, dt, n_hubs, Omega, edges, gamma, K):
    k1 = rhs(state, n_hubs, Omega, edges, gamma, K)
    k2 = rhs(state + 0.5*dt*k1, n_hubs, Omega, edges, gamma, K)
    k3 = rhs(state + 0.5*dt*k2, n_hubs, Omega, edges, gamma, K)
    k4 = rhs(state + dt*k3, n_hubs, Omega, edges, gamma, K)
    return state + (dt/6.0)*(k1+2*k2+2*k3+k4)

def lyapunov_benettin(n_hubs, Omega, edges, gamma, K, dt=0.01, n_steps=6000,
                        renorm_every=10, d0=1e-8, transient_steps=2000, seed=0):
    rng = np.random.default_rng(seed)
    state = np.concatenate([rng.uniform(0,2*np.pi,n_hubs), np.zeros(n_hubs)])
    # 과도상태 제거
    for _ in range(transient_steps):
        state = rk4_step(state, dt, n_hubs, Omega, edges, gamma, K)
    # 참조궤적 + 섭동궤적
    pert = rng.normal(size=2*n_hubs)
    pert = pert/np.linalg.norm(pert)*d0
    state2 = state + pert
    log_sum = 0.0
    n_renorm = 0
    for step in range(n_steps):
        state = rk4_step(state, dt, n_hubs, Omega, edges, gamma, K)
        state2 = rk4_step(state2, dt, n_hubs, Omega, edges, gamma, K)
        if (step+1) % renorm_every == 0:
            diff = state2 - state
            dist = np.linalg.norm(diff)
            if dist == 0 or not np.isfinite(dist):
                break
            log_sum += np.log(dist/d0)
            n_renorm += 1
            state2 = state + diff*(d0/dist)
    total_time = n_renorm*renorm_every*dt
    lle = log_sum/total_time if total_time>0 else np.nan
    return lle

if __name__ == "__main__":
    import time
    t0=time.time()
    nodes, idx, n_hubs, Omega, edges, H, cluster_idx = build_system(depth=0, buffer_radius=6)
    print(f"depth=0: 셀 {len(cluster_idx)}, 허브 {n_hubs}, 엣지 {len(edges)}  ({time.time()-t0:.1f}s)")
    print(f"Omega 분포: N6={2*np.pi:.4f}, N12={2*np.pi/np.sqrt(2):.4f}")

    print("\n=== K 스캔 (gamma=1 고정) : 최대 리아푸노프 지수(LLE) ===")
    print(f"{'K':>8} {'LLE':>12} {'판정':>10}")
    for K in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]:
        t1=time.time()
        lle = lyapunov_benettin(n_hubs, Omega, edges, gamma=1.0, K=K,
                                  dt=0.01, n_steps=4000, transient_steps=1500)
        verdict = "혼돈(양수)" if lle>0.01 else ("경계" if lle>-0.01 else "안정(음수)")
        print(f"{K:8.2f} {lle:12.5f} {verdict:>10}  ({time.time()-t1:.1f}s)")
