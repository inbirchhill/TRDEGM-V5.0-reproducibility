"""구조적 전파(도착시간 vs 홉거리) 재구성 검증.
mid150_matched_worker.py와 같은 인프라(K=0.09, gamma=0.5, dt=0.01,
hub_graph_mid150.npz, kick_magnitude=15.0)를 그대로 사용하되,
"도착시간"(phi_rel이 임계값을 처음 넘는 시각)을 명시적으로 기록한다.
"""
import numpy as np
import networkx as nx

d1 = np.load('hub_graph_mid150.npz')
A = d1['A']; omega = d1['omega']; n = len(omega)
G = nx.from_numpy_array(A)

K, gamma, noise_amp = 0.09, 0.5, 1.0
dt = 0.01
kick_magnitude = 15.0
T_total = 15.0
n_steps = int(T_total/dt)
kick_step = int(n_steps*0.5)
track_steps = 3000  # 킥 이후 30 model-time units 추적

def derivs(theta, phi):
    diff = theta[None,:]-theta[:,None]
    coupling = (A*np.sin(diff)).sum(axis=1)
    return phi, -gamma*phi+omega+K*coupling

def run_trial(source, seed):
    rng = np.random.default_rng(seed)
    theta = rng.uniform(-0.1,0.1,n)
    phi = omega.copy()
    arrival_step = np.full(n, -1)
    hop_dist = nx.shortest_path_length(G, source=source)
    for step in range(n_steps):
        k1t,k1p = derivs(theta,phi)
        k2t,k2p = derivs(theta+0.5*dt*k1t, phi+0.5*dt*k1p)
        k3t,k3p = derivs(theta+0.5*dt*k2t, phi+0.5*dt*k2p)
        k4t,k4p = derivs(theta+dt*k3t, phi+dt*k3p)
        theta = theta + (dt/6)*(k1t+2*k2t+2*k3t+k4t)
        phi = phi + (dt/6)*(k1p+2*k2p+2*k3p+k4p) + noise_amp*rng.normal(0,1,n)*np.sqrt(dt)
        if step == kick_step:
            phi_baseline = phi.copy()
            phi[source] += kick_magnitude
        if step > kick_step and step-kick_step <= track_steps:
            rel = np.abs(phi - phi_baseline)
            newly_arrived = (arrival_step < 0) & (rel > 2.0)  # 임계값=2.0(노이즈 대비 충분히 큼)
            arrival_step[newly_arrived] = step - kick_step
        if step-kick_step > track_steps:
            break
    return arrival_step, hop_dist

deg = dict(G.degree())
regular_candidates = sorted([i for i in range(n) if deg[i]==20])
rng_pick = np.random.default_rng(7)
sources = list(rng_pick.choice(regular_candidates, size=3, replace=False))

all_hops, all_times = [], []
for s in sources:
    arrival_step, hop_dist = run_trial(int(s), seed=42+s)
    for node in range(n):
        if node==s: continue
        if arrival_step[node] > 0:
            h = hop_dist.get(node, None)
            if h is not None and h <= 8:
                all_hops.append(h)
                all_times.append(arrival_step[node]*dt)

all_hops = np.array(all_hops); all_times = np.array(all_times)
print(f"수집된 (홉거리,도착시간) 쌍: {len(all_hops)}개")
from scipy.stats import linregress, pearsonr
if len(all_hops) > 5:
    res = linregress(all_hops, all_times)
    r,p = pearsonr(all_hops, all_times)
    print(f"선형회귀: 기울기(=1/속도)={res.slope:.4f}, R^2={res.rvalue**2:.4f}")
    print(f"피어슨 r={r:.4f}, p={p:.4e}")
    if res.slope != 0:
        print(f"유효 전파속도 = {1/res.slope:.4f} (hop/model-time)")
