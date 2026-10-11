"""[그리드 스케일링 문제, 우선순위1(c)] depth=0(94)과 depth=1(300)
사이 중간규모(174허브)에서 §3.296/302와 매칭되는 조건(K=0.09,
T=15, 창=30, 12소스, 시행6회)으로 홉거리 국소성 재검정.
"""
import numpy as np
import networkx as nx
import os
import json
import time

d1 = np.load('hub_graph_mid150.npz')
A = d1['A']
omega = d1['omega']
n = len(omega)
G = nx.from_numpy_array(A)
deg = dict(G.degree())

K, gamma, noise_amp = 0.09, 0.5, 1.0
dt = 0.01
kick_magnitude = 15.0
response_window = 30
n_trials = 6
T_total = 15.0

CKPT = 'mid150_matched_ckpt.npz'
QUEUE_FILE = 'mid150_matched_queue.json'

# 규칙패턴(고배위수20) 6개 + 불규칙패턴(무작위) 6개
regular_candidates = sorted([i for i in range(n) if deg[i] == 20])
rng_pick = np.random.default_rng(7)
regular_sources = list(rng_pick.choice(regular_candidates,
                                         size=min(6, len(regular_candidates)), replace=False))
irregular_candidates = [i for i in range(n) if deg[i] != 20]
irregular_sources = list(rng_pick.choice(irregular_candidates, size=6, replace=False))
all_sources = [int(s) for s in regular_sources + irregular_sources]

if not os.path.exists(QUEUE_FILE):
    queue = []
    for s in all_sources:
        queue.append({'source': s, 'condition': 'kick', 'done': False})
        queue.append({'source': s, 'condition': 'nokick', 'done': False})
    with open(QUEUE_FILE, 'w') as f:
        json.dump({'queue': queue, 'sources': all_sources}, f)
    print(f"작업큐 신규생성: 소스{all_sources}, 총{len(queue)}개 작업단위")
else:
    with open(QUEUE_FILE) as f:
        qdata = json.load(f)
    queue = qdata['queue']
    all_sources = qdata['sources']

if os.path.exists(CKPT):
    ckpt = dict(np.load(CKPT, allow_pickle=True))
else:
    ckpt = {}


def derivs(theta, phi):
    diff = theta[None, :] - theta[:, None]
    coupling = (A * np.sin(diff)).sum(axis=1)
    return phi, -gamma * phi + omega + K * coupling


def run_condition(source, do_kick, base_seed):
    n_steps = int(T_total / dt)
    kick_step = int(n_steps * 0.5)
    results = []
    for trial in range(n_trials):
        rng = np.random.default_rng(base_seed + trial)
        theta = rng.uniform(-0.1, 0.1, n)
        phi = omega.copy()
        recording = False
        rec_count = 0
        phi_track = []
        for step in range(n_steps):
            k1t, k1p = derivs(theta, phi)
            k2t, k2p = derivs(theta + 0.5 * dt * k1t, phi + 0.5 * dt * k1p)
            k3t, k3p = derivs(theta + 0.5 * dt * k2t, phi + 0.5 * dt * k2p)
            k4t, k4p = derivs(theta + dt * k3t, phi + dt * k3p)
            theta = theta + (dt / 6) * (k1t + 2 * k2t + 2 * k3t + k4t)
            phi = phi + (dt / 6) * (k1p + 2 * k2p + 2 * k3p + k4p) + \
                noise_amp * rng.normal(0, 1, n) * np.sqrt(dt)
            if do_kick and step == kick_step:
                phi[source] += kick_magnitude
                recording = True
            if not do_kick and step == kick_step:
                recording = True
            if recording:
                phi_track.append(phi.copy())
                rec_count += 1
                if rec_count >= response_window:
                    recording = False
        phi_track = np.array(phi_track)
        rel = phi_track - phi_track.mean(axis=1, keepdims=True)
        results.append(np.abs(rel).max(axis=0))
    return np.array(results)


t_start = time.time()
n_done_this_call = 0
for item in queue:
    if item['done']:
        continue
    if time.time() - t_start > 140:
        break
    source = item['source']
    cond = item['condition']
    key = f"s{source}_{cond}"
    base_seed = 9000 + source * 31 + (0 if cond == 'kick' else 500)
    result = run_condition(source, cond == 'kick', base_seed)
    ckpt[key] = result
    item['done'] = True
    n_done_this_call += 1
    print(f"완료: {key} (소요 {time.time()-t_start:.1f}초 누적)")

np.savez(CKPT, **ckpt)
with open(QUEUE_FILE, 'w') as f:
    json.dump({'queue': queue, 'sources': all_sources}, f)

n_remaining = sum(1 for item in queue if not item['done'])
print(f"\n이번 호출 완료: {n_done_this_call}개, 남은 작업: {n_remaining}개")
if n_remaining == 0:
    print(">>> 전체 작업큐 완료!")
else:
    print(">>> 이어서 진행")
