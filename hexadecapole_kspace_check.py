# Copyright (c) 2026 Hyeongsu Kim (김형수)
# Licensed under CC BY-NC-SA 4.0 (저작자표시-비영리-동일조건변경허락 4.0)
# https://creativecommons.org/licenses/by-nc-sa/4.0/deed.ko
# Part of: "Rhombic Dodecahedral Energy Grid Model V4.0"

"""Tier0 항목 0-2: §6.3 다중극-이탈 재시도. 기존 방법(실공간에서 같은 |R|의
방향별 그린함수를 비교)은 정류위상 진동에 묻혀 결론이 안 났다(50-51부).
대신 k공간에서 직접: 그래프 라플라시안 L(k)의 0으로 가는(acoustic) 분지
lambda(k)를 |k|->0 근방에서 방향별로 전개해, 등방항(k^2)과 4차 이방항
(k^4, Oh군의 육극자류 큐빅조화함수)의 계수비를 직접 뽑아낸다."""
import numpy as np

basis = np.load('basis_data.npy', allow_pickle=True).item()
basis_items = list(basis.items())
n_basis = len(basis_items)
L_const = 2.0
canon_pos = {fk: np.array(fk) for fk, _ in basis_items}
fk_list = [fk for fk, _ in basis_items]
fk_idx = {fk: i for i, fk in enumerate(fk_list)}

def frac_key(p, nd=4):
    return tuple(np.round(np.mod(p, L_const), nd))

hoppings = []
for a_idx, (fk_a, b) in enumerate(basis_items):
    pa = canon_pos[fk_a]
    for rel, nbr_type in b['template']:
        p_nbr = pa + rel
        fk_b = frac_key(p_nbr)
        b_idx = fk_idx[fk_b]
        pb_canon = canon_pos[fk_b]
        R_int = np.round((p_nbr - pb_canon) / L_const).astype(int)
        hoppings.append((a_idx, b_idx, R_int, np.linalg.norm(rel)))

print(f"기저 {n_basis}개, 홉핑 {len(hoppings)}개 (평균 배위수 {len(hoppings)/n_basis:.1f})")

hop_a = np.array([h[0] for h in hoppings])
hop_b = np.array([h[1] for h in hoppings])
hop_R = np.array([h[2] for h in hoppings]) * L_const  # 실공간 벡터

def L_of_k(kvec):
    """그래프 라플라시안 L(k) = D - A(k). 가중치는 전부 1(단순그래프)."""
    Lk = np.zeros((n_basis, n_basis), dtype=complex)
    phase = np.exp(1j * (hop_R @ kvec))
    for idx in range(len(hop_a)):
        a, b = hop_a[idx], hop_b[idx]
        Lk[a, a] += 1.0          # degree
        Lk[a, b] -= phase[idx]   # -adjacency(phase)
    return Lk

def acoustic_eigenvalue(kvec):
    """0으로 수렴하는 분지(가장 작은 고유값)를 반환."""
    if np.allclose(kvec, 0):
        return 0.0
    Lk = L_of_k(kvec)
    evals = np.linalg.eigvalsh(Lk)
    return evals[0].real  # 가장 작은 고유값 = acoustic 분지

# 여러 방향(구면 균일 샘플)에서 작은 |k|의 고유값을 측정
def fibonacci_sphere(n):
    i = np.arange(0, n, dtype=float) + 0.5
    phi = np.arccos(1 - 2*i/n)
    golden = np.pi*(1+5**0.5)
    theta = golden*i
    x = np.cos(theta)*np.sin(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(phi)
    return np.stack([x,y,z], axis=1)

dirs = fibonacci_sphere(400)

# 두 개의 작은 |k|에서 측정해 유한차분으로 4차항(anisotropic quartic)을 분리
# lambda(k) ~ c2*k^2 + c4_iso*k^4 + c4_aniso*k^4*K4(khat), K4 = (kx^4+ky^4+kz^4)/k^4 - 3/5
results = []
for kmag in [0.02, 0.04]:
    lam = np.array([acoustic_eigenvalue(kmag*d) for d in dirs])
    results.append(lam)
    print(f"|k|={kmag}: lambda 평균={lam.mean():.6e}  표준편차={lam.std():.3e}  변동계수={lam.std()/lam.mean()*100:.4f}%")

lam1, lam2 = results
k1, k2 = 0.02, 0.04
# lambda(k) = c2*k^2 + c4(khat)*k^4  (c4는 방향의존)
# 두 k에서 연립: c2 = (lam1/k1^2*k2^4 - lam2/k2^2*k1^4)/(k2^4-k1^4)*... 더 간단히 각 방향별로 직접 풀기
c4_dir = (lam2/k2**2 - lam1/k1**2) / (k2**2 - k1**2)   # per-direction quartic coefficient
c2_dir = lam1/k1**2 - c4_dir*k1**2

print(f"\nc2(등방 2차계수): 평균={c2_dir.mean():.6f}  변동계수={c2_dir.std()/c2_dir.mean()*100:.4f}%")
print(f"c4(방향의존 4차계수): 평균={c4_dir.mean():.6f}  표준편차={c4_dir.std():.6f}  변동계수={c4_dir.std()/abs(c4_dir.mean())*100:.2f}%")

# c4의 방향의존성을 Oh군 큐빅조화함수(K4)로 회귀분해
K4 = (dirs[:,0]**4 + dirs[:,1]**4 + dirs[:,2]**4) - 0.6  # (kx^4+ky^4+kz^4)/|k|^4 - 3/5, |k_hat|=1
A = np.vstack([np.ones_like(K4), K4]).T
coef, res, rank, sv = np.linalg.lstsq(A, c4_dir, rcond=None)
c4_iso, c4_aniso = coef
pred = A @ coef
ss_res = np.sum((c4_dir-pred)**2)
ss_tot = np.sum((c4_dir-c4_dir.mean())**2)
r2 = 1 - ss_res/ss_tot
print(f"\nc4 = {c4_iso:.6f} + {c4_aniso:.6f} * K4(k_hat)   (R^2={r2:.5f})")
print(f"이방성 비율(육극자/등방, |c4_aniso/c4_iso|) = {abs(c4_aniso/c4_iso)*100:.3f}%")

# 축방향 vs 체대각선 방향 직접 비교(참고용)
axis_dir = np.array([0,0,1.0])
body_dir = np.array([1,1,1.0])/np.sqrt(3)
for kmag in [0.02,0.04]:
    la = acoustic_eigenvalue(kmag*axis_dir)
    lb = acoustic_eigenvalue(kmag*body_dir)
    print(f"|k|={kmag}: 축방향 lambda={la:.6e}  체대각 lambda={lb:.6e}  비율={la/lb:.5f}")
