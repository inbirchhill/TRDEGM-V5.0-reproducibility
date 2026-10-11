"""
결합의 세기(K) 기하학적 유도 시도 1차: 다종 Finnis-Sinclair 포텐셜에서
결합타입별 '유효 결합 강성(local bond force constant)'을 해석적으로 계산.

방법론적 근거(외주 대체, 문헌 확인):
Riffe, Christensen & Wilson, "Vibrational Dynamics within the Embedded-Atom-Method
Formalism and the Relationship to Born-von-Karman Force Constants" (arXiv:1806.04244, 2018)
-- EAM류 포텐셜에서 결합별 BvK 힘상수를 유도하는 것이 표준적으로 확립된 절차임을 확인.
여기서는 그 절차의 최소 형태(다른 이웃은 고정한 "동결 환경" 근사, nearest-neighbor만)를 적용.

E_total = sum_i F(rho_i) + sum_<ij> V(r_ij)
F(rho) = -eps*sqrt(rho),  F'(rho) = -eps/(2 sqrt(rho)),  F''(rho) = eps/(4 rho^{3/2})
phi(r) = exp(-2q(r/r0-1)),  V(r) = A*exp(-p(r/r0-1))   (r0 = 결합타입의 이상거리)

한 결합(i,j)의 신장 r에 대해 (다른 이웃들의 rho 기여는 고정):
  rho_i(r) = rho_i^bulk - phi(r0) + phi(r) = rho_i^bulk - 1 + phi(r)   (phi(r0)=1)
  E_bond(r) = F(rho_i(r)) + F(rho_j(r)) + V(r)
  k_bond = d^2 E_bond / dr^2 |_{r=r0}
         = F''(rho_i)*phi'(r0)^2 + F'(rho_i)*phi''(r0)
         + F''(rho_j)*phi'(r0)^2 + F'(rho_j)*phi''(r0)
         + V''(r0)
"""
import numpy as np

q, p, eps = 4.0, 8.0, 1.0

def Fp(rho):  return -eps/(2*np.sqrt(rho))
def Fpp(rho): return eps/(4*rho**1.5)

# 결합타입별 r0, 그리고 양끝 허브의 벌크 배위수(rho^bulk = 배위수, phi(r0)=1이므로)
bonds = {
    "6-6":   dict(r0=1.0,       rho_i=10, rho_j=10),   # 짧은경로쪽-짧은경로쪽
    "6-12":  dict(r0=0.8660254, rho_i=10, rho_j=20),   # 짧은경로쪽-긴경로쪽
    "12-12": dict(r0=1.4142136, rho_i=20, rho_j=20),   # 긴경로쪽-긴경로쪽
}

# 결합타입별 캘리브레이션된 A_e (기존 §68 해석적 결과, q=4,p=8,eps=1 기준 재확인용으로 재계산)
def calibrate_A(r0, rho_i, rho_j):
    # A_e = -(2q/p) * ( F'(rho_i) + F'(rho_j) )   [완전배위 상태 결합력=0]
    return -(2*q/p) * (Fp(rho_i) + Fp(rho_j))

print("=== 캘리브레이션 재확인 ===")
for name, b in bonds.items():
    A = calibrate_A(b["r0"], b["rho_i"], b["rho_j"])
    print(f"  {name}: r0={b['r0']:.4f}  A_e={A:.4f}")

print("\n=== 결합타입별 유효 국소 힘상수 k_bond = d^2E/dr^2|_r0 (동결환경 근사) ===")
k_vals = {}
for name, b in bonds.items():
    r0, rho_i, rho_j = b["r0"], b["rho_i"], b["rho_j"]
    A = calibrate_A(r0, rho_i, rho_j)
    # phi, phi', phi'' at r0
    phi_p  = -2*q/r0                      # phi'(r0) = -2q/r0 * phi(r0) = -2q/r0 (phi(r0)=1)
    phi_pp = (2*q/r0)**2                  # phi''(r0) = (2q/r0)^2 * phi(r0) = (2q/r0)^2
    # V, V', V'' at r0
    Vpp = A*(p/r0)**2                      # V(r)=A*exp(-p(r/r0-1)); V''(r0)=A*(p/r0)^2
    k_embed = (Fpp(rho_i)+Fpp(rho_j))*phi_p**2 + (Fp(rho_i)+Fp(rho_j))*phi_pp
    k_pair  = Vpp
    k_total = k_embed + k_pair
    k_vals[name] = k_total
    print(f"  {name}: k_embed={k_embed:8.4f}  k_pair={k_pair:8.4f}  k_total={k_total:8.4f}")

kref = k_vals["6-6"]
print("\n=== 상대비 (6-6 기준=1) ===")
for name, k in k_vals.items():
    print(f"  {name}: {k/kref:.4f}")

# 허브종(N=6, N=12)별 '유효 국소강성' = 그 허브가 갖는 모든 결합의 k_bond 합
# (심부 배위: N=6허브는 6-6이웃6개+6-12이웃4개, N=12허브는 12-12이웃12개+6-12이웃8개)
k66, k612, k1212 = k_vals["6-6"], k_vals["6-12"], k_vals["12-12"]
K_N6  = 6*k66 + 4*k612
K_N12 = 12*k1212 + 8*k612
print(f"\n=== 허브종별 총 국소강성(모든 결합의 k 합) ===")
print(f"  N=6 허브:  {K_N6:.4f}")
print(f"  N=12 허브: {K_N12:.4f}   (N12/N6 비율 = {K_N12/K_N6:.4f})")
