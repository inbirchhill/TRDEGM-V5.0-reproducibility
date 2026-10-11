"""Track A/B 헤시안 안정성의 원인 점검 (V5.0 최종본2, §7.5 정정 / §13.3 / 부록 C.49)

복구된 96허브 격자(supercell_builder, m=2, 640결합, 위치 자유도 288)에서
정렬위상(theta_i = theta_j)의 위치 헤시안을 계산해 다음을 재현한다.

 (1) E_FS 단독: p = 2q(p=8, q=4)에서 결합별 반경 강성이 정확히 상쇄되어
     해석적 랭크 86, 1차 영모드 202개(평행이동 3개 포함)
 (2) E_FS + E_rphi(Track A, K0=9.724): 최소 고유값 약 -4597.79, 음의 모드 285개
     E_FS + E_rphi(Track B, K0=9.724): 음의 모드 0개, 최소 양의 고유값 약 2.0076
 (3) 안정성을 가르는 것은 K'(0)의 부호가 아니라 K''(0)의 부호
     (K'과 K''의 부호 조합을 바꾼 시험 함수 6개로 확인)

정렬위상에서 E_rphi = -sum K(s_ij) cos(dtheta) 는 추가 2체 퍼텐셜 -sum K(s_ij)로
작용하므로, 결합당 반경 강성에 -K''(s)/r0^2 가 더해진다.

필요 파일: supercell_builder.py, basis_data.npy (같은 폴더)
실행: python3 trackAB_hessian_curvature_check.py   (수 분 소요)
"""
import numpy as np
from supercell_builder import build_supercell

K0 = 9.724
q, p, eps = 4.0, 8.0, 1.0
SSTAR = 0.15

ns = build_supercell(2)
N = ns['Nsuper']; ei, ej, en = ns['edges_i'], ns['edges_j'], ns['edges_n']
r0 = ns['r0_edges']; box = ns['box0']; x0 = ns['x0']; typ = ns['atoms_type']
EG = ns['energy_grad']
shift = en @ box
print(f"격자: 허브 {N}개, 결합 {len(ei)}개, 이상격자 |grad|max = {ns['Gmax']:.1e}")


def geom(x):
    pos = x.reshape(N, 3)
    d = pos[ei] - pos[ej] - shift
    r = np.linalg.norm(d, axis=1)
    return r, d / r[:, None]


def grad_total(x, K, dK):
    """정렬위상에서 E = E_FS - sum K(s) 의 위치 그래디언트 (K=None 이면 E_FS 단독)"""
    _, G = EG(x)
    G = G.reshape(N, 3)
    if K is not None:
        r, u = geom(x)
        s = r / r0 - 1
        f = (-dK(s) / r0)[:, None] * u
        np.add.at(G, ei, f)
        np.add.at(G, ej, -f)
    return G.ravel()


def hessian(gfun, x, h=1e-5):
    n = len(x)
    H = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        H[:, k] = (gfun(x + e) - gfun(x - e)) / (2 * h)
    return 0.5 * (H + H.T)


def summary(name, H, tol=1e-6):
    w = np.sort(np.linalg.eigvalsh(H))
    pos = w[w > tol]
    print(f"  {name:40s} 최소={w[0]:11.4f}  음의모드={(w < -tol).sum():3d}  "
          f"최소양={pos.min() if len(pos) else float('nan'):.4f}")


# ---- (1) E_FS 단독: 해석적 헤시안 (p=2q 이면 sum_i F''(rho_i) grad(rho_i) grad(rho_i)^T 만 남음)
print("\n(1) E_FS 단독 헤시안")
r, u = geom(x0)
dph = -2 * q / r0 * np.exp(-2 * q * (r / r0 - 1))
Gr = np.zeros((N, 3 * N))
for k, (a, b) in enumerate(zip(ei, ej)):
    g = dph[k] * u[k]
    Gr[a, 3*a:3*a+3] += g; Gr[a, 3*b:3*b+3] -= g
    Gr[b, 3*a:3*a+3] += g; Gr[b, 3*b:3*b+3] -= g
rho = np.where(typ == 6, 10.0, 20.0)
H_an = (Gr.T * (eps * 0.25 * rho ** -1.5)) @ Gr
rank = np.linalg.matrix_rank(H_an, 1e-8)
H_FS = hessian(lambda x: grad_total(x, None, None), x0)
print(f"  해석적 랭크 = {rank}, 영모드 = {3*N - rank}개 (평행이동 3개 포함)")
print(f"  수치 헤시안과 해석 헤시안 차이 max = {np.abs(H_FS - H_an).max():.1e}")

# ---- (2) Track A/B
KA = lambda s: K0 * np.exp(-2 * q * s)
dKA = lambda s: -2 * q * K0 * np.exp(-2 * q * s)
def KB(s):
    x = (1 + s) / (1 + SSTAR); return K0 * x * np.exp(1 - x)
def dKB(s):
    x = (1 + s) / (1 + SSTAR); return K0 * np.exp(1 - x) * (1 - x) / (1 + SSTAR)

print("\n(2) 정렬위상 헤시안 (부록 C.49 재현)")
summary("E_FS + E_rphi (Track A)", hessian(lambda x: grad_total(x, KA, dKA), x0))
summary("E_FS + E_rphi (Track B)", hessian(lambda x: grad_total(x, KB, dKB), x0))

# ---- (3) K' 부호 vs K'' 부호
print("\n(3) 시험 함수: K(s) = K0 (1 + a s + b s^2)  ->  K'(0)=a K0, K''(0)=2b K0")
cases = [("K'<0, K''>0 (Track A형)", -8, 32), ("K'<0, K''<0", -8, -32),
         ("K'>0, K''>0", 8, 32), ("K'>0, K''<0 (Track B형)", 8, -32),
         ("K'=0, K''<0", 0, -32), ("K'<0, K''=0", -8, 0)]
for name, a, b in cases:
    K = lambda s, a=a, b=b: K0 * (1 + a * s + b * s ** 2)
    dK = lambda s, a=a, b=b: K0 * (a + 2 * b * s)
    summary(name, hessian(lambda x: grad_total(x, K, dK), x0))
print("\n=> 음의 모드는 K''>0 일 때만 나타난다(K'의 부호와 무관).")
