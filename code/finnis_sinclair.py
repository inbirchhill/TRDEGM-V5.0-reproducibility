"""Finnis-Sinclair류 다체 포텐셜(임베딩 함수 F(rho)=-sqrt(rho))로 실제 표면이완 재현 시도.
핵심: 각 결합의 실질적 세기가 그 원자의 '국소배위밀도' rho_i에 비선형적으로 의존하게 만들어,
2체/3체 모두에서 실패했던 문제(존재하는 관계는 항상 이미 만족됨)를 원천적으로 피한다."""
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from network_base import build_growing_shell
import time

def run_fs(depth, buffer_radius, cutoff, A=1.0, p=8.0, q=4.0, eps=1.0):
    all_pts, vor, cluster_idx = build_growing_shell(buffer_radius=buffer_radius, depth=depth)
    points0 = all_pts[cluster_idx].copy()
    N = len(points0)
    r0 = np.sqrt(2)  # 최근접이웃 거리(참조 스케일)

    tree = cKDTree(points0)
    pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
    i_idx, j_idx = pairs[:,0], pairs[:,1]
    print(f"depth={depth}: N={N}, 이웃쌍={len(pairs)}")

    def phi(r):   return np.exp(-2*q*(r/r0 - 1))          # 밀도 기여함수
    def dphi(r):  return -2*q/r0*np.exp(-2*q*(r/r0 - 1))
    def V(r):     return A*np.exp(-p*(r/r0 - 1))           # 2체 반발항
    def dV(r):    return -A*p/r0*np.exp(-p*(r/r0 - 1))

    def energy_grad(x):
        pos = x.reshape(N,3)
        d = pos[i_idx]-pos[j_idx]
        r = np.linalg.norm(d, axis=1)
        unit = d/r[:,None]

        # 국소밀도 rho_i = sum_j phi(r_ij)
        rho = np.zeros(N)
        ph = phi(r)
        np.add.at(rho, i_idx, ph); np.add.at(rho, j_idx, ph)
        rho = np.maximum(rho, 1e-8)
        F = -eps*np.sqrt(rho)
        dFdrho = -eps*0.5/np.sqrt(rho)

        Vr = V(r)
        E = np.sum(F) + 0.5*np.sum(Vr)*2/2  # V는 각 pair 1회씩만 등장(이미 unique pairs)
        E = np.sum(F) + np.sum(Vr)

        # dE/dr_ij = (dFdrho[i]+dFdrho[j])*dphi(r) + dV(r)
        dEdr = (dFdrho[i_idx]+dFdrho[j_idx])*dphi(r) + dV(r)
        fvec = dEdr[:,None]*unit
        grad = np.zeros((N,3))
        np.add.at(grad, i_idx, fvec)
        np.add.at(grad, j_idx, -fvec)
        return E, grad.ravel()

    center_idx = np.argmin(np.linalg.norm(points0, axis=1))
    free_mask = np.ones(N, dtype=bool); free_mask[center_idx]=False
    x0 = points0.ravel().copy()

    def egf(xfree):
        xfull = x0.copy(); xfull.reshape(N,3)[free_mask]=xfree.reshape(-1,3)
        E,G = energy_grad(xfull)
        return E, G.reshape(N,3)[free_mask].ravel()

    xfree0 = points0[free_mask].ravel()
    E0,G0 = egf(xfree0)
    print(f"초기 에너지={E0:.4f}, 초기 그래디언트 노름={np.linalg.norm(G0):.6f} (0이면 이미 평형)")

    t0=time.time()
    res = minimize(egf, xfree0, jac=True, method='L-BFGS-B', options=dict(maxiter=5000, ftol=1e-16, gtol=1e-11))
    print(f"수렴={res.success}, 최종에너지={res.fun:.4f}, 소요={time.time()-t0:.1f}s")

    pos_final = points0.copy(); pos_final[free_mask]=res.x.reshape(-1,3)
    disp = pos_final-points0
    dist0 = np.linalg.norm(points0,axis=1)
    radial_unit = points0/np.where(dist0[:,None]==0,1,dist0[:,None])
    radial_disp = np.sum(disp*radial_unit,axis=1)

    coef = np.polyfit(dist0, radial_disp,1)
    resid = radial_disp - np.polyval(coef,dist0)
    print(f"선형계수={coef[0]:.5f}, 잔차표준편차={resid.std():.6f}")
    print("거리별 방사변위:")
    for d in sorted(set(np.round(dist0,3))):
        m = np.isclose(dist0,d,atol=1e-3)
        print(f"  거리={d:6.3f}: n={m.sum():4d}, 방사변위={radial_disp[m].mean(): .6f}")
    return dist0, radial_disp, resid

print("=== 최근접이웃만(cutoff=1.5), Finnis-Sinclair ===")
run_fs(4, 18, cutoff=1.5)
