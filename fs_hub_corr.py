import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from network_base import build_growing_shell, build_multicell_graph, classify_hub_types
from scipy.stats import pearsonr

def run_fs_hub(depth, buffer_radius, cutoff, A=1.0, p=8.0, q=4.0, eps=1.0):
    all_pts, vor, cluster_idx = build_growing_shell(buffer_radius=buffer_radius, depth=depth)
    G, face_seen, vertex_coord = build_multicell_graph(all_pts, vor, cluster_idx)
    hub_pos, hub_N, n_mixed = classify_hub_types(G)

    hub_ids = list(hub_pos.keys())
    points0 = np.array([hub_pos[h] for h in hub_ids])
    hub_types = np.array([hub_N[h] for h in hub_ids])  # 10 or 20
    N = len(points0)
    r0 = np.median([np.linalg.norm(points0[i]-points0[j])
                     for i in range(min(50,N)) for j in range(i+1,min(50,N))
                     if np.linalg.norm(points0[i]-points0[j])<3])

    tree = cKDTree(points0)
    pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
    if len(pairs)==0:
        print(f"cutoff={cutoff} 너무 작음, 결합쌍 없음"); return
    i_idx, j_idx = pairs[:,0], pairs[:,1]
    print(f"depth={depth}: 허브 {N}개, 결합 {len(pairs)}개, r0={r0:.3f}")

    def phi(r):   return np.exp(-2*q*(r/r0 - 1))
    def dphi(r):  return -2*q/r0*np.exp(-2*q*(r/r0 - 1))
    def V(r):     return A*np.exp(-p*(r/r0 - 1))
    def dV(r):    return -A*p/r0*np.exp(-p*(r/r0 - 1))
    def energy_grad(x):
        pos = x.reshape(N,3)
        d = pos[i_idx]-pos[j_idx]; r = np.linalg.norm(d, axis=1)
        r = np.maximum(r, 1e-6)
        unit = d/r[:,None]
        rho = np.zeros(N); ph = phi(r)
        np.add.at(rho, i_idx, ph); np.add.at(rho, j_idx, ph)
        rho = np.maximum(rho, 1e-8)
        F = -eps*np.sqrt(rho); dFdrho = -eps*0.5/np.sqrt(rho)
        Vr = V(r); E = np.sum(F) + np.sum(Vr)
        dEdr = (dFdrho[i_idx]+dFdrho[j_idx])*dphi(r) + dV(r)
        fvec = dEdr[:,None]*unit
        grad = np.zeros((N,3))
        np.add.at(grad, i_idx, fvec); np.add.at(grad, j_idx, -fvec)
        return E, grad.ravel()

    center_idx = np.argmin(np.linalg.norm(points0, axis=1))
    free_mask = np.ones(N, dtype=bool); free_mask[center_idx]=False
    x0 = points0.ravel().copy()
    def egf(xfree):
        xfull = x0.copy(); xfull.reshape(N,3)[free_mask]=xfree.reshape(-1,3)
        E,G_ = energy_grad(xfull)
        return E, G_.reshape(N,3)[free_mask].ravel()
    xfree0 = points0[free_mask].ravel()
    res = minimize(egf, xfree0, jac=True, method='L-BFGS-B', options=dict(maxiter=3000, ftol=1e-14, gtol=1e-10))
    pos_final = points0.copy(); pos_final[free_mask]=res.x.reshape(-1,3)
    disp = pos_final-points0
    disp_mag = np.linalg.norm(disp, axis=1)

    dist0 = np.linalg.norm(points0, axis=1)
    core_mask = (dist0 < buffer_radius*0.5) & (dist0 > 0.1)
    print(f"수렴={res.success}, 중심부표본={core_mask.sum()}")

    # 배위수(coordination) 자체와의 상관
    coord = np.zeros(N, dtype=int)
    np.add.at(coord, i_idx, 1); np.add.at(coord, j_idx, 1)
    if core_mask.sum() > 5:
        r_c, p_c = pearsonr(coord[core_mask], disp_mag[core_mask])
        print(f"  배위수-변위크기 r={r_c:.4f}, p={p_c:.4e}")
        r_t, p_t = pearsonr(hub_types[core_mask], disp_mag[core_mask])
        print(f"  허브타입(N6=10/N12=20)-변위크기 r={r_t:.4f}, p={p_t:.4e}")

        # N12 허브만 따로
        n12_mask = core_mask & (hub_types==20)
        if n12_mask.sum() > 5:
            r_n12, p_n12 = pearsonr(coord[n12_mask], disp_mag[n12_mask])
            print(f"  [N12만] 배위수-변위 r={r_n12:.4f}, p={p_n12:.4e}, 표본={n12_mask.sum()}")
        n6_mask = core_mask & (hub_types==10)
        if n6_mask.sum() > 5:
            r_n6, p_n6 = pearsonr(coord[n6_mask], disp_mag[n6_mask])
            print(f"  [N6만] 배위수-변위 r={r_n6:.4f}, p={p_n6:.4e}, 표본={n6_mask.sum()}")

print("=== 허브그래프 기준, depth=2 ===")
run_fs_hub(2, 14, cutoff=2.0)

print("=== 허브그래프 기준, depth=2, cutoff 좁게(1.8) ===")
run_fs_hub(2, 14, cutoff=1.8)
