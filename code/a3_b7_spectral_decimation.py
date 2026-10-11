"""로드맵 8부 A-3(B7): depth=시간 치환, 스펙트럼 재규격화/데시메이션 재현성확인.
기존 예비증거(spectral_decimation_check.txt: depth0/1/2 고유값분포 비교,
schur_decimation.txt: Schur complement 유효라플라시안 비교)에 더해, 오늘 만든
Bloch L(k) 인프라(§macro/3부-8에서 이미 확립)로 "진짜 무한주기격자의 상태밀도
(DOS)"를 직접 계산하고, 유한 depth(0,1,2)의 라플라시안 고유값분포가 depth가
커질수록 이 무한격자 DOS로 수렴하는가를 독립적인 방법으로 확인한다.
"""
import numpy as np
import networkx as nx
from network_base import (build_growing_shell, build_multicell_graph,
                            classify_hub_types, build_hub_graph)

# ---- (1) 유한 depth 라플라시안 스펙트럼 ----
def get_laplacian_eigs(depth, buffer_radius):
    all_pts, vor, cluster_idx = build_growing_shell(buffer_radius=buffer_radius, depth=depth)
    G, face_seen, vertex_coord = build_multicell_graph(all_pts, vor, cluster_idx)
    hub_pos, hub_N, n_mixed = classify_hub_types(G)
    H = build_hub_graph(G, hub_pos, hub_N)
    L = nx.laplacian_matrix(H).toarray().astype(float)
    eigvals = np.linalg.eigvalsh(L)
    return eigvals, H.number_of_nodes()

eig0, n0 = get_laplacian_eigs(0, 6)
eig1, n1 = get_laplacian_eigs(1, 9)
eig2, n2 = get_laplacian_eigs(2, 12)
print(f"depth0: {n0}노드, depth1: {n1}노드, depth2: {n2}노드")

# ---- (2) 무한주기격자 진짜 DOS (Bloch L(k), 12기저원자, k메시 조밀) ----
basis = np.load('basis_data.npy', allow_pickle=True).item()
basis_items = list(basis.items())
n_basis = len(basis_items)
L_const = 2.0
canon_pos = {fk: np.array(fk) for fk, _ in basis_items}
fk_list = [fk for fk, _ in basis_items]
fk_idx = {fk: i for i, fk in enumerate(fk_list)}
def frac_key(p, nd=4): return tuple(np.round(np.mod(p, L_const), nd))

hoppings = []; degree = np.zeros(n_basis)
for a_idx, (fk_a, b) in enumerate(basis_items):
    pa = canon_pos[fk_a]
    degree[a_idx] = len(b['template'])
    for rel, nbr_type in b['template']:
        p_nbr = pa + rel
        fk_b = frac_key(p_nbr)
        b_idx = fk_idx[fk_b]
        pb_canon = canon_pos[fk_b]
        R_int = np.round((p_nbr - pb_canon) / L_const).astype(int)
        hoppings.append((a_idx, b_idx, R_int))

def build_Lk(k):
    Hk = np.zeros((n_basis, n_basis), dtype=complex)
    for a_idx, b_idx, R_int in hoppings:
        Hk[a_idx, b_idx] += -1.0 * np.exp(1j*np.dot(k, R_int*L_const))
    return np.diag(degree).astype(complex) + Hk

Nk = 20
ks = (2*np.pi/(Nk*L_const)) * (np.arange(Nk) - Nk//2)
bloch_eigs = []
for kx in ks:
    for ky in ks:
        for kz in ks:
            Lk = build_Lk(np.array([kx, ky, kz]))
            bloch_eigs.extend(np.linalg.eigvalsh(Lk).tolist())
bloch_eigs = np.array(bloch_eigs)
print(f"\n무한격자 Bloch DOS 표본: {len(bloch_eigs)}개 고유값(Nk={Nk}^3 k점 x {n_basis}밴드)")
print(f"Bloch 고유값 범위: {bloch_eigs.min():.4f}~{bloch_eigs.max():.4f}")

# ---- (3) 정규화된 DOS 비교(코모고로프-스미르노프 거리로 depth별 Bloch극한과의 근접도 정량화) ----
from scipy.stats import ks_2samp
print("\n=== depth별 유한 스펙트럼 vs 무한격자 Bloch DOS: KS(콜모고로프-스미르노프) 거리 ===")
print("(KS거리가 depth 증가에 따라 감소하면 '무한격자 극한으로 수렴'을 정량 지지)")
for name, eig in [("depth0", eig0), ("depth1", eig1), ("depth2", eig2)]:
    stat, pval = ks_2samp(eig, bloch_eigs)
    print(f"  {name} (n={len(eig)}): KS거리={stat:.4f}, p={pval:.4e}")

print("\n=== 참고: 각 분포의 기초통계량 ===")
for name, eig in [("depth0", eig0), ("depth1", eig1), ("depth2", eig2), ("Bloch(무한)", bloch_eigs)]:
    print(f"  {name}: mean={eig.mean():.4f}, median={np.median(eig):.4f}, std={eig.std():.4f}")

print("\n=== 재현성 재확인 ===")
eig0_b, _ = get_laplacian_eigs(0, 6)
print("depth0 재계산 최대차이:", np.abs(np.sort(eig0)-np.sort(eig0_b)).max(), "(0이어야 정상)")
