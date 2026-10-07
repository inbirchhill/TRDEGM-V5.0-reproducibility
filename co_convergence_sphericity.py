# Copyright (c) 2026 Hyeongsu Kim (김형수)
# Licensed under CC BY-NC-SA 4.0 (저작자표시-비영리-동일조건변경허락 4.0)
# https://creativecommons.org/licenses/by-nc-sa/4.0/deed.ko
# Part of: "Rhombic Dodecahedral Energy Grid Model V4.0"

"""[V5.0 준비/4.0 편입, 항목B, 재설계] 1차 시도(격자점 중심의 볼록껍질)는
자명하게 매 depth마다 완전히 동일(psi=0.905, 변동계수=0 고정)했다 —
이는 §2.2가 논하는 "RD 고체 합집합의 실제 표면"이 아니라 "격자점
자체의 볼록껍질"(육팔면체류 성장수열의 수학적 성질)을 잰 것이었기
때문. 이번엔 실제 RD 다면체 정점들(vertex_coord)의 볼록껍질로
재계산해 §2.2가 말한 "진짜 표면"에 더 가깝게 만든다."""
import numpy as np
from scipy.spatial import ConvexHull
from network_base import build_growing_shell, build_multicell_graph

def sphericity_metrics_v2(depth, buffer_radius):
    all_pts, vor, cluster_idx = build_growing_shell(buffer_radius=buffer_radius, depth=depth)
    G, face_seen, vertex_coord = build_multicell_graph(all_pts, vor, cluster_idx)
    verts = np.array(list(vertex_coord.values()))
    hull = ConvexHull(verts)
    V = hull.volume
    A = hull.area
    psi = (np.pi**(1/3) * (6*V)**(2/3)) / A
    centroid = verts.mean(axis=0)
    hull_pts = verts[hull.vertices]
    radii = np.linalg.norm(hull_pts - centroid, axis=1)
    cv_radius = radii.std() / radii.mean()
    return dict(depth=depth, n_verts=len(verts), n_hull_verts=len(hull.vertices),
                V=V, A=A, psi=psi, cv_radius=cv_radius,
                r_mean=radii.mean(), r_min=radii.min(), r_max=radii.max())

if __name__ == '__main__':
    import time
    results = []
    for depth, buf in [(0,4),(1,6),(2,8),(3,10),(4,12)]:
        t0=time.time()
        m = sphericity_metrics_v2(depth, buf)
        results.append(m)
        print(f"depth={depth}: 전체정점={m['n_verts']:5d} 껍질정점={m['n_hull_verts']:4d}  "
              f"구형도(psi)={m['psi']:.5f}  반지름변동계수={m['cv_radius']:.5f}  "
              f"r_min/r_max={m['r_min']:.3f}/{m['r_max']:.3f}  t={time.time()-t0:.1f}s", flush=True)

    print("\n=== 요약 ===")
    print("psi:", [round(m['psi'],5) for m in results])
    print("cv_radius:", [round(m['cv_radius'],5) for m in results])
