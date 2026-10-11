"""
96허브 초격자(m=2, 12기저x8=96)에서 국소화된 단일입자 파동함수를
슈뢰딩거 방정식으로 시간발전시킨 뒤, 서로 다른 거리 r의 두 허브
사이 진폭 상관(|psi_i psi_j*|, "상호정보"의 대리지표)이 이미
확립된 유효 중력퍼텐셜 G(r)~r^-0.964와 같은 방식으로 감쇠하는지
대조한다. supercell_builder.py(96허브 생성)와 hexadecapole_kspace_
check.py의 홉핑 구성 방식을 재사용해 구성.
"""
import numpy as np

def build_supercell(m):
    basis = np.load('basis_data.npy', allow_pickle=True).item()
    L = 2.0
    Lsup = L*m
    basis_items = list(basis.items())
    canon_pos = {fk: np.array(fk) for fk,_ in basis_items}

    atoms_pos, atoms_basis_key = [], []
    for i in range(m):
        for j in range(m):
            for k in range(m):
                off = np.array([L*i, L*j, L*k])
                for fk, b in basis_items:
                    p = canon_pos[fk] + off
                    atoms_pos.append(p)
                    atoms_basis_key.append(fk)
    atoms_pos = np.array(atoms_pos)
    Nsuper = len(atoms_pos)

    # 좌표 -> 전역 인덱스 매핑(주기경계, mod Lsup)
    def wrap_key(p, nd=4):
        return tuple(np.round(np.mod(p, Lsup), nd))
    pos_to_idx = {wrap_key(p): i for i, p in enumerate(atoms_pos)}

    edges = []
    for gi in range(Nsuper):
        fk = atoms_basis_key[gi]
        b = basis[fk]
        pi = atoms_pos[gi]
        for rel, nbr_type in b['template']:
            p_nbr = np.mod(pi + rel, Lsup)
            gj = pos_to_idx[wrap_key(p_nbr)]
            if gi < gj:
                edges.append((gi, gj))
    return atoms_pos, edges, Nsuper

pos, edges, N = build_supercell(2)
print(f"초격자 구성: {N}개 허브, {len(edges)}개 결합")

# 홉핑 해밀토니안(단일입자, tight-binding, hopping=-1)
H = np.zeros((N, N), dtype=complex)
for i, j in edges:
    H[i, j] = -1.0
    H[j, i] = -1.0

# 국소화된 초기상태(원점 근방 허브 하나에 100% 집중)
origin_idx = np.argmin(np.linalg.norm(pos, axis=1))
psi0 = np.zeros(N, dtype=complex)
psi0[origin_idx] = 1.0

# 슈뢰딩거 시간발전 U(t)=exp(-iHt), 고유분해 사용
evals, evecs = np.linalg.eigh(H)
c0 = evecs.conj().T @ psi0

def psi_at(t):
    phase = np.exp(-1j*evals*t)
    return evecs @ (phase*c0)

# 원점으로부터의 실제거리(최근접 이미지, 주기경계 고려)
Lsup = 4.0
def min_image_dist(p, q):
    d = p - q
    d = d - Lsup*np.round(d/Lsup)
    return np.linalg.norm(d)

dists = np.array([min_image_dist(pos[origin_idx], pos[i]) for i in range(N)])

for T in [2.0, 5.0, 10.0, 20.0, 30.0]:
    psi_T = psi_at(T)
    prob = np.abs(psi_T)**2
    # "상호정보" 대리지표: 원점-i 사이 결합확률(0이면 완전 독립, 클수록 상관)
    # 여기서는 |psi_i|^2 자체를 원점과의 잔류 상관으로 사용(초기 원점집중에서
    # 얼마나 퍼졌는지가 곧 원점과 i의 "연결" 정도)
    r_bins = np.array([1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0])
    print(f"\n--- T={T} ---")
    vals = []
    for r in r_bins:
        mask = np.abs(dists - r) < 0.3
        if mask.sum() > 0:
            avg_p = prob[mask].mean()
            vals.append(avg_p)
            print(f"  r~{r}: 평균 |psi|^2 = {avg_p:.6f} (표본 {mask.sum()}개)")
    if len(vals) >= 2:
        print(f"  T={T} 전체 범위: {min(vals):.6f} ~ {max(vals):.6f}")
