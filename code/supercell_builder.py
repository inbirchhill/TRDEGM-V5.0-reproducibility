import numpy as np

def build_supercell(m):
    basis = np.load('basis_data.npy', allow_pickle=True).item()
    L = 2.0
    Lsup = L*m
    basis_items = list(basis.items())
    canon_pos = {fk: np.array(fk) for fk,_ in basis_items}

    atoms_pos, atoms_type, atoms_basis_key, cell_map = [], [], [], {}
    for i in range(m):
        for j in range(m):
            for k in range(m):
                off = np.array([L*i, L*j, L*k])
                for fk, b in basis_items:
                    p = canon_pos[fk] + off
                    key = tuple(np.round(p, 4))
                    idx_g = len(atoms_pos)
                    atoms_pos.append(p); atoms_type.append(b['type']); atoms_basis_key.append(fk)
                    cell_map[key] = idx_g
    atoms_pos = np.array(atoms_pos)
    atoms_type = np.array(atoms_type)
    Nsuper = len(atoms_pos)

    edges_i, edges_j, edges_n, edges_rtype = [], [], [], []
    seen = set()
    for gi in range(Nsuper):
        fk = atoms_basis_key[gi]; b = basis[fk]; pi = atoms_pos[gi]
        for rel, nbr_type in b['template']:
            p_wrapped = np.mod(pi+rel, Lsup)
            key = tuple(np.round(p_wrapped,4))
            gj = cell_map[key]
            n_vec = np.round((pi+rel-atoms_pos[gj])/Lsup).astype(int)
            pk = (min(gi,gj), max(gi,gj), tuple(n_vec if gi<gj else -n_vec))
            if pk in seen: continue
            seen.add(pk)
            edges_i.append(gi); edges_j.append(gj); edges_n.append(n_vec)
            ta,tb = atoms_type[gi], atoms_type[gj]
            rtype = "66" if (ta==6 and tb==6) else ("1212" if (ta==12 and tb==12) else "612")
            edges_rtype.append(rtype)
    edges_i=np.array(edges_i); edges_j=np.array(edges_j); edges_n=np.array(edges_n)

    q,p_exp,eps = 4.0,8.0,1.0
    RHO6,RHO12 = 10.0,20.0
    def dFdrho(rho): return -eps*0.5/np.sqrt(rho)
    def calibrate_A(ri,rj): return -(2*q/p_exp)*(dFdrho(ri)+dFdrho(rj))
    R0 = {"66":1.0,"612":0.8660254,"1212":1.4142136}
    Amap = {"66":calibrate_A(RHO6,RHO6),"612":calibrate_A(RHO6,RHO12),"1212":calibrate_A(RHO12,RHO12)}
    r0_edges = np.array([R0[t] for t in edges_rtype])
    A_edges = np.array([Amap[t] for t in edges_rtype])

    box0 = Lsup*np.eye(3)

    def phi(r,r0): return np.exp(-2*q*(r/r0-1))
    def dphi(r,r0): return -2*q/r0*np.exp(-2*q*(r/r0-1))
    def V(r,r0,A): return A*np.exp(-p_exp*(r/r0-1))
    def dV(r,r0,A): return -A*p_exp/r0*np.exp(-p_exp*(r/r0-1))

    def energy_grad(pos_flat):
        pos = pos_flat.reshape(Nsuper,3)
        shift = edges_n @ box0
        d = (pos[edges_i]-pos[edges_j]) - shift
        r = np.linalg.norm(d,axis=1)
        rho = np.zeros(Nsuper)
        np.add.at(rho, edges_i, phi(r,r0_edges))
        np.add.at(rho, edges_j, phi(r,r0_edges))
        F = -eps*np.sqrt(rho)
        E_embed = F.sum()
        E_pair = V(r,r0_edges,A_edges).sum()
        E = E_embed + E_pair
        grad = np.zeros((Nsuper,3))
        unit = d/r[:,None]
        dF_drho_i = dFdrho(rho[edges_i]); dF_drho_j = dFdrho(rho[edges_j])
        dEmbed_dr = (dF_drho_i+dF_drho_j)*dphi(r,r0_edges)
        dPair_dr = dV(r,r0_edges,A_edges)
        dE_dr = dEmbed_dr + dPair_dr
        fvec = dE_dr[:,None]*unit
        np.add.at(grad, edges_i, fvec)
        np.add.at(grad, edges_j, -fvec)
        return E, grad.ravel()

    x0 = atoms_pos.ravel()
    E0, G0 = energy_grad(x0)
    return dict(x0=x0, Nsuper=Nsuper, edges_i=edges_i, edges_j=edges_j, edges_n=edges_n,
                r0_edges=r0_edges, box0=box0, energy_grad=energy_grad, atoms_type=atoms_type,
                E0=E0, Gmax=np.abs(G0).max())

if __name__ == '__main__':
    import time
    for m in [6, 8]:
        t0 = time.time()
        ns = build_supercell(m)
        print(f"m={m}: 허브={ns['Nsuper']}, E0={ns['E0']:.4f}, |grad|max={ns['Gmax']:.2e}, "
              f"구성시간={time.time()-t0:.2f}s")
