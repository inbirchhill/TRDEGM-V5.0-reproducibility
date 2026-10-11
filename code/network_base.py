"""공통 네트워크 구성 함수 (프로젝트 network_base.txt 그대로 복사)"""
import numpy as np
from scipy.spatial import Voronoi
import networkx as nx

def fcc_points(shell_radius):
    pts = []
    r = shell_radius
    for i in range(-r, r + 1):
        for j in range(-r, r + 1):
            for k in range(-r, r + 1):
                if (i + j + k) % 2 == 0:
                    pts.append((i, j, k))
    return np.array(pts, dtype=float)

def key(p, nd=6):
    return tuple(np.round(p, nd))

def build_growing_shell(buffer_radius, depth):
    all_pts = fcc_points(buffer_radius)
    idx_by_key = {key(p): i for i, p in enumerate(all_pts)}
    dists0 = np.linalg.norm(all_pts, axis=1)
    center_idx = np.where(dists0 < 1e-9)[0]
    cluster = set(center_idx.tolist())
    nn_offsets = []
    for i in (-1,0,1):
        for j in (-1,0,1):
            for k in (-1,0,1):
                if (i,j,k)==(0,0,0): continue
                if (i+j+k)%2==0 and i*i+j*j+k*k==2:
                    nn_offsets.append(np.array([i,j,k],dtype=float))
    for _ in range(depth+1):
        new_cluster = set(cluster)
        for pidx in cluster:
            p = all_pts[pidx]
            for off in nn_offsets:
                nk = key(p+off)
                if nk in idx_by_key:
                    new_cluster.add(idx_by_key[nk])
        cluster = new_cluster
    cluster_idx = np.array(sorted(cluster))
    vor = Voronoi(all_pts)
    return all_pts, vor, cluster_idx

def build_multicell_graph(all_pts, vor, cluster_idx):
    face_seen, vertex_coord = {}, {}
    for pidx in cluster_idx:
        for ridge_pts, ridge_verts in zip(vor.ridge_points, vor.ridge_vertices):
            if pidx not in ridge_pts or -1 in ridge_verts: continue
            other = ridge_pts[0] if ridge_pts[1]==pidx else ridge_pts[1]
            d = np.linalg.norm(all_pts[pidx]-all_pts[other])
            if not np.isclose(d, np.sqrt(2), atol=1e-6): continue
            vcoords = vor.vertices[ridge_verts]
            if len(vcoords)!=4: continue
            fc = vcoords.mean(axis=0)
            fkey = key(fc)
            if fkey in face_seen: continue
            vkeys=[key(v) for v in vcoords]
            for vk,vc in zip(vkeys,vcoords): vertex_coord[vk]=vc
            face_seen[fkey]=(fc,vkeys,tuple(sorted([pidx,other])))
    G = nx.Graph()
    for fkey,(fc,vkeys,pair) in face_seen.items():
        fnode=("F",fkey); G.add_node(fnode,pos=fc,kind="face")
        vcs=np.array([vertex_coord[vk] for vk in vkeys])
        dfc=np.linalg.norm(vcs-fc,axis=1)
        order=np.argsort(dfc)
        short_pair=[vkeys[order[0]],vkeys[order[1]]]
        for vk in vkeys:
            vnode=("V",vk)
            if vnode not in G: G.add_node(vnode,pos=vertex_coord[vk],kind="vertex")
            etype = "short" if vk in short_pair else "long"
            w = np.linalg.norm(vertex_coord[vk]-fc)
            G.add_edge(fnode,vnode,etype=etype,weight=w)
    return G, face_seen, vertex_coord

def classify_hub_types(G):
    hub_pos, hub_N = {}, {}
    n_mixed = 0
    for node,data in G.nodes(data=True):
        if data.get("kind")!="vertex": continue
        etypes=set(G.edges[node,f]["etype"] for f in G.neighbors(node))
        if etypes=={"short"}: hub_N[node]=6
        elif etypes=={"long"}: hub_N[node]=12
        else:
            n_mixed += 1
            continue
        hub_pos[node]=data["pos"]
    return hub_pos, hub_N, n_mixed

def build_hub_graph(G, hub_pos, hub_N):
    H = nx.Graph()
    for hn in hub_pos:
        H.add_node(hn, pos=hub_pos[hn], N=hub_N[hn])
    for fnode, data in G.nodes(data=True):
        if data.get("kind")!="face": continue
        neighbors = [v for v in G.neighbors(fnode) if v in hub_pos]
        for i in range(len(neighbors)):
            for j in range(i+1, len(neighbors)):
                a,b = neighbors[i], neighbors[j]
                d = np.linalg.norm(hub_pos[a]-hub_pos[b])
                if H.has_edge(a,b): continue
                H.add_edge(a,b, dist=d)
    return H
