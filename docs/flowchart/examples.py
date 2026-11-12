"""Worked examples for the flowchart page: 2D toy molecules run through a
Python transcription of min_rmsd (principal axes + LSAP per rank) and of the
centre walk (irmsd_walk.f90), recorded step by step as animation frames.

The toys are flat, so the 3D mirror (z -> -z) changes nothing and the four
180-degree orientations of min_rmsd become identity, x-flip, y-flip and
180-degree rotation in the plane.

Writes examples.json next to this file.
"""

import itertools
import json
from collections import deque
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------
# molecule helpers
# ---------------------------------------------------------------------


def ranks_1wl(elem, bonds):
    n = len(elem)
    nb = [[] for _ in range(n)]
    for i, j in bonds:
        nb[i].append(j)
        nb[j].append(i)
    col = [(elem[i], len(nb[i])) for i in range(n)]
    col = dense(col)
    while True:
        key = [(col[i], tuple(sorted(col[j] for j in nb[i]))) for i in range(n)]
        new = dense(key)
        if len(set(new)) == len(set(col)):
            return new, nb
        col = new


def dense(keys):
    order = {k: r + 1 for r, k in enumerate(sorted(set(keys)))}
    return [order[k] for k in keys]


def kabsch(X, Y):
    """R, cm, cr with R @ (x - cm) + cr ~ y (proper rotation, 2D)."""
    cm, cr = X.mean(0), Y.mean(0)
    A, B = X - cm, Y - cr
    s = np.sum(A[:, 0] * B[:, 1] - A[:, 1] * B[:, 0])
    c = np.sum(A[:, 0] * B[:, 0] + A[:, 1] * B[:, 1])
    t = np.arctan2(s, c)
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return R, cm, cr


def apply(R, cm, cr, P):
    return (P - cm) @ R.T + cr


def fit_rmsd(ref, mol, mapping):
    """mapping[i] = mol atom matched to ref atom i (all atoms)."""
    X = mol[mapping]
    R, cm, cr = kabsch(X, ref)
    Z = apply(R, cm, cr, mol)
    return float(np.sqrt(np.mean(np.sum((Z[mapping] - ref) ** 2, axis=1)))), Z


def violations(ref_bonds, mol_nb, mapping):
    return [(i, j) for i, j in ref_bonds if mapping[j] not in mol_nb[mapping[i]]]


def principal_frame(P):
    c = P.mean(0)
    Q = P - c
    w, v = np.linalg.eigh(Q.T @ Q)
    ax = v[:, 1]
    t = -np.arctan2(ax[1], ax[0])
    R = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
    return Q @ R.T


def lsap_per_rank(ref, mol, rank_r, rank_m):
    mapping = [-1] * len(ref)
    cost = 0.0
    for r in sorted(set(rank_r)):
        ri = [i for i in range(len(ref)) if rank_r[i] == r]
        mi = [j for j in range(len(mol)) if rank_m[j] == r]
        best = None
        for p in itertools.permutations(mi):
            c = sum(np.sum((ref[a] - mol[b]) ** 2) for a, b in zip(ri, p))
            if best is None or c < best[0]:
                best = (c, p)
        cost += best[0]
        for a, b in zip(ri, best[1]):
            mapping[a] = b
    return mapping, cost


# ---------------------------------------------------------------------
# centre walk (transcription of irmsd_walk.f90, heavy atoms only)
# ---------------------------------------------------------------------


class Walk:
    def __init__(self, ref, mol, rank, nb_r, nb_m, framethr=0.3):
        self.ref, self.mol, self.rank = ref, mol, rank
        self.nbr, self.nbm = nb_r, nb_m
        self.n = len(ref)
        self.framethr = framethr
        cls = {}
        for i, r in enumerate(rank):
            cls.setdefault(r, []).append(i)
        self.members = cls
        self.fixed = [len(cls[rank[i]]) == 1 for i in range(self.n)]
        ecc = []
        for i in range(self.n):
            d = self.bfs(i)
            ecc.append(max(d.values()))
        self.ecc = ecc
        c = min(range(self.n), key=lambda i: (ecc[i], len(cls[rank[i]]), rank[i], i))
        self.centre = c
        d = self.bfs(c)
        self.order = [i for i in sorted(d, key=lambda k: (d[k], list(d).index(k))) if not self.fixed[i]]
        self.bfsorder = list(d)
        self.order = [i for i in self.bfsorder if not self.fixed[i]]

    def bfs(self, s):
        d = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for v in self.nbr[u]:
                if v not in d:
                    d[v] = d[u] + 1
                    q.append(v)
        return d

    def start(self):
        return {i: i for i in range(self.n) if self.fixed[i]}

    def candidates(self, x, mp):
        inv = set(mp.values())
        anchor = next((p for p in self.nbr[x] if p in mp), None)
        pool = self.nbm[mp[anchor]] if anchor is not None else self.members[self.rank[x]]
        out = []
        for y in pool:
            if y in inv or self.rank[y] != self.rank[x] or len(self.nbm[y]) != len(self.nbr[x]):
                continue
            mx = [mp[z] for z in self.nbr[x] if z in mp]
            if any(z not in self.nbm[y] for z in mx):
                continue
            ny = sum(1 for z in self.nbm[y] if z in inv)
            if ny == len(mx):
                out.append(y)
        return out, anchor

    def defines_frame(self, mp):
        idx = list(mp)
        if len(idx) < 3:
            return False
        P = self.ref[idx]
        p0 = P[0]
        far = P[np.argmax(np.sum((P - p0) ** 2, axis=1))]
        e = (far - p0) / np.linalg.norm(far - p0)
        for p in P:
            v = p - p0
            v = v - np.dot(v, e) * e
            if np.linalg.norm(v) > self.framethr:
                return True
        return False

    def frame_of(self, mp):
        idx = list(mp)
        return kabsch(self.mol[[mp[i] for i in idx]], self.ref[idx])

    def walk_fixed(self, frame, force=None, log=None):
        """greedy walk in a fixed frame; force=(pos, k) takes the k-th nearest at pos"""
        mp = self.start()
        R, cm, cr = frame
        Z = apply(R, cm, cr, self.mol)
        decisions = []
        for pos, x in enumerate(self.order):
            c, anchor = self.candidates(x, mp)
            if not c:
                return None, decisions
            c.sort(key=lambda y: np.sum((Z[y] - self.ref[x]) ** 2))
            k = 0
            if force and force[0] == pos:
                k = min(force[1], len(c) - 1)
            if len(c) > 1 and anchor is not None:
                decisions.append((pos, len(c)))
            mp[x] = c[k]
            if log is not None:
                log.append((x, list(c), c[k], anchor))
        return mp, decisions


# ---------------------------------------------------------------------
# frame recording
# ---------------------------------------------------------------------


def r3(P):
    return [[round(float(a), 3), round(float(b), 3)] for a, b in P]


class Rec:
    def __init__(self, title, elem, ref, ref_bonds, mol_bonds, rank_r, rank_m, mol_labels):
        self.d = {
            "title": title, "elem": elem, "ref": r3(ref), "refBonds": ref_bonds,
            "molBonds": mol_bonds, "rankRef": rank_r, "rankMol": rank_m,
            "molLabels": mol_labels, "frames": [],
        }

    def frame(self, stage, caption, mol, **kw):
        f = {"stage": stage, "caption": caption, "mol": r3(mol)}
        f.update(kw)
        self.d["frames"].append(f)


def mapping_pairs(mp, status=None):
    out = []
    for i, j in (mp.items() if isinstance(mp, dict) else enumerate(mp)):
        if j is None or j < 0:
            continue
        out.append([int(i), int(j), (status or {}).get(i, "ok")])
    return out


# ---------------------------------------------------------------------
# example 1 + 2: octane-like chain
# ---------------------------------------------------------------------


def chain(n, angles, bond=1.5):
    P = [np.zeros(2)]
    d = 0.0
    for k in range(n - 1):
        d += angles[k]
        P.append(P[-1] + bond * np.array([np.cos(d), np.sin(d)]))
    return np.array(P)


def build_chain_examples():
    n = 8
    elem = ["C"] * n
    bonds = [(i, i + 1) for i in range(n - 1)]
    zig = [0.0] + [(-1) ** k * 1.0 for k in range(1, n - 1)]
    ref0 = chain(n, [0.5] + [(-1) ** k * 1.0 for k in range(1, n - 1)])
    rank, nb = ranks_1wl(elem, bonds)
    rng = np.random.default_rng(11)
    best = None
    for trial in range(40000):
        ang = [0.5] + [(-1) ** k * 1.0 + rng.normal(0, 1.3) for k in range(1, n - 1)]
        mol0 = chain(n, ang)
        # min_rmsd analogue
        A = principal_frame(ref0)
        B = principal_frame(mol0)
        res = []
        for o in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            Bo = B * np.array(o)
            mp, cost = lsap_per_rank(A, Bo, rank, rank)
            res.append((cost, o, mp))
        cost, o, mp = min(res, key=lambda t: t[0])
        viol = violations(bonds, nb, mp)
        r_lsap, _ = fit_rmsd(A, B * np.array(o), mp)
        # walk
        Bo = B * np.array(o)
        wk = Walk(A, Bo, rank, nb, nb)
        leaves = []
        for y in wk.members[rank[wk.order[0]]]:
            m0 = {wk.order[0]: y}
            mp2, _ = walk_from(wk, m0)
            if mp2:
                leaves.append(mp2)
        rw = min(fit_rmsd(A, Bo, [l[i] for i in range(n)])[0] for l in leaves)
        # self-avoiding: no two atoms too close
        D = np.linalg.norm(mol0[:, None] - mol0[None], axis=-1)
        far = np.abs(np.arange(n)[:, None] - np.arange(n)[None]) >= 2
        if D[far].min() < 1.8:
            continue
        if len(viol) >= 1 and r_lsap < rw - 0.1 and rw < 2.2:
            score = (rw - r_lsap) - 0.2 * abs(len(viol) - 2)
            if best is None or score > best[0]:
                best = (score, ang)
    ang = best[1]
    mol0 = chain(n, ang)
    return chain_frames(ref0, mol0, elem, bonds, rank, nb)


def walk_from(wk, m0):
    """complete a walk whose first order position is fixed by m0 (chains only need this)"""
    mp = wk.start()
    mp.update(m0)
    for x in wk.order[1:]:
        c, _ = wk.candidates(x, mp)
        if not c:
            return None, None
        mp[x] = c[0]
    return mp, None


def chain_frames(ref0, mol0, elem, bonds, rank, nb):
    n = len(elem)
    rng = np.random.default_rng(5)
    perm = rng.permutation(n)  # input order of mol: mol atom k is chain atom perm[k]
    inv = np.argsort(perm)
    th = 2.2
    Rr = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    A = principal_frame(ref0)
    mol_in = (mol0 - mol0.mean(0)) @ Rr.T + np.array([0.0, 0.0])
    mol_labels = [int(inv[i]) + 1 for i in range(n)]  # input index shown on chain atom i

    ex1 = Rec("Chain without unique atoms: min_rmsd", elem, A, bonds, bonds, rank, rank, mol_labels)
    side = mol_in + np.array([A[:, 0].max() - mol_in[:, 0].min() + 2.5, 0.0])
    ex1.frame("Input", "Reference (filled) and a second conformer (rings) of the same eight-carbon chain. "
              "The second structure comes in a different atom order: its labels are its input indices.", side,
              showRanks=False, labels="input")
    ex1.frame("Ranks", "Canonical ranks colour the atoms. The two chain ends share a rank, and so do each "
              "pair of atoms at equal distance from the ends: four ranks, two atoms each. "
              "No atom is unique, so min_rmsd cannot align on unique atoms.", side, showRanks=True, labels="input")
    B = principal_frame(mol0)
    ex1.frame("Principal axes", "Both structures are moved to their centre of mass and turned so that "
              "the long principal axis lies along x. This fixes the orientation only up to flips.", B,
              showRanks=True, labels="input", axis=True)
    res = []
    names = ["as aligned", "flipped about x", "flipped about y", "rotated 180 deg"]
    for name, o in zip(names, [(1, 1), (1, -1), (-1, 1), (-1, -1)]):
        Bo = B * np.array(o)
        mp, cost = lsap_per_rank(A, Bo, rank, rank)
        res.append((cost, o, mp, name))
        ex1.frame("LSAP per rank", f"Orientation '{name}': within each rank, the assignment with the smallest "
                  f"sum of squared distances (LSAP). Total cost {cost:.2f}.", Bo, showRanks=True,
                  labels="input", pairs=mapping_pairs(mp), metrics={"LSAP cost": round(cost, 2)}, axis=True)
    cost, o, mp, name = min(res, key=lambda t: t[0])
    Bo = B * np.array(o)
    viol = violations(bonds, nb, mp)
    r, Z = fit_rmsd(A, Bo, mp)
    st = {}
    for i, j in viol:
        st[i] = "bad"
        st[j] = "bad"
    ex1.frame("Lowest cost", f"The lowest total cost wins: '{name}'. Orientations that differ by a "
              "180-degree turn tie here, because reversing the whole chain is allowed by the ranks. Each rank "
              "was assigned on its own, without looking at bonds.", Bo, showRanks=True, labels="input",
              pairs=mapping_pairs(mp), metrics={"LSAP cost": round(cost, 2)})
    ex1.frame("Bond check", f"Fitted with this assignment the RMSD is {r:.2f}, but {len(viol)} bond(s) of the "
              "reference map onto atoms that are not bonded in the second structure (red). "
              "The low value comes from an impossible correspondence.", Z, showRanks=True, labels="input",
              pairs=mapping_pairs(mp, st), badBonds=[list(map(int, b)) for b in viol],
              metrics={"RMSD": round(r, 2), "broken bonds": len(viol)})

    # ----- example 2: centre walk on the same pair
    wk = Walk(A, Bo, rank, nb, nb)
    ex2 = Rec("Same pair: centre walk", elem, A, bonds, bonds, rank, rank, mol_labels)
    ex2.frame("Start", "The walk starts from min_rmsd's output frame. It keeps no assignment from it.",
              Bo, showRanks=True, labels="input")
    order_lbl = {int(x): k + 1 for k, x in enumerate(wk.order)}
    ex2.frame("Walk order", "Atoms are visited in BFS order from the graph centre (smallest eccentricity). "
              "Numbers on the reference show the order.", Bo, showRanks=True, labels="input",
              refOrder=order_lbl, current=int(wk.order[0]))
    root = wk.order[0]
    cands = wk.members[rank[root]]
    ex2.frame("No frame yet", "The first atom has no matched neighbour and nothing is matched yet, so there "
              "is no frame to judge distances in. Both atoms of its rank are tried, one after the other.",
              Bo, showRanks=True, labels="input", refOrder=order_lbl, current=int(root),
              candidates=[int(c) for c in cands])
    leaves = []
    for bi, y in enumerate(cands):
        mp2 = wk.start()
        mp2[root] = y
        ex2.frame(f"Start {bi + 1}", f"Start {bi + 1}: the centre atom is matched to atom "
                  f"{mol_labels[y]}.", Bo, showRanks=True, labels="input", refOrder=order_lbl,
                  pairs=mapping_pairs(mp2), current=int(root))
        for k, x in enumerate(wk.order[1:], start=2):
            c, anchor = wk.candidates(x, mp2)
            mp2[x] = c[0]
            if bi == 0:
                if k == 2:
                    txt = (f"Atom {k} in the order may only go to a free neighbour of its matched "
                           f"neighbour's partner, with the same rank. There is exactly one: atom "
                           f"{mol_labels[c[0]]}.")
                else:
                    txt = f"Atom {k}: again a single candidate, atom {mol_labels[c[0]]}. In a chain every step is forced."
                ex2.frame(f"Start {bi + 1}", txt, Bo, showRanks=True, labels="input", refOrder=order_lbl,
                          pairs=mapping_pairs(mp2), current=int(x), candidates=[int(v) for v in c])
        full = [mp2[i] for i in range(len(elem))]
        r2, Z2 = fit_rmsd(A, Bo, full)
        leaves.append((r2, full, Z2, y))
        ex2.frame(f"Start {bi + 1}", f"Complete. Kabsch fit on all atoms: RMSD {r2:.2f}. Every bond is kept, "
                  "because each step extended the assignment along a bond.", Z2, showRanks=True,
                  labels="input", pairs=mapping_pairs(full), metrics={"RMSD": round(r2, 2), "broken bonds": 0})
    r2, full, Z2, y = min(leaves, key=lambda t: t[0])
    ex2.frame("Selection", f"Both starts give {r2:.2f}: they differ by reversing the chain, a true symmetry. "
              f"min_rmsd's {r:.2f} broke {len(viol)} bonds, so refine_assignment keeps the walk result, "
              "the lowest RMSD among bond-preserving candidates.", Z2, showRanks=True, labels="input",
              pairs=mapping_pairs(full), metrics={"walk": round(r2, 2), "min_rmsd": round(r, 2),
                                                  "chosen": "walk"})
    return ex1.d, ex2.d


# ---------------------------------------------------------------------
# example 3: branched molecule, local search
# ---------------------------------------------------------------------


def arm(base, phi, psi, bond=1.5):
    p1 = base + bond * np.array([np.cos(phi), np.sin(phi)])
    p2 = p1 + bond * np.array([np.cos(phi + psi), np.sin(phi + psi)])
    return p1, p2


def build_branch_example():
    # 0 N centre, 1 C, 2 O, 3 Cl (unique tail); 4 a1, 5 a2, 6 b1, 7 b2 (two equivalent arms)
    elem = ["N", "C", "O", "Cl", "C", "C", "C", "C"]
    bonds = [(0, 1), (1, 2), (1, 3), (0, 4), (4, 5), (0, 6), (6, 7)]
    rank, nb = ranks_1wl(elem, bonds)
    n = len(elem)
    bonded = np.eye(n, dtype=bool)
    for i, j in bonds:
        bonded[i, j] = bonded[j, i] = True

    def make(ph, tau=0.0):
        """arms at angles ph; the C(=O)Cl tail turned by tau about N"""
        P = np.zeros((n, 2))
        R = np.array([[np.cos(tau), -np.sin(tau)], [np.sin(tau), np.cos(tau)]])
        P[1] = R @ [1.45, 0.0]
        P[2] = P[1] + 1.25 * (R @ [np.cos(1.0), np.sin(1.0)])
        P[3] = P[1] + 1.7 * (R @ [np.cos(-1.0), np.sin(-1.0)])
        P[4], P[5] = arm(P[0], ph[0], ph[1])
        P[6], P[7] = arm(P[0], ph[2], ph[3])
        return P

    ref = make([np.deg2rad(125), np.deg2rad(25), np.deg2rad(235), np.deg2rad(-25)])
    rng = np.random.default_rng(3)
    best = None
    for trial in range(30000):
        ph = [rng.uniform(1.4, 3.5), rng.uniform(-1.8, 1.8), rng.uniform(2.8, 4.9), rng.uniform(-1.8, 1.8)]
        tau = rng.uniform(-1.8, 1.8)
        mol = make(ph, tau)
        D = np.linalg.norm(mol[:, None] - mol[None], axis=-1)
        if D[~bonded].min() < 1.7:
            continue
        wk = Walk(ref, mol, rank, nb, nb)
        mp, dec = wk.walk_fixed(wk.frame_of(wk.start()))
        full = [mp[i] for i in range(n)]
        cur = full
        for _ in range(10):
            fr2 = kabsch(mol[cur], ref)
            mp2, _ = wk.walk_fixed(fr2)
            nxt = [mp2[i] for i in range(n)]
            if nxt == cur:
                break
            cur = nxt
        if cur != full:
            continue
        rp, _ = fit_rmsd(ref, mol, cur)
        mpa, _ = wk.walk_fixed(kabsch(mol[cur], ref), force=(dec[0][0], 1))
        ra, _ = fit_rmsd(ref, mol, [mpa[i] for i in range(n)])
        if ra < rp - 0.1 and (best is None or rp - ra > best[0]):
            best = (rp - ra, mol.copy())
    return branch_frames(ref, best[1], elem, bonds, rank, nb)


def branch_frames(ref, mol0, elem, bonds, rank, nb):
    n = len(elem)
    th = -1.1
    Rr = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    mol_in = (mol0 - mol0.mean(0)) @ Rr.T
    perm = np.random.default_rng(9).permutation(n)
    mol_labels = [int(np.argsort(perm)[i]) + 1 for i in range(n)]
    side = mol_in + np.array([ref[:, 0].max() - mol_in[:, 0].min() + 2.2, 0.0])
    ex = Rec("Branched molecule: local search", elem, ref, bonds, bonds, rank, rank, mol_labels)
    ex.frame("Input", "A tail of three different atoms (C, O, Cl) on a nitrogen, plus two identical "
             "two-carbon arms. In the second conformer the tail is turned and the arms point elsewhere.",
             side, showRanks=False, labels="elem")
    ex.frame("Ranks", "N, C, O and Cl are unique. The arms give two ranks with two atoms each: the inner "
             "carbons, and the outer ones.", side, showRanks=True, labels="elem")
    wk = Walk(ref, mol_in, rank, nb, nb)
    st = wk.start()
    fr = wk.frame_of(st)
    Z = apply(*fr, mol_in)
    ex.frame("Frame from unique atoms", "The four unique atoms are matched up front and define a frame: "
             "a Kabsch fit on them moves the second structure onto the reference.", Z, showRanks=True,
             labels="elem", pairs=mapping_pairs(st))
    log = []
    mp, dec = wk.walk_fixed(fr, log=log)
    first = log[0]
    x, c, chosen, anchor = first
    d = {int(y): round(float(np.linalg.norm(Z[y] - ref[x])), 2) for y in c}
    part = dict(st)
    part[x] = chosen
    ex.frame("Nearest candidate", f"First inner carbon: two candidates, at {d[c[0]]:.2f} and {d[c[1]]:.2f}. "
             "With a frame in place the walk takes the nearer one. The frame came from the turned tail "
             "only, so it misleads here.", Z, showRanks=True, labels="elem",
             pairs=mapping_pairs(part), current=int(x), candidates=[int(v) for v in c], dist=d)
    full = [mp[i] for i in range(n)]
    ex.frame("Forced", "The other inner carbon and both outer carbons have exactly one candidate each "
             "once that choice is made.", Z, showRanks=True, labels="elem", pairs=mapping_pairs(full))
    rg, Zg = fit_rmsd(ref, mol_in, full)
    ex.frame("Polish", f"Kabsch fit on all atoms: RMSD {rg:.2f}. Re-walking in this new frame picks the same "
             "candidates, so the alternation stops here.", Zg, showRanks=True, labels="elem",
             pairs=mapping_pairs(full), metrics={"RMSD": round(rg, 2)})
    fr2 = kabsch(mol_in[full], ref)
    mpa, _ = wk.walk_fixed(fr2, force=(dec[0][0], 1))
    alt = [mpa[i] for i in range(n)]
    ra, Za = fit_rmsd(ref, mol_in, alt)
    ex.frame("Local search", "Local search: take the other candidate at that choice and let the walk "
             "finish. The arms swap as a whole.", Zg, showRanks=True, labels="elem",
             pairs=mapping_pairs(alt, {4: "try", 5: "try", 6: "try", 7: "try"}), current=int(x))
    ex.frame("Accepted", f"Fitted, the swapped assignment gives RMSD {ra:.2f} instead of {rg:.2f}: the outer "
             "carbons fit far better. The change is accepted and the polish restarts from it; nothing "
             "improves further.", Za, showRanks=True, labels="elem", pairs=mapping_pairs(alt),
             metrics={"RMSD": round(ra, 2), "before": round(rg, 2)})
    return ex.d


if __name__ == "__main__":
    ex1, ex2 = build_chain_examples()
    ex3 = build_branch_example()
    out = Path(__file__).with_name("examples.json")
    out.write_text(json.dumps([ex1, ex2, ex3]))
    for e in (ex1, ex2, ex3):
        print(e["title"], len(e["frames"]), "frames")
        for f in e["frames"]:
            print("   ", f["stage"], "|", f.get("metrics", ""))
