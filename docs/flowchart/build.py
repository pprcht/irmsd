"""Build the iRMSD flowchart page (irmsd-flowchart.html) from node/edge specs.

Nodes are placed by centre on a hand-picked grid; edges are orthogonal
polylines between node sides, with optional waypoints and a label.
"""

import re
from html import escape
from pathlib import Path

LH = 16  # line height
PADY = 9
CW_SANS = 6.7  # approx. glyph widths at 12.5px
CW_MONO = 7.4

# inline formatting: code spans (`...`, calls, snake_case names, a few
# known identifiers) in monospace, the structure names ref/mol in bold
CODE_WORDS = r"molatomsort|checkranks|fallbackranks|iinversion|nunique|rinv"
TOKEN = re.compile(
    r"`[^`]+`"
    r"|[A-Za-z_][A-Za-z0-9_]*\((?!s\))[^()]*\)"
    r"|\b[a-z][a-z0-9]*(?:_[a-z0-9]+)+\b"
    rf"|\b(?:{CODE_WORDS})\b"
    r"|\b(?:ref|mol)\b"
)


def segments(line):
    """split a line into (text, kind) with kind in plain / code / bold;
    a line starting with a single backtick is code as a whole"""
    if line.startswith("`") and line.count("`") == 1:
        return [(line[1:], "code")]
    out, pos = [], 0
    for m in TOKEN.finditer(line):
        if m.start() > pos:
            out.append((line[pos:m.start()], "plain"))
        t = m.group(0)
        if t.startswith("`"):
            out.append((t[1:-1], "code"))
        elif t in ("ref", "mol"):
            out.append((t, "bold"))
        else:
            out.append((t, "code"))
        pos = m.end()
    if pos < len(line):
        out.append((line[pos:], "plain"))
    return out


def seg_width(segs):
    w = {"plain": CW_SANS, "bold": CW_SANS * 1.07, "code": CW_MONO}
    return sum(len(t) * w[k] for t, k in segs)


def rich(segs):
    cls = {"plain": None, "bold": "b", "code": "code"}
    out = []
    for t, k in segs:
        if cls[k]:
            out.append(f'<tspan class="{cls[k]}">{escape(t)}</tspan>')
        else:
            out.append(escape(t))
    return "".join(out)


# page order of the sections (ids in template.html), and link labels
ORDER = ["call", "pipeline", "ranks", "min-rmsd", "refine", "walk", "ex"]
SECTION = {"call": "call", "pipe": "pipeline", "rank": "ranks", "min": "min-rmsd",
           "refine": "refine", "walk": "walk"}
TITLES = {
    "call": "Go to section 0, what a get_irmsd call decides",
    "pipeline": "Go to section 1, one comparison end to end",
    "ranks": "Go to section 2, atom ranks",
    "min-rmsd": "Go to section 3, min_rmsd",
    "refine": "Go to section 4, bond-preserving selection",
    "walk": "Go to section 5, centre walk",
    "ex-a": "Go to worked example A, min_rmsd on a chain",
    "ex-b": "Go to worked example B, centre walk on the same chain",
    "ex-c": "Go to worked example C, local search",
}


class Fig:
    def __init__(self, fid, w, h, label):
        self.fid, self.w, self.h, self.label = fid, w, h, label
        self.nodes = {}
        self.parts = []
        self.deferred = []
        self.nodeparts = {}
        self.links = {}

    # ---- nodes -------------------------------------------------------
    def node(self, key, cx, cy, lines, kind="proc", w=None, new=False):
        """kind: proc | dec | term | note.  lines: str; '`' prefix = mono."""
        lines = [lines] if isinstance(lines, str) else lines
        widths = [seg_width(segments(t)) for t in lines]
        need = max(widths) + (44 if kind == "dec" else 26)
        w = w or max(need, 120)
        if need > w + 1:
            print(f"warning: {self.fid}:{key} text {need:.0f} > box {w}")
        h = len(lines) * LH + 2 * PADY
        x0, y0 = cx - w / 2, cy - h / 2
        cls = f"{kind}{' new' if new else ''}"
        start = len(self.parts)
        if kind == "dec":
            k = 14
            pts = [
                (x0 + k, y0), (x0 + w - k, y0), (x0 + w, cy), (x0 + w - k, y0 + h),
                (x0 + k, y0 + h), (x0, cy),
            ]
            self.parts.append(
                f'<polygon class="{cls}" points="{" ".join(f"{a:.1f},{b:.1f}" for a, b in pts)}"/>'
            )
        else:
            r = h / 2 if kind == "term" else 5
            self.parts.append(
                f'<rect class="{cls}" x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{r:.1f}"/>'
            )
        ty = cy - (len(lines) - 1) * LH / 2 + 4.3
        for i, t in enumerate(lines):
            self.parts.append(
                f'<text class="t" x="{cx:.1f}" y="{ty + i * LH:.1f}" text-anchor="middle">{rich(segments(t))}</text>'
            )
        self.nodes[key] = (cx, cy, w, h, kind)
        self.nodeparts[key] = self.parts[start:]
        del self.parts[start:]

    def link(self, key, target):
        """make node `key` a link to the section/player with id `target`"""
        self.links[key] = target

    def anchor(self, key, side):
        cx, cy, w, h, _ = self.nodes[key]
        return {
            "t": (cx, cy - h / 2), "b": (cx, cy + h / 2),
            "l": (cx - w / 2, cy), "r": (cx + w / 2, cy),
        }[side]

    # ---- edges -------------------------------------------------------
    def edge(self, *args, **kw):
        self.deferred.append((args, kw))

    def _edge(self, a, sa, b, sb, via=(), label=None, lpos=None, anchor="start", dashed=False):
        p0 = self.anchor(a, sa)
        p1 = self.anchor(b, sb)
        pts = [p0, *via, p1]
        cls = "e dash" if dashed else "e"
        self.parts.append(
            f'<polyline class="{cls}" points="{" ".join(f"{x:.1f},{y:.1f}" for x, y in pts)}" marker-end="url(#ah-{self.fid})"/>'
        )
        if label:
            if lpos is None:
                (x0, y0), (x1, y1) = pts[0], pts[1]
                lpos = ((x0 + x1) / 2 + 6, (y0 + y1) / 2 + 4) if x0 == x1 else ((x0 + x1) / 2, y0 - 6)
                anchor = "start" if x0 == x1 else "middle"
            self.parts.append(
                f'<text class="lbl" x="{lpos[0]:.1f}" y="{lpos[1]:.1f}" text-anchor="{anchor}">{rich(segments(label))}</text>'
            )

    def marker(self, x, y, num):
        self.parts.append(f'<circle class="mk" cx="{x:.1f}" cy="{y:.1f}" r="10"/>')
        self.parts.append(
            f'<text class="mkt" x="{x:.1f}" y="{y + 4.3:.1f}" text-anchor="middle">{num}</text>'
        )

    def text(self, x, y, s, cls="lbl", anchor="start"):
        body = escape(s) if cls == "grpt" else rich(segments(s))
        self.parts.append(
            f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}">{body}</text>'
        )

    def frame(self, x, y, w, h, title):
        self.parts.append(
            f'<rect class="grp" x="{x}" y="{y}" width="{w}" height="{h}" rx="8"/>'
        )
        self.text(x + 12, y + 18, title, "grpt")

    def svg(self):
        for args, kw in self.deferred:
            self._edge(*args, **kw)
        self.deferred = []
        defs = (
            f'<defs><marker id="ah-{self.fid}" viewBox="0 0 10 10" refX="9" refY="5" '
            f'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            f'<path class="ah" d="M0,1 L9,5 L0,9 z"/></marker></defs>'
        )
        # groups first (background), then edges, then nodes: sort by type
        grp = [p for p in self.parts if 'class="grp' in p or 'class="grpt' in p]
        edges = [p for p in self.parts if "<polyline" in p or 'class="lbl"' in p]
        rest = [p for p in self.parts if p not in grp and p not in edges]
        nodes = []
        here = ORDER.index(SECTION[self.fid])
        for key, parts in self.nodeparts.items():
            if key not in self.links:
                nodes.extend(parts)
                continue
            target = self.links[key]
            cx, cy, w, h, kind = self.nodes[key]
            down = ORDER.index(target.split("-")[0] if target.startswith("ex-") else target) > here
            bx = cx + w / 2 - (14 if kind == "dec" else 4)
            by = cy - h / 2 + 4
            badge = (
                f'<circle class="lkb" cx="{bx:.1f}" cy="{by:.1f}" r="8"/>'
                f'<path class="lka" d="{arrow(bx, by, down)}"/>'
            )
            nodes.append(
                f'<a href="#{target}" class="lnk" data-target="{target}" '
                f'aria-label="{escape(TITLES[target])}">' + "".join(parts) + badge + "</a>"
            )
        body = "\n".join(grp + edges + rest + nodes)
        role = 'role="group"' if self.links else 'role="img"'
        return (
            f'<svg viewBox="0 0 {self.w} {self.h}" {role} aria-label="{escape(self.label)}" '
            f'xmlns="http://www.w3.org/2000/svg">{defs}\n{body}\n</svg>'
        )


def arrow(x, y, down):
    """small chevron inside the link badge, pointing to where the target is"""
    s = 1 if down else -1
    return f"M{x - 3.5:.1f},{y - 1.5 * s:.1f} L{x:.1f},{y + 2 * s:.1f} L{x + 3.5:.1f},{y - 1.5 * s:.1f}"


# =====================================================================
# Figure 0: what one get_irmsd call decides on its own
# =====================================================================
f0 = Fig("call", 860, 1340, "Automatic decisions in one get_irmsd call, from input to result")
L, M, R = 125, 400, 700
f0.node("s", M, 40, ["`get_irmsd(ref, mol)", "any atom order; iinversion = 0 (auto)"], "term", w=320)
f0.node("ids", M, 120, ["canonical IDs supplied", "with the structures?"], "dec", w=230)
f0.node("idy", R, 120, ["use the supplied IDs", "(recomputed if inconsistent)"], w=230)
f0.node("rr", M, 205, ["canonical ranks for ref"], w=230)
f0.node("same", M, 285, ["mol has ref's atom", "order and bonds?"], "dec", w=210, new=True)
f0.node("cp", R, 285, ["copy ref's ranks"], w=200, new=True)
f0.node("rm", M, 365, ["canonical ranks for mol"], w=230)
f0.node("cons", M, 450, ["ranks consistent?", "(same atoms per rank)"], "dec", w=220)
f0.node("fb", L, 450, ["fall back to", "element types"], w=160)
f0.node("topo", M, 535, ["same composition?"], "dec", w=200)
f0.node("bad", R, 535, ["warning; plain RMSD", "without reordering"], "term", w=210)
f0.node("st", M, 625, ["reference has a", "stereocentre?"], "dec", w=200)
f0.node("nomir", L, 625, ["mirror image", "not tried"], w=160)
f0.node("mir", R, 625, ["mirror image", "also tried"], w=160)
f0.node("uq", M, 715, [">= 3 atoms with", "a unique rank?"], "dec", w=200)
f0.node("ua", L, 805, ["align on the unique", "atoms; LSAP per rank"], w=200)
f0.node("pa", R, 805, ["principal axes, 4/8/16", "orientations; LSAP per rank"], w=230)
f0.node("ca", M, 895, ["candidate A: that assignment", "+ terminal-atom refinement"], w=290, new=True)
f0.node("clu", M, 980, ["largest bonded fragment", ">= half the heavy atoms?"], "dec", w=240, new=True)
f0.node("conly", L, 1065, ["cluster:", "candidate A only"], w=170, new=True)
f0.node("walk", M, 1065, ["centre walk (+ mirror if", "tried) -> candidates B, C"], w=260, new=True)
f0.node("val", M, 1155, ["any candidate keeps", "every bond?"], "dec", w=220, new=True)
f0.node("keepa", R, 1155, ["keep A"], w=140, new=True)
f0.node("low", M, 1235, ["lowest RMSD among", "the bond-preserving ones"], w=250, new=True)
f0.node("out", M, 1310, ["iRMSD + mol reordered and", "aligned onto ref"], "term", w=290)
f0.edge("s", "b", "ids", "t")
f0.edge("ids", "r", "idy", "l", label="yes")
f0.edge("ids", "b", "rr", "t", label="no")
f0.edge("rr", "b", "same", "t")
f0.edge("same", "r", "cp", "l", label="yes")
f0.edge("same", "b", "rm", "t", label="no")
f0.edge("idy", "b", "cons", "r", via=[(R, 450)])
f0.edge("cp", "b", "cons", "r", via=[(R, 450)])
f0.edge("rm", "b", "cons", "t")
f0.edge("cons", "l", "fb", "r", label="no")
f0.edge("cons", "b", "topo", "t", label="yes")
f0.edge("fb", "b", "topo", "l", via=[(L, 535)])
f0.edge("topo", "r", "bad", "l", label="no")
f0.edge("topo", "b", "st", "t", label="yes")
f0.edge("st", "l", "nomir", "r", label="yes")
f0.edge("st", "r", "mir", "l", label="no")
f0.edge("nomir", "b", "uq", "t", via=[(L, 670), (M, 670)])
f0.edge("mir", "b", "uq", "t", via=[(R, 670), (M, 670)])
f0.edge("uq", "l", "ua", "t", via=[(L, 715)], label="yes", lpos=(260, 708), anchor="middle")
f0.edge("uq", "r", "pa", "t", via=[(R, 715)], label="no", lpos=(540, 708), anchor="middle")
f0.edge("ua", "b", "ca", "t", via=[(L, 855), (M, 855)])
f0.edge("pa", "b", "ca", "t", via=[(R, 855), (M, 855)])
f0.edge("ca", "b", "clu", "t")
f0.edge("clu", "l", "conly", "t", via=[(L, 980)], label="no", lpos=(255, 973), anchor="middle")
f0.edge("clu", "b", "walk", "t", label="yes")
f0.edge("conly", "b", "val", "l", via=[(L, 1155)])
f0.edge("walk", "b", "val", "t")
f0.edge("val", "b", "low", "t", label="yes")
f0.edge("val", "r", "keepa", "l", label="no")
f0.edge("low", "b", "out", "t")
f0.edge("keepa", "b", "out", "r", via=[(R, 1310)])
f0.marker(M + 132, 962, 1)
f0.marker(M + 122, 432, 2)
f0.marker(M + 112, 607, 3)
f0.text(R, 592, "iinversion = 1 / 2 overrides this", "lbl", "middle")

# =====================================================================
# Figure 1: pipeline of one pairwise comparison
# =====================================================================
f1 = Fig("pipe", 900, 760, "Pipeline of one pairwise iRMSD comparison from input to output")
X = 330
f1.node("in", X, 40, ["ref, mol  (atomic numbers, xyz)", "optional canonical IDs per atom"], "term", w=300)
f1.node("axis", X, 120, ["`axis(ref)", "ref to centre of mass, principal axes"], w=300)
f1.node("ranks", X, 205, ["`setup_irmsd_ranks", "ranks r1, r2 + stereo flag", "see figure 2"], w=300)
f1.node("inv", X, 300, ["iinversion?"], "dec", w=170)
f1.node("inv0", 110, 385, ["0: auto", "mirror tried if ref", "has no stereocentre"], w=170)
f1.node("inv1", X, 385, ["1: mirror always", "tried"], w=170)
f1.node("inv2", 550, 385, ["2: never mirror"], w=170)
f1.node("min", X, 480, ["`min_rmsd", "ranks -> orientation -> LSAP per rank", "see figure 3"], w=300)
f1.node("io", X, 565, ["io = 0?"], "dec", w=150)
f1.node("topo", 690, 565, ["topology mismatch", "plain quaternion RMSD", "or huge value"], w=220)
f1.node("ref", X, 650, ["`refine_assignment", "bond-preserving candidates", "see figures 4 and 5"], w=300, new=True)
f1.node("out", X, 725, ["iRMSD (Angstrom) + mol reordered", "and aligned onto ref"], "term", w=300)
f1.edge("in", "b", "axis", "t")
f1.edge("axis", "b", "ranks", "t")
f1.edge("ranks", "b", "inv", "t")
f1.edge("inv", "l", "inv0", "t", via=[(110, 300)])
f1.edge("inv", "b", "inv1", "t")
f1.edge("inv", "r", "inv2", "t", via=[(550, 300)])
f1.edge("inv0", "b", "min", "l", via=[(110, 480)])
f1.edge("inv1", "b", "min", "t")
f1.edge("inv2", "b", "min", "r", via=[(550, 480)])
f1.edge("min", "b", "io", "t")
f1.edge("io", "r", "topo", "l", label="no")
f1.edge("io", "b", "ref", "t", label="yes")
f1.edge("ref", "b", "out", "t")
f1.edge("topo", "b", "out", "r", via=[(690, 725)])
f1.frame(700, 90, 185, 150, "Sorters (ensembles)")
f1.text(712, 132, "same core per pair:", "lbl")
f1.text(712, 150, "ranks once for structure 1", "lbl")
f1.text(712, 168, "(or per structure, copied", "lbl")
f1.text(712, 186, "when bonds and order match),", "lbl")
f1.text(712, 204, "then min_rmsd +", "lbl")
f1.text(712, 222, "refine_assignment", "lbl")

# =====================================================================
# Figure 2: rank setup
# =====================================================================
f2 = Fig("rank", 900, 820, "Decision tree of setup_irmsd_ranks")
L, R = 220, 640
f2.node("rid", L, 40, ["ref has IDs?"], "dec", w=160)
f2.node("ruse", 80, 125, ["r1 = ref IDs", "stereo from IDs"], w=150)
f2.node("rcan", L + 110, 125, ["canonical ranks for ref", "+ stereo flag"], w=200)
f2.node("mid", L, 215, ["mol has IDs?"], "dec", w=160)
f2.node("muse", 80, 300, ["r2 = mol IDs"], w=150)
f2.node("same", L + 110, 300, ["same atom types and", "same bonds as ref?", "`same_graph"], "dec", w=210, new=True)
f2.node("copy", L + 10, 400, ["r2 = r1", "(no 2nd canonicalization)"], w=200, new=True)
f2.node("mcan", L + 250, 400, ["canonical ranks for mol"], w=196)
f2.node("chk", L, 495, ["`checkranks(r1, r2)", "same atom count per rank?"], "dec", w=240)
f2.node("prov", L + 260, 590, ["IDs were provided?"], "dec", w=180)
f2.node("rec", L + 260, 670, ["recompute the provided", "side(s), check again"], w=200)
f2.node("fb", L + 260, 760, ["`fallbackranks", "rank = element"], w=170)
f2.node("ok", 80, 590, ["use r1, r2"], "term", w=130)
f2.edge("rid", "l", "ruse", "t", via=[(80, 40)], label="yes", lpos=(110, 33), anchor="middle")
f2.edge("rid", "r", "rcan", "t", via=[(L + 110, 40)], label="no", lpos=(L + 100, 33), anchor="middle")
f2.edge("ruse", "b", "mid", "t", via=[(80, 175), (L, 175)])
f2.edge("rcan", "b", "mid", "t", via=[(L + 110, 175), (L, 175)])
f2.edge("mid", "l", "muse", "t", via=[(80, 215)], label="yes", lpos=(110, 208), anchor="middle")
f2.edge("mid", "r", "same", "t", via=[(L + 110, 215)], label="no", lpos=(L + 100, 208), anchor="middle")
f2.edge("same", "b", "copy", "t", via=[(L + 110, 360), (L + 10, 360)], label="yes", lpos=(L + 40, 353), anchor="middle")
f2.edge("same", "r", "mcan", "t", via=[(L + 250, 300)], label="no", lpos=(L + 245, 293), anchor="middle")
f2.edge("muse", "b", "chk", "l", via=[(80, 495), (95, 495)])
f2.edge("copy", "b", "chk", "t", via=[(L + 10, 455), (L, 455)])
f2.edge("mcan", "b", "chk", "t", via=[(L + 250, 455), (L, 455)])
f2.edge("chk", "b", "ok", "t", via=[(L, 545), (80, 545)], label="yes", lpos=(150, 538), anchor="middle")
f2.edge("chk", "r", "prov", "t", via=[(L + 260, 495)], label="no", lpos=(L + 180, 488), anchor="middle")
f2.edge("prov", "b", "rec", "t", label="yes")
f2.edge("prov", "r", "fb", "r", via=[(640, 590), (640, 760)], label="no", lpos=(600, 583), anchor="middle")
f2.edge("rec", "b", "fb", "t", label="still fails")
f2.edge("rec", "l", "ok", "b", via=[(80, 670)], label="passes", lpos=(230, 663), anchor="middle")
# canonical ranks box
f2.frame(560, 90, 325, 470, "Canonical ranks (canonical.f90)")
f2.node("c1", 722, 150, ["CN bonds -> neighbour lists", "H atoms: isolated nodes"], w=280)
f2.node("c2", 722, 230, ["BFS from every heavy atom", "shell sums rinv over distances", "other fragments = last shell"], w=280, new=True)
f2.node("c3", 722, 320, ["initial invariant: rank of rinv,", "element, number of H"], w=280)
f2.node("c4", 722, 405, ["1-WL refinement: split ranks by", "product of neighbour primes", "(sort; overflow-safe)"], w=280, new=True)
f2.node("c5", 722, 500, ["heavy ranks renumbered,", "H take their parent's class"], w=280)
f2.edge("c1", "b", "c2", "t")
f2.edge("c2", "b", "c3", "t")
f2.edge("c3", "b", "c4", "t")
f2.edge("c4", "b", "c5", "t")
f2.edge("rcan", "r", "c1", "l", via=[(500, 125), (500, 150)], dashed=True)

# =====================================================================
# Figure 3: min_rmsd
# =====================================================================
f3 = Fig("min", 960, 1180, "Flowchart and decision tree of min_rmsd")
M = 480
f3.node("s", M, 35, ["ref, mol, ranks, stereo flag"], "term", w=260)
f3.node("topo", M, 110, ["topology check: same size,", "same max rank, same sorted", "(element, rank) pairs?"], "dec", w=280)
f3.node("tfail", 820, 110, ["warning; plain RMSD", "without reordering", "(or huge), io > 0"], "term", w=200)
f3.node("ord", M, 205, ["mol ranks in ref order?"], "dec", w=230)
f3.node("sort", 790, 205, ["sort mol atoms by rank", "`rank_2_order + molatomsort"], w=250)
f3.node("grp", M, 285, ["count atoms per rank;", "ranks with 1 atom are assigned"], w=280)
f3.node("nu", M, 370, ["nunique >= 3 ?", "(atoms with a unique rank)"], "dec", w=250)
# left branch: unique atoms
A = 220
f3.frame(30, 420, 380, 560, "Unique-atom alignment")
f3.node("a1", A, 470, ["mol to centre of mass"], w=240)
f3.node("a2", A, 540, ["Kabsch on the unique atoms", "rotate mol"], w=260)
f3.node("a3", A, 620, ["LSAP for every rank with", ">1 atom (cost: squared", "distance) -> cost c1"], w=260)
f3.node("a4", A, 710, ["mirror allowed?"], "dec", w=170)
f3.node("a5", A, 800, ["mirror mol (z -> -z),", "Kabsch uniques, LSAP", "-> cost c2, unmirror"], w=260)
f3.node("a6", A, 905, ["keep lower of c1, c2", "(mirror mol again if c2)"], w=260)
f3.edge("a1", "b", "a2", "t")
f3.edge("a2", "b", "a3", "t")
f3.edge("a3", "b", "a4", "t")
f3.edge("a4", "b", "a5", "t", label="yes")
f3.edge("a4", "l", "a6", "l", via=[(60, 710), (60, 905)], label="no", lpos=(68, 760))
f3.edge("a5", "b", "a6", "t")
# right branch: inertia axes
B = 715
f3.frame(450, 420, 480, 560, "Principal-axis alignment")
f3.node("b1", B, 470, ["`axis(mol)", "centre of mass, principal axes,", "rotational constants"], w=290)
f3.node("b2", B, 562, ["distinct principal", "moments?"], "dec", w=210)
f3.node("b3", 545, 650, ["all 3", "4 orientations"], w=140)
f3.node("b4", B, 650, ["A or C only", "8 (with 90 deg tilt)"], w=160)
f3.node("b5", 870, 650, ["none", "16 (3 tilts)"], w=112)
f3.node("b6", B, 745, ["per orientation: 180 deg flips", "about x, y, z; LSAP per rank", "-> costs 1..16"], w=290)
f3.node("b7", B, 835, ["mirror allowed?"], "dec", w=170)
f3.node("b8", B, 915, ["mirror z, re-axis, repeat", "-> costs 17..32; take minimum"], w=290)
f3.node("fin", M, 1040, ["apply the chosen orientation, mirror", "and assignment (molatomsort)"], w=330)
f3.node("rm", M, 1130, ["Kabsch RMSD; align mol if asked", "-> refine_assignment (figure 4)"], "term", w=330)
f3.edge("b1", "b", "b2", "t")
f3.edge("b2", "l", "b3", "t", via=[(545, 562)])
f3.edge("b2", "b", "b4", "t")
f3.edge("b2", "r", "b5", "t", via=[(870, 562)])
f3.edge("b3", "b", "b6", "t", via=[(545, 700), (B, 700)])
f3.edge("b4", "b", "b6", "t")
f3.edge("b5", "b", "b6", "t", via=[(870, 700), (B, 700)])
f3.edge("b6", "b", "b7", "t")
f3.edge("b7", "b", "b8", "t", label="yes")
f3.edge("b7", "r", "fin", "r", via=[(905, 835), (905, 1040)], label="no", lpos=(912, 860))
# branch edges
f3.edge("s", "b", "topo", "t")
f3.edge("topo", "r", "tfail", "l", label="no")
f3.edge("topo", "b", "ord", "t", label="yes")
f3.edge("ord", "r", "sort", "l", label="no")
f3.edge("ord", "b", "grp", "t", label="yes")
f3.edge("sort", "b", "grp", "r", via=[(790, 285)])
f3.edge("grp", "b", "nu", "t")
f3.edge("nu", "l", "a1", "t", via=[(A, 370)], label="yes", lpos=(330, 363), anchor="middle")
f3.edge("nu", "r", "b1", "t", via=[(B, 370)], label="no", lpos=(640, 363), anchor="middle")
f3.edge("a6", "b", "fin", "l", via=[(A, 1040)])
f3.edge("b8", "b", "fin", "t", via=[(B, 990), (M, 990)])
f3.edge("fin", "b", "rm", "t")

# =====================================================================
# Figure 4: refine_assignment
# =====================================================================
f4 = Fig("refine", 900, 640, "Candidates and selection in refine_assignment")
f4.node("s", 450, 35, ["min_rmsd output (matched, aligned)"], "term", w=300)
f4.node("adj", 450, 105, ["CN bond matrices of ref and mol"], w=300, new=True)
CA, CB, CC = 125, 450, 750
f4.node("ca", CA, 200, ["candidate A", "min_rmsd's assignment"], w=200)
f4.node("clu", CB, 200, ["largest fragment >= half", "of the heavy atoms?"], "dec", w=230, new=True)
f4.node("cb", CB, 300, ["candidate B", "`walk_assignment(mol)", "figure 5"], w=210, new=True)
f4.node("mir", CC, 300, ["mirror allowed?"], "dec", w=170)
f4.node("cc", CC, 395, ["candidate C", "`walk_assignment", "(z-mirrored mol)"], w=210, new=True)
f4.node("term", 450, 490, ["each candidate: terminal-atom refinement", "(siblings on one parent, brute force <= 8,",
                           "4 starts, alternating with Kabsch)", "then bond check on the result"], w=380, new=True)
f4.node("pick", 450, 590, ["lowest RMSD among bond-preserving", "candidates; if none: candidate A"], "term", w=330, new=True)
f4.edge("s", "b", "adj", "t")
f4.edge("adj", "b", "ca", "t", via=[(450, 150), (CA, 150)])
f4.edge("adj", "b", "clu", "t")
f4.edge("clu", "b", "cb", "t", label="yes")
f4.edge("clu", "l", "ca", "r", label="no: A only")
f4.edge("cb", "r", "mir", "l")
f4.edge("mir", "b", "cc", "t", label="yes")
f4.edge("ca", "b", "term", "l", via=[(CA, 490)])
f4.edge("cb", "b", "term", "t")
f4.edge("cc", "b", "term", "r", via=[(CC, 490)])
f4.edge("term", "b", "pick", "t")

# =====================================================================
# Figure 5: the centre walk
# =====================================================================
f5 = Fig("walk", 960, 1100, "Flowchart of walk_assignment: setup, depth-first walk, polish and local search")
W = 360
RAIL = 40
f5.node("s1", W, 40, ["neighbour lists; heavy atoms with a", "unique rank are matched up front"], w=320, new=True)
f5.node("s2", W, 120, ["walk order: BFS from each fragment's", "centre (min eccentricity), largest first"], w=320, new=True)
f5.frame(20, 170, 920, 530, "Depth-first walk over the order")
f5.node("x", W, 215, ["next atom x in the order"], w=230, new=True)
f5.node("anc", W, 295, ["x has a matched", "neighbour p?"], "dec", w=190, new=True)
f5.node("cand1", 170, 385, ["free same-rank neighbours", "of p's partner"], w=210, new=True)
f5.node("cand2", 550, 385, ["all free same-rank atoms", "(fragment start)"], w=200, new=True)
f5.node("vf2", W, 470, ["keep y if degree matches and x's matched", "neighbours map onto exactly y's"], w=340, new=True)
f5.node("k", W, 555, ["how many?"], "dec", w=150, new=True)
f5.node("bt", 140, 555, ["0: backtrack to", "the last choice with", "untried candidates"], w=170, new=True)
f5.node("one", W, 640, ["1: take it"], w=120, new=True)
f5.node("fr", 590, 555, [">1: matched atoms", "define a frame?"], "dec", w=170, new=True)
f5.node("near", 590, 640, ["Kabsch on matched;", "take nearest, keep rest"], w=190, new=True)
f5.node("br", 810, 640, ["try each in turn", "(up to 64 complete", "assignments)"], w=160, new=True)
f5.edge("s1", "b", "s2", "t")
f5.edge("s2", "b", "x", "t")
f5.edge("x", "b", "anc", "t")
f5.edge("anc", "l", "cand1", "t", via=[(170, 295)], label="yes", lpos=(215, 288), anchor="middle")
f5.edge("anc", "r", "cand2", "t", via=[(550, 295)], label="no", lpos=(505, 288), anchor="middle")
f5.edge("cand1", "b", "vf2", "t", via=[(170, 430), (W, 430)])
f5.edge("cand2", "b", "vf2", "t", via=[(550, 430), (W, 430)])
f5.edge("vf2", "b", "k", "t")
f5.edge("k", "l", "bt", "r")
f5.edge("k", "b", "one", "t")
f5.edge("k", "r", "fr", "l")
f5.edge("fr", "b", "near", "t", label="yes")
f5.edge("fr", "r", "br", "t", via=[(810, 555)], label="no", lpos=(740, 548), anchor="middle")
# loop-backs share one rail on the left
f5.edge("bt", "l", "x", "l", via=[(RAIL, 555), (RAIL, 215)], dashed=True)
f5.edge("one", "b", "x", "l", via=[(W, 680), (RAIL, 680), (RAIL, 215)], dashed=True)
f5.edge("near", "b", "x", "l", via=[(590, 680), (RAIL, 680), (RAIL, 215)], dashed=True)
f5.edge("br", "b", "x", "l", via=[(810, 680), (RAIL, 680), (RAIL, 215)], dashed=True)
f5.text(RAIL + 8, 695, "next position", "lbl")
f5.edge("x", "r", "leaf", "r", via=[(905, 215), (905, 760)], label="order complete", lpos=(560, 208), anchor="middle")
# polish
f5.frame(20, 708, 920, 372, "Polish (per start), then local search (best start only)")
f5.node("leaf", 360, 760, ["complete heavy-atom assignment"], "term", w=260, new=True)
f5.node("p1", 360, 835, ["H atoms to the nearest H", "on the parent's partner"], w=260, new=True)
f5.node("p2", 360, 910, ["Kabsch on all atoms", "-> frame, RMSD"], w=260, new=True)
f5.node("p3", 360, 990, ["re-walk, every choice", "nearest in that frame"], w=260, new=True)
f5.node("ch", 660, 990, ["assignment", "changed?"], "dec", w=150, new=True)
f5.node("ls", 660, 880, ["best start: try each other candidate", "at each bond-anchored choice;", "accept if total RMSD drops"], w=280, new=True)
f5.node("done", 885, 880, ["best", "assignment"], "term", w=100, new=True)
f5.edge("leaf", "b", "p1", "t")
f5.edge("p1", "b", "p2", "t")
f5.edge("p2", "b", "p3", "t")
f5.edge("p3", "r", "ch", "l")
f5.edge("ch", "b", "p1", "l", via=[(660, 1050), (190, 1050), (190, 835)], label="yes", lpos=(670, 1035))
f5.edge("ch", "t", "ls", "b", label="no")
f5.edge("ls", "t", "p1", "r", via=[(660, 800), (540, 800), (540, 835)], label="improved", lpos=(600, 792), anchor="middle")
f5.edge("ls", "r", "done", "l")
f5.text(885, 849, "no gain", "lbl", "middle")


# ---------------------------------------------------------------------
# links between the figures and to the worked examples
# ---------------------------------------------------------------------
for k in ("ids", "rr", "same", "cp", "rm", "cons", "fb"):
    f0.link(k, "ranks")
f0.link("st", "pipeline")
f0.link("topo", "min-rmsd")
for k in ("uq", "ua", "pa"):
    f0.link(k, "min-rmsd")
for k in ("ca", "clu", "conly", "val", "keepa", "low"):
    f0.link(k, "refine")
f0.link("walk", "walk")
f1.link("ranks", "ranks")
f1.link("min", "min-rmsd")
f1.link("topo", "min-rmsd")
f1.link("ref", "refine")
f3.link("a3", "ex-a")
f3.link("b6", "ex-a")
f3.link("rm", "refine")
f4.link("ca", "min-rmsd")
f4.link("cb", "walk")
f4.link("cc", "walk")
f4.link("pick", "ex-b")
f5.link("anc", "ex-b")
f5.link("k", "ex-b")
f5.link("ls", "ex-c")
f5.link("near", "ex-c")

# ---------------------------------------------------------------------
TEMPLATE = Path(__file__).with_name("template.html").read_text()
html = (
    TEMPLATE.replace("{{FIG0}}", f0.svg())
    .replace("{{FIG1}}", f1.svg())
    .replace("{{FIG2}}", f2.svg())
    .replace("{{FIG3}}", f3.svg())
    .replace("{{FIG4}}", f4.svg())
    .replace("{{FIG5}}", f5.svg())
    .replace("{{EXAMPLES}}", Path(__file__).with_name("examples.json").read_text().replace("</", "<\\/"))
)
out = Path(__file__).with_name("irmsd-flowchart.html")
out.write_text(html)
print("wrote", out)
