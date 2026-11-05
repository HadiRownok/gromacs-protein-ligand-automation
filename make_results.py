#!/usr/bin/env python3
"""make_results.py <result folder>

Builds Excel workbooks, 600-dpi PNG figures, CSV copies and a QC report from the
GROMACS outputs written by run_md_analysis.sh. Every number comes from those outputs;
nothing is interpolated or invented. Exit code 0 = no QC FAIL, 2 = QC FAIL present.
"""
import sys, os, re, json, datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

OUT = os.path.abspath(sys.argv[1])
W, RAW, XL, FIG, CSV, REP, STR = [os.path.join(OUT, x) for x in
    ("work", "Raw_XVG", "Excel", "Figures_PNG", "Data_CSV", "Reports", "Structures")]
M = json.load(open(os.path.join(W, "meta.json")))
SIM, D, TNS = M["SIM"], M["D"], int(M["TNS"])
NPROT, NCA = int(M["NPROT"]), int(M["NCA"])
grid = np.arange(0.0, TNS * 1000.0 + 0.5, 250.0)
NG = len(grid); tns = grid / 1000.0
QC = []
def qc(level, check, detail=""):
    QC.append((level, check, str(detail))); print(f"  [{level}] {check}: {detail}")
norm = lambda s: re.sub(r"\s+", "-", s.strip())

def read_xvg(path):
    sets, cur, leg = [], [], []
    for line in open(path, errors="ignore"):
        if line.startswith("#"): continue
        if line.startswith("@"):
            m = re.match(r'@\s+s(\d+)\s+legend\s+"(.*)"', line)
            if m: leg.append(m.group(2))
            continue
        if line.startswith("&"):
            if cur: sets.append(np.array(cur)); cur = []
            continue
        p = line.split()
        if p:
            try: cur.append([float(x) for x in p])
            except ValueError: pass
    if cur: sets.append(np.array(cur))
    return sets, leg

S = {}
def load(key, fname, col=1, setno=0, strict=True, legend=None, required=True):
    p = os.path.join(RAW, fname)
    if not os.path.exists(p):
        qc("FAIL" if required else "WARN", "missing output", fname); return None
    sets, leg = read_xvg(p)
    if len(sets) <= setno or len(sets[setno]) == 0:
        qc("FAIL", "empty xvg", fname); return None
    a = sets[setno]
    if legend is not None:
        nl = [norm(l) for l in leg]
        if legend not in nl: return None
        col = nl.index(legend) + 1
    if a.ndim != 2 or a.shape[1] <= col:
        qc("FAIL", "column missing", f"{fname} column {col}"); return None
    t = a[:, 0]
    if strict:
        if len(t) != NG or np.max(np.abs(t - grid)) > 1e-3:
            qc("FAIL", "time grid", f"{fname}: {len(t)} rows, t={t[0]:.3f}..{t[-1]:.3f} ps (expected {NG} rows 0..{grid[-1]:.0f} ps every 250 ps)")
            return None
        v = a[:, col].astype(float)
    else:
        tr = np.round(t, 3); first = {}
        for i, x in enumerate(tr): first.setdefault(float(x), i)
        idx = [first.get(round(float(g), 3)) for g in grid]
        miss = sum(i is None for i in idx)
        if miss:
            qc("FAIL", "time grid", f"{fname}: {miss} of {NG} grid times (every 250 ps) not present"); return None
        uq, cn = np.unique(tr, return_counts=True)
        dup = set(float(x) for x in uq[cn > 1])
        nd = sum(1 for g in grid if round(float(g), 3) in dup)
        if nd:
            qc("FAIL", "duplicate time stamps", f"{fname}: {nd} grid times occur more than once"); return None
        v = a[idx, col].astype(float)
    S[key] = v
    return v

print(f"== building results for {SIM}: {NG} time points (0-{TNS} ns every 0.25 ns)")
load("rmsd_ca", "rmsd_calpha.xvg"); load("rmsd_bb", "rmsd_backbone.xvg"); load("rmsd_cplx", "rmsd_complex_heavy.xvg")
load("rmsd_ligmove", "rmsd_ligand_after_backbone_fit.xvg"); load("rmsd_ligconf", "rmsd_ligand_selffit.xvg")
load("rg_prot", "gyrate_protein.xvg"); load("rg_cplx", "gyrate_complex.xvg")
load("hb_pl", "hbnum_protein_ligand.xvg"); load("hb_intra", "hbnum_intra_protein.xvg")
load("hb_solute", "hbnum_solute.xvg"); load("hb_sltslv", "hbnum_solute_solvent.xvg")
load("sasa_prot", "sasa_protein.xvg"); load("sasa_cplx", "sasa_complex.xvg", 1)
load("sasa_cplx_prot", "sasa_complex.xvg", 2); load("sasa_lig_in", "sasa_complex.xvg", 3)
load("sasa_lig_free", "sasa_ligand_alone.xvg"); load("vdw_cplx", "vdwsurf_complex.xvg")
load("mindist", "mindist_ligand_protein.xvg"); load("ncont", "contacts_ligand_protein_0.4nm.xvg")
load("dssp_t", "dssp_num.xvg")

ELEG = []
pE = os.path.join(RAW, "energy_md.xvg")
if os.path.exists(pE):
    ELEG = [norm(l) for l in read_xvg(pE)[1]]
    for l in ELEG: load("E:" + l, "energy_md.xvg", strict=False, legend=l)
else: qc("FAIL", "missing output", "energy_md.xvg")
pI = os.path.join(RAW, "interaction_energy_protein_ligand.xvg")
if os.path.exists(pI):
    for l in [norm(x) for x in read_xvg(pI)[1]]: load("IE:" + l, "interaction_energy_protein_ligand.xvg", strict=False, legend=l)
coul_pl, lj_pl = S.get("IE:Coul-SR:Protein-Ligand"), S.get("IE:LJ-SR:Protein-Ligand")
eint = coul_pl + lj_pl if (coul_pl is not None and lj_pl is not None) else None
if eint is None: qc("FAIL", "interaction energy", "Coul-SR / LJ-SR Protein-Ligand not available")

# box lengths (vector lengths a, b, c)
comp, box_src = {}, ""
pB = os.path.join(RAW, "box.xvg")
if os.path.exists(pB):
    sets, leg = read_xvg(pB); nc = sets[0].shape[1] - 1
    names = [norm(l).replace("Box-", "") for l in leg] if len(leg) == nc else ["XX", "YY", "ZZ", "YX", "ZX", "ZY"][:nc]
    for i, nm in enumerate(names):
        v = load("BOX:" + nm, "box.xvg", col=i + 1)
        if v is not None: comp[nm if len(nm) == 2 else nm * 2] = v
    box_src = "gmx traj -ob (trajectory frames)"
else:
    for nm in ["X", "Y", "Z", "XX", "YY", "ZZ", "YX", "ZX", "ZY"]:
        v = S.get("E:Box-" + nm)
        if v is not None: comp[nm if len(nm) == 2 else nm * 2] = v
    box_src = "gmx energy Box terms (edr)"
z = np.zeros(NG)
if "XX" in comp:
    bx_a = comp["XX"]
    bx_b = np.sqrt(comp.get("YX", z) ** 2 + comp.get("YY", z) ** 2)
    bx_c = np.sqrt(comp.get("ZX", z) ** 2 + comp.get("ZY", z) ** 2 + comp.get("ZZ", z) ** 2)
    S["box_a"], S["box_b"], S["box_c"] = bx_a, bx_b, bx_c
else:
    bx_a = bx_b = bx_c = None; qc("FAIL", "box lengths", "no box data")

# secondary structure from dssp.dat
SS = {k: None for k in ["Helix", "Sheet", "Turn", "Coil", "Helix310", "HelixPi"]}
pD = os.path.join(W, "dssp.dat")
if os.path.exists(pD):
    L = [l.strip() for l in open(pD) if l.strip()]
    if len(L) != NG: qc("FAIL", "DSSP frames", f"{len(L)} lines in dssp.dat, expected {NG}")
    else:
        arr = {k: [] for k in SS}; nres = set()
        for s in L:
            r = s.replace("=", ""); n = len(r); nres.add(n)
            h, e, t = r.count("H"), r.count("E") + r.count("B"), r.count("T")
            g3, pi = r.count("G"), r.count("I")
            arr["Helix"].append(100 * h / n); arr["Sheet"].append(100 * e / n); arr["Turn"].append(100 * t / n)
            arr["Helix310"].append(100 * g3 / n); arr["HelixPi"].append(100 * pi / n)
            arr["Coil"].append(100 * (n - h - e - t - g3 - pi) / n)
        SS = {k: np.array(v) for k, v in arr.items()}
        qc("PASS" if nres == {NCA} else "WARN", "DSSP residues per frame", f"{sorted(nres)} (protein C-alpha: {NCA})")
else: qc("FAIL", "missing output", "dssp.dat")

# PCA
pc1 = pc2 = None
pP = os.path.join(RAW, "pca_projection_pc1_pc2.xvg")
if os.path.exists(pP):
    sets, _ = read_xvg(pP)
    if len(sets) >= 2 and all(len(s) == NG and np.max(np.abs(s[:, 0] - grid)) < 1e-3 for s in sets[:2]):
        pc1, pc2 = sets[0][:, 1], sets[1][:, 1]; S["pc1"], S["pc2"] = pc1, pc2
    else: qc("FAIL", "time grid", "pca_projection_pc1_pc2.xvg")
eig = None
pEv = os.path.join(RAW, "pca_eigenvalues.xvg")
if os.path.exists(pEv): eig = read_xvg(pEv)[0][0][:, 1]

# residues (consecutive 1..N) from the t=0 complex PDB
def parse_pdb(path):
    at = []
    for l in open(path, errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            occ = l[54:60].strip()
            at.append(dict(name=l[12:16].strip(), resn=l[17:21].strip(), chain=l[21].strip(),
                           resi=l[22:26].strip(), icode=l[26].strip(), occ=float(occ) if occ else 0.0))
    return at
def residue_list(atoms):
    res, prev = [], None
    for i, a in enumerate(atoms):
        k = (a["chain"], a["resi"], a["icode"])
        if k != prev: res.append(dict(chain=a["chain"], resi=a["resi"], resn=a["resn"], atoms=[])); prev = k
        res[-1]["atoms"].append(i)
    ch, last_chain, last_num = 1, None, None
    for r in res:
        try: num = int(r["resi"])
        except ValueError: num = None
        if last_chain is not None and ((r["chain"] and r["chain"] != last_chain) or (num is not None and last_num is not None and num <= last_num)):
            ch += 1
        r["chainidx"] = ch; last_chain, last_num = r["chain"], num
    return res
ref_atoms = parse_pdb(os.path.join(STR, "complex_0ns.pdb"))[:NPROT]
RES = residue_list(ref_atoms)
ca_res = []
for ri, r in enumerate(RES):
    for i in r["atoms"]:
        if ref_atoms[i]["name"] == "CA": ca_res.append(ri)
rmsf_df = None
pR = os.path.join(RAW, "rmsf_calpha.xvg")
if os.path.exists(pR):
    a = read_xvg(pR)[0][0]
    if len(a) != len(ca_res) or len(a) != NCA:
        qc("FAIL", "RMSF rows", f"{len(a)} rows vs {len(ca_res)} C-alpha in PDB / {NCA} in index")
    else:
        rmsf_df = pd.DataFrame({
            "Residue number": np.arange(1, len(a) + 1),
            "Original residue number": [RES[i]["resi"] for i in ca_res],
            "Residue name": [RES[i]["resn"] for i in ca_res],
            "Chain": [RES[i]["chain"] or f"Chain{RES[i]['chainidx']}" for i in ca_res],
            "Chain index": [RES[i]["chainidx"] for i in ca_res],
            "RMSF (nm)": a[:, 1], "RMSF (Å)": a[:, 1] * 10})
else: qc("FAIL", "missing output", "rmsf_calpha.xvg")

# per-residue ligand contact occupancy (gmx select -ofpdb writes the fraction of frames in the occupancy column)
occ_df = None
pO = os.path.join(W, "contact_occupancy.pdb")
if os.path.exists(pO):
    oat = parse_pdb(pO)[:NPROT]; ores = residue_list(oat)
    if len(ores) != len(RES): qc("WARN", "contact occupancy", "residue count differs from reference PDB - table skipped")
    else:
        occ_df = pd.DataFrame({
            "Residue number": np.arange(1, len(ores) + 1),
            "Original residue number": [r["resi"] for r in ores], "Residue name": [r["resn"] for r in ores],
            "Chain": [r["chain"] or f"Chain{r['chainidx']}" for r in ores],
            "Contact occupancy (% frames, any atom within 0.4 nm of ligand)": [100 * max(oat[i]["occ"] for i in r["atoms"]) for r in ores]})
        if rmsf_df is not None:
            occ_by_res = occ_df.iloc[:, 4].values
            rmsf_df["Ligand contact occupancy (%)"] = [occ_by_res[i] for i in ca_res]

# ---------------- QC ----------------
qc("PASS", "time stamps", f"all trajectory-derived series checked against the same {NG}-point 250 ps grid")
for k in ["rmsd_ca", "rmsd_bb", "rmsd_cplx", "rmsd_ligmove", "rmsd_ligconf"]:
    v = S.get(k)
    if v is not None: qc("PASS" if abs(v[0]) < 0.01 else "FAIL", f"RMSD(t=0) {k}", f"{v[0]:.4f} nm")
for k, lim in [("rg_prot", 0.3), ("rg_cplx", 0.3)]:
    v = S.get(k)
    if v is not None:
        j = np.max(np.abs(np.diff(v))); qc("PASS" if j < lim else "FAIL", f"frame-to-frame jump {k}", f"max {j:.4f} nm (limit {lim})")
for k in ["rmsd_cplx", "rmsd_ca", "rmsd_bb"]:
    v = S.get(k)
    if v is not None:
        j = np.max(np.abs(np.diff(v))); qc("PASS" if j < 0.5 else "FAIL", f"frame-to-frame jump {k} (chain split test)", f"max {j:.4f} nm")
for k in ["hb_pl", "hb_intra", "hb_solute", "hb_sltslv"]:
    v = S.get(k)
    if v is not None:
        ok = np.all(v >= 0) and np.all(np.abs(v - np.round(v)) < 1e-9)
        qc("PASS" if ok else "FAIL", f"H-bond counts {k}", "non-negative integers" if ok else "non-integer or negative values")
for k in ["sasa_prot", "sasa_cplx", "sasa_cplx_prot", "sasa_lig_in", "sasa_lig_free", "vdw_cplx"]:
    v = S.get(k)
    if v is not None: qc("PASS" if np.all(v >= 0) else "FAIL", f"SASA non-negative {k}", f"min {v.min():.4f} nm2")
lm = S.get("rmsd_ligmove")
if lm is not None:
    if lm.max() > 1.0: qc("WARN", "ligand RMSD > 1 nm", f"max {lm.max():.3f} nm; first > 1 nm at {tns[np.argmax(lm > 1.0)]:.2f} ns (check mindist)")
    else: qc("PASS", "ligand RMSD (after backbone fit)", f"max {lm.max():.3f} nm")
    j = np.abs(np.diff(lm))
    if j.max() >= 0.5: qc("WARN", "large ligand jump", f"max {j.max():.3f} nm between frames at {tns[np.argmax(j) + 1]:.2f} ns")
md = S.get("mindist")
if md is not None:
    away = md > 0.6
    if away.any(): qc("WARN", "ligand left contact", f"min distance > 0.6 nm in {away.sum()} frames (first at {tns[np.argmax(away)]:.2f} ns) - real result, not an error")
    else: qc("PASS", "ligand stays in contact", f"max ligand-protein min distance {md.max():.3f} nm")
if rmsf_df is not None and rmsf_df["RMSF (nm)"].max() > 1.0:
    qc("WARN", "RMSF > 1 nm", f"max {rmsf_df['RMSF (nm)'].max():.3f} nm at residue {int(rmsf_df['Residue number'][rmsf_df['RMSF (nm)'].idxmax()])} (often a terminus)")

# ---------------- tables ----------------
def col(v, f=1.0):
    if isinstance(v, str): v = S.get(v)
    return np.asarray(v, float) * f if v is not None else np.full(NG, np.nan)
def esum(names):
    used = [n for n in names if "E:" + n in S]
    return (np.sum([S["E:" + n] for n in used], axis=0) if used else None), used
TotE, _ = esum(["Total-Energy"]); BondE, bondT = esum(["Bond"])
AngE, angT = esum(["Angle", "U-B"]); DihE, dihT = esum(["Proper-Dih.", "Ryckaert-Bell.", "CMAP-Dih.", "Fourier-Dih.", "Restricted-Dih."])
PlaE, plaT = esum(["Improper-Dih.", "Per.-Imp.-Dih."]); CouE, couT = esum(["Coulomb-(SR)", "Coulomb-14", "Coul.-recip."])
VdwE, vdwT = esum(["LJ-(SR)", "LJ-14", "Disper.-corr.", "Buck.ham-(SR)"])
for nm, v in [("TotalEnergy", TotE), ("Bond", BondE), ("Angle", AngE), ("Dihedral", DihE), ("Planarity", PlaE), ("Coulomb", CouE), ("VdW", VdwE)]:
    if v is None: qc("WARN", "energy column empty", f"{nm}: no matching GROMACS term in edr")

C27 = ["Time[ps]", "Time[ns]", "CellLengthX", "CellLengthY", "CellLengthZ", "TotalEnergy", "Bond", "Angle", "Dihedral",
       "Planarity", "Coulomb", "VdW", "SurfVdW", "SurfMol", "SurfAcc", "SoluteHBonds", "SltSlvHBonds", "Helix", "Sheet",
       "Turn", "Coil", "Helix310", "HelixPi", "RadGyration", "RMSDCa", "RMSDBb", "RMSDAll"]
d27 = {"Time[ps]": grid, "Time[ns]": tns, "CellLengthX": col(bx_a, 10), "CellLengthY": col(bx_b, 10), "CellLengthZ": col(bx_c, 10),
       "TotalEnergy": col(TotE), "Bond": col(BondE), "Angle": col(AngE), "Dihedral": col(DihE), "Planarity": col(PlaE),
       "Coulomb": col(CouE), "VdW": col(VdwE), "SurfVdW": col("vdw_cplx", 100), "SurfMol": np.full(NG, np.nan),
       "SurfAcc": col("sasa_cplx", 100), "SoluteHBonds": col("hb_solute"), "SltSlvHBonds": col("hb_sltslv"),
       "Helix": col(SS["Helix"]), "Sheet": col(SS["Sheet"]), "Turn": col(SS["Turn"]), "Coil": col(SS["Coil"]),
       "Helix310": col(SS["Helix310"]), "HelixPi": col(SS["HelixPi"]), "RadGyration": col("rg_cplx", 10),
       "RMSDCa": col("rmsd_ca", 10), "RMSDBb": col("rmsd_bb", 10), "RMSDAll": col("rmsd_cplx", 10)}
df27 = pd.DataFrame({c: d27[c] for c in C27})
buried = (col("sasa_prot") + col("sasa_lig_free") - col("sasa_cplx"))
ext = {"RMSDLigMove": col("rmsd_ligmove", 10), "RMSDLigConf": col("rmsd_ligconf", 10), "ProtLigHBonds": col("hb_pl"),
       "MinDistLigProt": col("mindist", 10), "LigContacts<4A": col("ncont"), "SurfAccLigInComplex": col("sasa_lig_in", 100),
       "SurfAccLigAlone": col("sasa_lig_free", 100), "BuriedSASA": buried * 100, "SurfAccProtein": col("sasa_prot", 100),
       "CoulSR_ProtLig": col(coul_pl), "LJSR_ProtLig": col(lj_pl), "Einter_ProtLig": col(eint),
       "RadGyrProtein": col("rg_prot", 10), "IntraProtHBonds": col("hb_intra"), "PotentialEnergy": col("E:Potential"),
       "Temperature": col("E:Temperature"), "Pressure": col("E:Pressure"), "Density": col("E:Density")}
df2 = pd.concat([df27, pd.DataFrame(ext)], axis=1)
DEF = [("Time[ps]", "ps", "Simulation time (0.25 ns grid)", "trajectory"), ("Time[ns]", "ns", "Simulation time", "trajectory"),
 ("CellLengthX/Y/Z", "Å", "Lengths of box vectors a, b, c", box_src),
 ("TotalEnergy", "kJ/mol", "GROMACS Total Energy (potential + kinetic), whole system", "gmx energy"),
 ("Bond", "kJ/mol", "+".join(bondT) or "n/a", "gmx energy"), ("Angle", "kJ/mol", "+".join(angT) or "n/a", "gmx energy"),
 ("Dihedral", "kJ/mol", "+".join(dihT) or "n/a", "gmx energy"), ("Planarity", "kJ/mol", "+".join(plaT) or "n/a", "gmx energy"),
 ("Coulomb", "kJ/mol", "+".join(couT) or "n/a", "gmx energy"), ("VdW", "kJ/mol", "+".join(vdwT) or "n/a", "gmx energy"),
 ("SurfVdW", "Å²", "van der Waals surface of protein-ligand complex (probe 0)", "gmx sasa -probe 0"),
 ("SurfMol", "Å²", "Molecular (SES) surface - NOT computable with GROMACS, left empty", "-"),
 ("SurfAcc", "Å²", "Solvent-accessible surface of protein-ligand complex (probe 0.14 nm)", "gmx sasa"),
 ("SoluteHBonds", "count", "H-bonds within the solute (protein + ligand)", "gmx hbond -r Complex -t Complex"),
 ("SltSlvHBonds", "count", "H-bonds between solute (protein + ligand) and water", "gmx hbond -r Complex -t Water (full system)"),
 ("Helix", "%", "alpha-helix (DSSP H)", "gmx dssp -hmode dssp"), ("Sheet", "%", "beta strand + bridge (E+B)", "gmx dssp"),
 ("Turn", "%", "turn (T)", "gmx dssp"), ("Coil", "%", "everything else: loop, bend, PPII (~, S, P)", "gmx dssp"),
 ("Helix310", "%", "3-10 helix (G)", "gmx dssp"), ("HelixPi", "%", "pi helix (I)", "gmx dssp"),
 ("RadGyration", "Å", "Radius of gyration of protein-ligand complex (mass weighted)", "gmx gyrate"),
 ("RMSDCa", "Å", "C-alpha RMSD vs t=0 after C-alpha fit", "gmx rms -nomw -nopbc"),
 ("RMSDBb", "Å", "Backbone RMSD vs t=0 after backbone fit", "gmx rms"),
 ("RMSDAll", "Å", "Heavy-atom RMSD of protein-ligand complex after heavy-atom fit", "gmx rms"),
 ("RMSDLigMove", "Å", "Ligand heavy-atom RMSD after fitting protein backbone (ligand movement in site)", "gmx rms"),
 ("RMSDLigConf", "Å", "Ligand heavy-atom RMSD after ligand self-fit (conformational change)", "gmx rms"),
 ("ProtLigHBonds", "count", "Protein-ligand H-bonds (D-A <= 0.35 nm, H-D-A <= 30 deg)", "gmx hbond"),
 ("MinDistLigProt", "Å", "Minimum ligand-protein atom distance", "gmx mindist"),
 ("LigContacts<4A", "count", "Ligand-protein atom pairs closer than 0.4 nm", "gmx mindist -on"),
 ("SurfAccLigInComplex", "Å²", "Ligand SASA while bound", "gmx sasa"), ("SurfAccLigAlone", "Å²", "Ligand SASA alone (same coordinates)", "gmx sasa"),
 ("BuriedSASA", "Å²", "SASA(protein) + SASA(ligand alone) - SASA(complex)", "derived from gmx sasa"),
 ("SurfAccProtein", "Å²", "Protein SASA alone", "gmx sasa"),
 ("CoulSR_ProtLig / LJSR_ProtLig / Einter_ProtLig", "kJ/mol", "Short-range protein-ligand Coulomb, LJ and their sum (not a binding free energy)", "mdrun -rerun with energygrps"),
 ("RadGyrProtein", "Å", "Protein Rg", "gmx gyrate"), ("IntraProtHBonds", "count", "H-bonds within protein", "gmx hbond"),
 ("PotentialEnergy / Temperature / Pressure / Density", "kJ/mol, K, bar, kg/m³", "Whole system", "gmx energy")]

def statrows(df, skip=("Time[ps]", "Time[ns]")):
    rows = []
    for lab, fn in [("_ Mean", np.nanmean), ("_ Min", np.nanmin), ("_ Max", np.nanmax)]:
        r = {c: None for c in df.columns}; r["Time[ps]"] = lab
        for c in df.columns:
            if c in skip: continue
            v = df[c].values.astype(float)
            if np.isfinite(v).any(): r[c] = float(fn(v))
        rows.append(r)
    return pd.DataFrame(rows, columns=df.columns)

HL = PatternFill("solid", fgColor="FFF2CC")
def cv(x):
    if x is None or isinstance(x, str): return x
    try:
        if pd.isna(x): return None
    except Exception: pass
    if isinstance(x, np.integer): return int(x)
    if isinstance(x, (np.floating, float)): return float(x)
    return x
def put_df(ws, df, highlight=(), fmt=None):
    ws.append(list(df.columns))
    for c in ws[1]:
        c.font = Font(name="Calibri", size=14, bold=True)
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for row in df.itertuples(index=False): ws.append([cv(x) for x in row])
    for j, cname in enumerate(df.columns, 1):
        L = get_column_letter(j); ws.column_dimensions[L].width = max(12, min(45, len(str(cname)) + 4))
        f = (fmt or {}).get(cname)
        for cell in ws[L][1:]:
            if f and isinstance(cell.value, (int, float)): cell.number_format = f
            if cname in highlight: cell.fill = HL
        if cname in highlight: ws[L][0].fill = HL
    ws.freeze_panes = "A2"
def fmt_for(cols):
    f = {}
    for c in cols:
        if c == "Time[ps]": f[c] = "0"
        elif c == "Time[ns]": f[c] = "0.00"
        elif "HBonds" in c or "Contacts" in c: f[c] = "0"
        else: f[c] = "0.000"
    return f
HCOLS = ["Time[ns]", "SurfAcc", "SoluteHBonds", "RadGyration", "RMSDAll"]
os.makedirs(XL, exist_ok=True)
def roundtrip(path, sheet, df):
    ws = load_workbook(path, data_only=True)[sheet]
    hdr = [c.value for c in ws[1]]
    if hdr != list(df.columns): qc("FAIL", "Excel round-trip", f"{os.path.basename(path)}: header mismatch"); return
    bad = 0
    for j, c in enumerate(df.columns, 1):
        src = df[c].values[:NG]
        for i in range(NG):
            g = ws.cell(row=i + 2, column=j).value; s = src[i]
            if s is None or (not isinstance(s, str) and pd.isna(s)):
                bad += g is not None
            elif g is None or abs(float(g) - float(s)) > 1e-12 * max(1.0, abs(float(s))): bad += 1
    qc("PASS" if bad == 0 else "FAIL", "Excel round-trip",
       f"{os.path.basename(path)}: {bad} mismatching cells (relative tolerance 1e-12; Excel stores 15-16 significant digits)")

f1 = os.path.join(XL, f"analysis_{SIM}.xlsx")
wb = Workbook(); ws = wb.active; ws.title = "analysis"
put_df(ws, pd.concat([df27.astype(object), statrows(df27)], ignore_index=True), HCOLS, fmt_for(C27)); wb.save(f1)
roundtrip(f1, "analysis", df27)
f2 = os.path.join(XL, f"Analysis_2_{SIM}.xlsx")
wb = Workbook(); ws = wb.active; ws.title = "analysis"
put_df(ws, pd.concat([df2.astype(object), statrows(df2)], ignore_index=True), HCOLS, fmt_for(df2.columns))
ws = wb.create_sheet("Column_definitions")
put_df(ws, pd.DataFrame(DEF, columns=["Column", "Unit", "Definition", "Source / GROMACS tool"])); wb.save(f2)
roundtrip(f2, "analysis", df2)

T = {"Time (ns)": tns, "Time (ps)": grid}
sheets = {}
sheets["RMSD"] = pd.DataFrame({**T, "RMSD C-alpha (nm)": col("rmsd_ca"), "RMSD backbone (nm)": col("rmsd_bb"),
    "RMSD complex heavy atoms (nm)": col("rmsd_cplx"), "RMSD ligand after backbone fit (nm)": col("rmsd_ligmove"),
    "RMSD ligand self-fit (nm)": col("rmsd_ligconf")})
if rmsf_df is not None: sheets["RMSF"] = rmsf_df
sheets["Radius_of_Gyration"] = pd.DataFrame({**T, "Rg protein (nm)": col("rg_prot"), "Rg protein-ligand complex (nm)": col("rg_cplx")})
sheets["Solute_H_Bonds"] = pd.DataFrame({**T, "Protein-ligand H-bonds": col("hb_pl"), "Intra-protein H-bonds": col("hb_intra"),
    "Solute (protein+ligand) internal H-bonds": col("hb_solute"), "Solute-solvent H-bonds": col("hb_sltslv")})
sheets["SASA"] = pd.DataFrame({**T, "SASA protein alone (nm2)": col("sasa_prot"), "SASA complex (nm2)": col("sasa_cplx"),
    "SASA protein in complex (nm2)": col("sasa_cplx_prot"), "SASA ligand in complex (nm2)": col("sasa_lig_in"),
    "SASA ligand alone (nm2)": col("sasa_lig_free"), "Buried SASA (nm2)": buried, "vdW surface complex, probe 0 (nm2)": col("vdw_cplx")})
sheets["Ligand_Binding"] = pd.DataFrame({**T, "Min distance ligand-protein (nm)": col("mindist"), "Atom contacts < 0.4 nm": col("ncont"),
    "Protein-ligand H-bonds": col("hb_pl"), "RMSD ligand after backbone fit (nm)": col("rmsd_ligmove"),
    "RMSD ligand self-fit (nm)": col("rmsd_ligconf"), "Coul-SR protein-ligand (kJ/mol)": col(coul_pl),
    "LJ-SR protein-ligand (kJ/mol)": col(lj_pl), "Interaction energy SR total (kJ/mol)": col(eint)})
sheets["Secondary_Structure"] = pd.DataFrame({**T, **{f"{k} (%)": col(v) for k, v in SS.items()}})
sheets["Energy"] = pd.DataFrame({**T, **{k[2:]: S[k] for k in S if k.startswith("E:")},
    "Box a (nm)": col(bx_a), "Box b (nm)": col(bx_b), "Box c (nm)": col(bx_c)})
if pc1 is not None: sheets["PCA"] = pd.DataFrame({**T, "PC1 projection (nm)": pc1, "PC2 projection (nm)": pc2})
if eig is not None:
    sheets["PCA_Eigenvalues"] = pd.DataFrame({"Eigenvector": np.arange(1, len(eig) + 1), "Eigenvalue (nm2)": eig,
        "Variance (%)": 100 * eig / eig.sum(), "Cumulative (%)": 100 * np.cumsum(eig) / eig.sum()})
if occ_df is not None: sheets["Contact_Occupancy"] = occ_df

# ---------------- methods parameters (read from files) ----------------
def findf(name):
    for f in os.listdir(D):
        if f.lower() == name.lower(): return os.path.join(D, f)
    return None
def mdp(name):
    p = findf(name); d = {}
    if p:
        for l in open(p, errors="ignore"):
            s = l.split(";")[0].strip()
            if "=" in s:
                k, v = [x.strip() for x in s.split("=", 1)]; d[k.lower().replace("_", "-")] = v
    return d
def length(d):
    try: return f"{int(d['nsteps']) * float(d.get('dt', '0.001')) / 1000:g} ns ({d['nsteps']} steps x {d.get('dt', '0.001 (default)')} ps)"
    except Exception: return "n/a"
EMd, NVTd, NPTd, MDd = mdp("EM.mdp"), mdp("NVT.mdp"), mdp("NPT.mdp"), mdp("MD.mdp")
top = open(M["TOP"], errors="ignore").read()
ff = re.findall(r'#include\s+"[^"]*?([^"/]+)\.ff/forcefield\.itp"', top)
wat = re.findall(r'#include\s+"([^"]*(?:tip3p|tip4p|tip5p|spce|spc|tips3p|opc)[^"]*\.itp)"', top, re.I)
mols = open(os.path.join(W, "molecules.txt")).read().split("\n")
def box_info(gro):
    if not gro or not os.path.exists(gro): return "n/a"
    v = [float(x) for x in open(gro).read().strip().splitlines()[-1].split()]
    if len(v) == 3: return f"{v} nm; rectangular" + (" (cubic)" if max(v) - min(v) < 1e-3 else "")
    a = v[0]; tol = 2e-3 * a
    v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = v
    if abs(v2x) < tol and abs(v2y - a) < tol and abs(v3x - a / 2) < tol and abs(v3y - a / 2) < tol and abs(v3z - a * np.sqrt(2) / 2) < tol: s = "rhombic dodecahedron (xy-square)"
    elif abs(v2x - a / 2) < tol and abs(v2y - a * np.sqrt(3) / 2) < tol and abs(v3x - a / 2) < tol and abs(v3y - a * np.sqrt(3) / 6) < tol and abs(v3z - a * np.sqrt(6) / 3) < tol: s = "rhombic dodecahedron (xy-hexagon)"
    else: s = "triclinic"
    return f"vectors {v} nm; {s}"
emlog = findf("EM.log"); emtxt = open(emlog, errors="ignore").read() if emlog else ""
m_conv = re.search(r"converged to Fmax\s*<\s*(\S+)\s+in\s+(\d+)\s+steps", emtxt)
m_fmax = re.search(r"Maximum force\s*=\s*(\S+)", emtxt); m_epot = re.search(r"Potential Energy\s*=\s*(\S+)", emtxt)
mdver = re.search(r"GROMACS version:\s*(\S+)", open(M["MDLOG"], errors="ignore").read()) if os.path.exists(M["MDLOG"]) else None
lhead = "".join(open(M["LITP"], errors="ignore").readlines()[:6]).strip().replace("\n", " | ")
P = [("GROMACS version (production run, md.log)", mdver.group(1) if mdver else "n/a", "MD.log"),
     ("GROMACS version (this analysis)", M["GMXVER"], "gmx --version"),
     ("Force field", ", ".join(ff) or "not found in topol.top", "topol.top #include"),
     ("Water model include", ", ".join(wat) or "not found", "topol.top"),
     ("Ligand topology tool", M["LTOOL"] or "not stated in itp header", os.path.basename(M["LITP"])),
     ("Ligand itp header", lhead, os.path.basename(M["LITP"])),
     ("Ligand", f"{M['LIGM']} x{M['NLIG']}, resname {M['LRES']}, {M['LNAT']} atoms", "topol.top / itp"),
     ("[ molecules ]", "; ".join(x for x in mols if x), "topol.top"),
     ("Initial box (box.gro)", box_info(findf("box.gro")), "box.gro"),
     ("Production start box", box_info(os.path.join(W, "t0_whole.gro")), "MD frame 0"),
     ("Box-solute distance (editconf -d)", "not recorded in any file - take from your notes", "-"),
     ("Salt concentration", "not recorded in files (genion -conc); ion counts listed in [ molecules ]", "-")]
for stage, d in [("EM", EMd), ("NVT", NVTd), ("NPT", NPTd), ("MD", MDd)]:
    if not d: P.append((f"{stage} mdp", "file not found", "-")); continue
    if stage != "EM": P.append((f"{stage} length", length(d), f"{stage}.mdp"))
    for k in ["integrator", "emtol", "emstep", "nsteps", "dt", "tcoupl", "tc-grps", "tau-t", "ref-t", "pcoupl", "pcoupltype",
              "tau-p", "ref-p", "compressibility", "gen-vel", "gen-temp", "constraints", "constraint-algorithm", "cutoff-scheme",
              "coulombtype", "rcoulomb", "vdwtype", "vdw-modifier", "rvdw", "rvdw-switch", "dispcorr", "fourierspacing",
              "pme-order", "nstxout-compressed", "define", "refcoord-scaling"]:
        if k in d: P.append((f"{stage}: {k}", d[k], f"{stage}.mdp"))
if m_conv: P.append(("EM convergence", f"Fmax < {m_conv.group(1)} in {m_conv.group(2)} steps", "EM.log"))
if m_fmax: P.append(("EM final maximum force", m_fmax.group(1), "EM.log"))
if m_epot: P.append(("EM final potential energy", m_epot.group(1), "EM.log"))
dfP = pd.DataFrame(P, columns=["Parameter", "Value", "Source file"])
with open(os.path.join(REP, "Methods_parameters.txt"), "w") as fh:
    fh.write(f"Methods parameters read from files - {SIM}\n\n")
    for a, b, c in P: fh.write(f"{a}: {b}   [{c}]\n")
notes = [("Trajectory processing", "t=0 frame -pbc whole -> protein+ligand clustered (-pbc cluster) -> complex atoms extracted every 250 ps -> -pbc nojump (reference = clustered frame) -> -pbc mol -center Protein -ur compact for ligand/complex analyses"),
         ("RMSD reference", "t=0 complex (Structures/complex_0ns.pdb); gmx rms -nomw -nopbc"),
         ("H-bond criteria", "gmx hbond defaults: donor-acceptor <= 0.35 nm, hydrogen-donor-acceptor angle <= 30 deg"),
         ("SASA", "gmx sasa probe 0.14 nm; van der Waals surface with probe 0"),
         ("Secondary structure", "gmx dssp -hmode dssp; Coil = loop + bend + PPII"),
         ("Interaction energy", "mdrun -rerun of full system (frames every 250 ps) with energygrps = Protein Ligand; short-range Coulomb + LJ only; not a binding free energy"),
         ("PCA / FEL", f"C-alpha covariance (gmx covar/anaeig); FEL = -kT ln P(PC1,PC2), 30x30 bins, {NG} frames"),
         ("SurfMol", "Molecular surface (SES) cannot be computed with GROMACS - column left empty on purpose"),
         ("Units", "Standard analysis files use Å, Å², kJ/mol, %; detail workbook uses nm, nm²")]
sheets["Methods_Notes"] = pd.concat([pd.DataFrame(notes, columns=["Parameter", "Value"]).assign(**{"Source file": "this pipeline"}), dfP], ignore_index=True)

# ---------------- summary statistics ----------------
METRICS = [("RMSD C-alpha", "nm", "rmsd_ca"), ("RMSD backbone", "nm", "rmsd_bb"), ("RMSD complex heavy atoms", "nm", "rmsd_cplx"),
    ("RMSD ligand after backbone fit", "nm", "rmsd_ligmove"), ("RMSD ligand self-fit", "nm", "rmsd_ligconf"),
    ("Rg protein", "nm", "rg_prot"), ("Rg complex", "nm", "rg_cplx"), ("H-bonds protein-ligand", "count", "hb_pl"),
    ("H-bonds intra-protein", "count", "hb_intra"), ("H-bonds solute internal", "count", "hb_solute"),
    ("H-bonds solute-solvent", "count", "hb_sltslv"), ("SASA protein", "nm2", "sasa_prot"), ("SASA complex", "nm2", "sasa_cplx"),
    ("SASA ligand in complex", "nm2", "sasa_lig_in"), ("SASA ligand alone", "nm2", "sasa_lig_free"), ("Buried SASA", "nm2", buried),
    ("vdW surface complex", "nm2", "vdw_cplx"), ("Min distance ligand-protein", "nm", "mindist"), ("Contacts < 0.4 nm", "count", "ncont"),
    ("Coul-SR protein-ligand", "kJ/mol", coul_pl), ("LJ-SR protein-ligand", "kJ/mol", lj_pl), ("Interaction energy SR", "kJ/mol", eint)] + \
    [(f"{k}", "%", SS[k]) for k in SS] + \
    [("Potential energy", "kJ/mol", "E:Potential"), ("Total energy", "kJ/mol", "E:Total-Energy"), ("Temperature", "K", "E:Temperature"),
     ("Pressure", "bar", "E:Pressure"), ("Density", "kg/m3", "E:Density"), ("Box a", "nm", bx_a), ("Box b", "nm", bx_b), ("Box c", "nm", bx_c),
     ("PC1 projection", "nm", pc1), ("PC2 projection", "nm", pc2)]
LASTW = min(50.0, tns[-1] / 2)
last = tns >= tns[-1] - LASTW
rows = []
for name, unit, src in METRICS:
    v = col(src)
    if not np.isfinite(v).any(): continue
    rows.append({"Metric": name, "Unit": unit, "N frames": int(np.isfinite(v).sum()), "Mean": np.nanmean(v), "SD": np.nanstd(v, ddof=1),
                 "Median": np.nanmedian(v), "Min": np.nanmin(v), "Time of min (ns)": tns[np.nanargmin(v)], "Max": np.nanmax(v),
                 "Time of max (ns)": tns[np.nanargmax(v)], f"Mean last {LASTW:g} ns": np.nanmean(v[last]), f"SD last {LASTW:g} ns": np.nanstd(v[last], ddof=1)})
dfS = pd.DataFrame(rows)
bind = []
if S.get("hb_pl") is not None: bind.append(("Frames with >=1 protein-ligand H-bond (%)", 100 * np.mean(S["hb_pl"] >= 1)))
if S.get("mindist") is not None:
    bind.append(("Frames with ligand within 0.4 nm of protein (%)", 100 * np.mean(S["mindist"] < 0.4)))
    bind.append(("Frames with ligand min distance > 0.6 nm (%)", 100 * np.mean(S["mindist"] > 0.6)))
    bind.append(("Maximum ligand-protein min distance (nm)", float(np.max(S["mindist"]))))
if occ_df is not None:
    top_occ = occ_df.sort_values(occ_df.columns[4], ascending=False).head(15)
    for _, r in top_occ.iterrows():
        if r.iloc[4] > 0: bind.append((f"Contact occupancy {r['Residue name']}{r['Original residue number']} chain {r['Chain']} (%)", r.iloc[4]))
dfB = pd.DataFrame(bind, columns=["Quantity", "Value"])
dfRS = None
if rmsf_df is not None:
    r = rmsf_df["RMSF (nm)"]
    dfRS = pd.DataFrame([{"Mean RMSF (nm)": r.mean(), "SD (nm)": r.std(ddof=1), "Median (nm)": r.median(),
        "Min (nm)": r.min(), "Residue of min": f"{rmsf_df['Residue name'][r.idxmin()]}{rmsf_df['Original residue number'][r.idxmin()]} (#{int(rmsf_df['Residue number'][r.idxmin()])})",
        "Max (nm)": r.max(), "Residue of max": f"{rmsf_df['Residue name'][r.idxmax()]}{rmsf_df['Original residue number'][r.idxmax()]} (#{int(rmsf_df['Residue number'][r.idxmax()])})"}])
wb = Workbook(); ws = wb.active; ws.title = "Time_series_stats"; put_df(ws, dfS, fmt={c: "0.000" for c in dfS.columns})
if dfRS is not None: put_df(wb.create_sheet("RMSF_stats"), dfRS)
if len(dfB): put_df(wb.create_sheet("Binding_summary"), dfB, fmt={"Value": "0.0"})
wb.save(os.path.join(XL, "Summary_Statistics.xlsx"))

# ---------------- figures (600 dpi) ----------------
plt.rcParams.update({"font.size": 11, "axes.linewidth": 1.0, "savefig.dpi": 600})
def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(FIG, name), dpi=600); plt.close(fig)
def ts(name, items, ylab, title, step=False):
    items = [(l, np.asarray(v, float)) for l, v in items if v is not None]
    if not items: qc("WARN", "figure skipped", name); return
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for l, v in items:
        (ax.step(tns, v, where="mid", lw=0.7, label=l) if step else ax.plot(tns, v, lw=0.8, label=l))
    ax.set(xlabel="Time (ns)", ylabel=ylab, title=f"{title} - {SIM}", xlim=(0, tns[-1]))
    if len(items) > 1: ax.legend(frameon=False, fontsize=9)
    save(fig, name)
ts("01_RMSD_protein_complex.png", [("C-alpha", S.get("rmsd_ca")), ("Backbone", S.get("rmsd_bb")), ("Complex heavy atoms", S.get("rmsd_cplx"))], "RMSD (nm)", "RMSD")
ts("02_RMSD_ligand.png", [("Ligand after backbone fit", S.get("rmsd_ligmove")), ("Ligand self-fit", S.get("rmsd_ligconf"))], "RMSD (nm)", "Ligand RMSD")
if rmsf_df is not None:
    fig, ax = plt.subplots(figsize=(9, 4.2))
    ax.plot(rmsf_df["Residue number"], rmsf_df["RMSF (nm)"], lw=0.9, color="tab:blue")
    starts = rmsf_df["Residue number"][rmsf_df["Chain index"].diff().fillna(0) != 0]
    for s in starts: ax.axvline(s - 0.5, color="grey", ls="--", lw=0.8)
    if "Ligand contact occupancy (%)" in rmsf_df:
        m = rmsf_df["Ligand contact occupancy (%)"] >= 30
        ax.scatter(rmsf_df["Residue number"][m], rmsf_df["RMSF (nm)"][m], s=12, color="red", zorder=3, label="ligand contact >= 30% of frames")
        if m.any(): ax.legend(frameon=False, fontsize=9)
    ax.set(xlabel="Residue number (consecutive)", ylabel="RMSF (nm)", title=f"C-alpha RMSF - {SIM}", xlim=(1, len(rmsf_df)))
    save(fig, "03_RMSF_Calpha.png")
ts("04_Radius_of_gyration.png", [("Protein", S.get("rg_prot")), ("Protein-ligand complex", S.get("rg_cplx"))], "Rg (nm)", "Radius of gyration")
ts("05_HBonds_protein_ligand.png", [("Protein-ligand", S.get("hb_pl"))], "Number of H-bonds", "Protein-ligand H-bonds", step=True)
ts("06_HBonds_intra_protein_solute.png", [("Intra-protein", S.get("hb_intra")), ("Solute internal", S.get("hb_solute"))], "Number of H-bonds", "Solute H-bonds")
ts("07_HBonds_solute_solvent.png", [("Solute-water", S.get("hb_sltslv"))], "Number of H-bonds", "Solute-solvent H-bonds")
ts("08_SASA_protein_complex.png", [("Protein", S.get("sasa_prot")), ("Complex", S.get("sasa_cplx"))], "SASA (nm²)", "SASA")
ts("09_SASA_ligand_buried.png", [("Ligand in complex", S.get("sasa_lig_in")), ("Ligand alone", S.get("sasa_lig_free")), ("Buried SASA", buried if np.isfinite(buried).any() else None)], "Area (nm²)", "Ligand SASA / buried SASA")
if S.get("mindist") is not None:
    fig, ax = plt.subplots(figsize=(7.5, 4.2)); ax.plot(tns, S["mindist"], lw=0.8, color="tab:blue")
    ax.set(xlabel="Time (ns)", ylabel="Min. distance (nm)", xlim=(0, tns[-1]), title=f"Ligand-protein distance and contacts - {SIM}")
    if S.get("ncont") is not None:
        a2 = ax.twinx(); a2.plot(tns, S["ncont"], lw=0.6, color="tab:orange"); a2.set_ylabel("Contacts < 0.4 nm", color="tab:orange")
    save(fig, "10_MinDistance_contacts.png")
ts("11_Secondary_structure.png", [(k, SS[k]) for k in SS], "Residues (%)", "Secondary structure (DSSP)")
ts("12_Interaction_energy.png", [("Coulomb SR", coul_pl), ("LJ SR", lj_pl), ("Total SR", eint)], "Energy (kJ/mol)", "Protein-ligand interaction energy")
ts("13_Energy_total_potential.png", [("Total", S.get("E:Total-Energy")), ("Potential", S.get("E:Potential"))], "Energy (kJ/mol)", "System energy")
tpd = [(n, S.get("E:" + n), u) for n, u in [("Temperature", "K"), ("Pressure", "bar"), ("Density", "kg/m³")] if S.get("E:" + n) is not None]
if tpd:
    fig, axs = plt.subplots(len(tpd), 1, figsize=(7.5, 2.4 * len(tpd)), sharex=True); axs = np.atleast_1d(axs)
    for ax, (n, v, u) in zip(axs, tpd): ax.plot(tns, v, lw=0.6); ax.set_ylabel(f"{n} ({u})")
    axs[-1].set_xlabel("Time (ns)"); axs[0].set_title(f"Thermodynamic properties - {SIM}"); save(fig, "14_Temperature_Pressure_Density.png")
ts("15_Box_lengths.png", [("a", bx_a), ("b", bx_b), ("c", bx_c)], "Box vector length (nm)", "Box")
if eig is not None:
    n = min(20, len(eig)); fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(np.arange(1, n + 1), 100 * eig[:n] / eig.sum()); ax.plot(np.arange(1, n + 1), 100 * np.cumsum(eig[:n]) / eig.sum(), "o-", color="k", ms=3)
    ax.set(xlabel="Eigenvector", ylabel="Variance (%)", title=f"PCA (C-alpha) - {SIM}"); save(fig, "16_PCA_eigenvalues.png")
if pc1 is not None:
    fig, ax = plt.subplots(figsize=(6, 5)); sc = ax.scatter(pc1, pc2, c=tns, s=6, cmap="viridis")
    fig.colorbar(sc, label="Time (ns)"); ax.set(xlabel="PC1 (nm)", ylabel="PC2 (nm)", title=f"PC1 vs PC2 - {SIM}"); save(fig, "17_PCA_PC1_PC2.png")
    refT = None
    try: refT = float(MDd.get("ref-t", "").split()[0])
    except Exception: pass
    H, xe, ye = np.histogram2d(pc1, pc2, bins=30); Pp = H / H.sum()
    with np.errstate(divide="ignore"): G = -np.log(Pp / Pp.max())
    unit = "kT"
    if refT: G = G * 0.0083144626 * refT; unit = f"kJ/mol (T = {refT:g} K)"
    G = np.ma.masked_invalid(np.where(np.isfinite(G), G, np.nan))
    fig, ax = plt.subplots(figsize=(6, 5)); cf = ax.contourf(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), G.T, levels=20, cmap="jet")
    fig.colorbar(cf, label=f"Gibbs free energy ({unit})"); ax.set(xlabel="PC1 (nm)", ylabel="PC2 (nm)", title=f"Free energy landscape - {SIM}")
    save(fig, "18_Free_energy_landscape_PC1_PC2.png")
if occ_df is not None:
    oc = occ_df.columns[4]; sel = occ_df[occ_df[oc] >= 10]
    if len(sel):
        fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(sel) + 2), 4.2))
        multi = occ_df["Chain"].nunique() > 1
        ax.bar([f"{a}{b}" + (f":{c}" if multi else "") for a, b, c in zip(sel["Residue name"], sel["Original residue number"], sel["Chain"])], sel[oc], color="tab:green")
        ax.set(ylabel="Contact occupancy (%)", title=f"Residues within 0.4 nm of ligand (>=10% frames) - {SIM}", ylim=(0, 100))
        plt.xticks(rotation=90, fontsize=8); save(fig, "19_Ligand_contact_occupancy.png")

# ---------------- QC report, detail workbook, CSV, README ----------------
nF = sum(1 for q in QC if q[0] == "FAIL"); nW = sum(1 for q in QC if q[0] == "WARN")
sheets["QC"] = pd.DataFrame(QC, columns=["Level", "Check", "Detail"])
wb = Workbook(); wb.remove(wb.active)
for nm, df in sheets.items(): put_df(wb.create_sheet(nm[:31]), df)
wb.save(os.path.join(XL, "Protein_Ligand_MD_Analysis.xlsx"))
for nm, df in sheets.items(): df.to_csv(os.path.join(CSV, f"{nm}.csv"), index=False)
df2.to_csv(os.path.join(CSV, f"Analysis_2_{SIM}.csv"), index=False); dfS.to_csv(os.path.join(CSV, "Summary_Statistics.csv"), index=False)
with open(os.path.join(REP, "QC_REPORT.txt"), "w") as fh:
    fh.write(f"QC REPORT - {SIM}\nGenerated {datetime.datetime.now():%Y-%m-%d %H:%M:%S}\n")
    fh.write(f"Grid: {NG} points, 0-{TNS} ns every 0.25 ns\nRESULT: {nF} FAIL, {nW} WARN, {len(QC) - nF - nW} PASS\n\n")
    for l, c, d in QC: fh.write(f"[{l}] {c}: {d}\n")
with open(os.path.join(OUT, "README.txt"), "w") as fh:
    fh.write(f"""MD analysis results - {SIM}
Source: {D}  (original files untouched)
Production: {TNS} ns analysed on a 0.25 ns grid ({NG} frames). QC: {nF} FAIL, {nW} WARN (Reports/QC_REPORT.txt)

Excel/
  analysis_{SIM}.xlsx             27-column standard layout (Å, Å², kJ/mol, %); SurfMol empty (not computable in GROMACS)
  Analysis_2_{SIM}.xlsx           same + ligand columns, interaction energy, extras + Column_definitions
  Protein_Ligand_MD_Analysis.xlsx detail sheets in nm (RMSD, RMSF, Rg, H-bonds, SASA, ligand binding, DSSP, energy, PCA, contacts, methods, QC)
  Summary_Statistics.xlsx         mean/SD/median/min/max (+time), last-window stats, RMSF stats, binding summary
Figures_PNG/   600 dpi figures         Data_CSV/  CSV copies of all tables
Raw_XVG/       raw GROMACS outputs     Structures/ complex at 0 ns and {TNS} ns, input PDBs
Reports/       QC_REPORT.txt, gmx_commands.txt (exact commands), Methods_parameters.txt, pipeline_log.txt
work/          intermediate files (processed trajectories, index files, helper tpr)
""")
print(f"== {SIM}: {nF} FAIL, {nW} WARN")
sys.exit(2 if nF else 0)
