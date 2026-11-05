#!/usr/bin/env bash
# setup_md_swissparam.sh "<simulation folder>" ["<reference folder>"]
#
# Builds a solvated, neutralised protein-ligand system from a docked complex:
#   SwissParam ligand (unk.itp/unk.pdb) + receptor REC.pdb  ->  topology, box, water, ions,
#   ligand position restraints, index groups and EM.tpr; then test-runs grompp for NVT/NPT/MD.
#
# Required in the folder: REC.pdb, unk.itp, unk.pdb, ions.mdp, EM.mdp, NVT.mdp, NPT.mdp, MD.mdp
# Optional reference folder: force field and water model are copied from its topol.top,
#                            and its .mdp files are compared with yours.
# Settings (environment variables, defaults in brackets):
#   FF [charmm27]  WATER [tip3p]  CONC [0.1]  BOXD [1.0]  BOXTYPE [triclinic]  LIG_SRC [unk]
# Nothing that already exists is overwritten: the script stops if the folder is already set up.
set -uo pipefail; set +H; export LC_NUMERIC=C LC_ALL=C
die(){ echo; echo "!! STOP: $*"; exit 1; }
command -v gmx > /dev/null || die "gmx not found - activate your GROMACS environment first"
D="$(cd "${1:?usage: setup_md_swissparam.sh <folder> [reference folder]}" 2>/dev/null && pwd)" || die "folder not found: $1"
REFD="${2:-}"
FF="${FF:-charmm27}"; WAT="${WATER:-tip3p}"; CONC="${CONC:-0.1}"; BOXD="${BOXD:-1.0}"; BT="${BOXTYPE:-triclinic}"; SRC="${LIG_SRC:-unk}"
[ "$SRC" != "LIG" ] || die "LIG_SRC must not be LIG (LIG.itp/LIG.pdb are the files this script writes)"
cd "$D" || exit 1
mkdir -p setup_logs setup_check
exec > >(tee -a "$D/setup_logs/setup_log.txt") 2>&1
echo "=== setup $(date '+%F %T')  folder: $D"
g(){ local L="$1" IN="$2"; shift 2; echo "printf '$IN' | gmx $*" >> setup_logs/commands.txt
  gmx "$@" > "setup_logs/$L.log" 2>&1 <<< "$(printf '%b' "$IN")" || { echo "!! gmx step FAILED: $L"; tail -40 "setup_logs/$L.log"; exit 1; }
  if grep -E '^WARNING [0-9]' "setup_logs/$L.log" > /dev/null; then echo "  warnings in $L:"; grep -A5 -E '^WARNING [0-9]' "setup_logs/$L.log" | head -30; fi
  echo "  ok  $L"; }
for f in REC.pdb "$SRC.itp" "$SRC.pdb" ions.mdp EM.mdp NVT.mdp NPT.mdp MD.mdp; do [ -f "$f" ] || die "missing $f"; done
for f in conf.gro topol.top box.gro box_sol.gro box_sol_ion.gro EM.tpr LIG.itp LIG.pdb LIG.gro posre_LIG.itp index.ndx; do
  [ -e "$f" ] && die "$f already exists - this folder looks already set up; not overwriting anything"; done

if [ -n "$REFD" ] && [ -f "$REFD/topol.top" ]; then
  f=$(grep -oE '[A-Za-z0-9_.+-]+\.ff/forcefield\.itp' "$REFD/topol.top" | head -1 | sed 's/\.ff\/forcefield\.itp//')
  w=$(grep -oiE '(tip3p|tip4p|spce|spc|tips3p)\.itp' "$REFD/topol.top" | head -1 | sed 's/\.itp//')
  [ -n "$f" ] && FF="$f"; [ -n "$w" ] && WAT="$w"
  echo "  force field / water copied from reference $(basename "$REFD"): $FF / $WAT"
  for m in ions EM NVT NPT MD; do
    [ -f "$REFD/$m.mdp" ] && { cmp -s "$REFD/$m.mdp" "$m.mdp" && echo "  $m.mdp identical to reference" || echo "  NOTE: $m.mdp differs from reference $(basename "$REFD")"; }
  done
fi
echo "  settings: FF=$FF  WATER=$WAT  CONC=$CONC M  BOXD=$BOXD nm  BOXTYPE=$BT  ligand source=$SRC.itp/$SRC.pdb"
awk '/^ *(nsteps|dt) *=/{gsub(/;.*/,""); split($0,a,"="); gsub(/ /,"",a[1]); gsub(/ /,"",a[2]); v[a[1]]=a[2]} END{printf "  MD.mdp: nsteps=%s dt=%s -> %.1f ns\n", v["nsteps"], v["dt"], v["nsteps"]*v["dt"]/1000}' MD.mdp

echo "[1] ligand: $SRC.itp/$SRC.pdb -> LIG.itp/LIG.pdb (moleculetype and resname LIG)"
python3 - "$SRC" <<'PY'
import sys
src = sys.argv[1]
itp = open(f"{src}.itp", errors="ignore").read().splitlines()
out, sec, named, names, q = [], None, False, [], 0.0
for l in itp:
    s = l.strip()
    if s.startswith('['):
        sec = s.split(';')[0].strip('[] ').lower(); out.append(l); continue
    if sec == 'moleculetype' and s and not s.startswith(';') and not named:
        p = s.split(); out.append(f"LIG   {p[1] if len(p) > 1 else 3}"); named = True; continue
    if sec == 'atoms' and s and not s.startswith(';'):
        body, com = (l.split(';', 1) + [""])[:2]
        p = body.split()
        if len(p) >= 7 and p[0].isdigit():
            names.append(p[4]); p[3] = "LIG"; q += float(p[6])
            l = "  ".join(p) + ((" ;" + com) if com else "")
    out.append(l)
if not named or not names: sys.exit(f"!! could not find [ moleculetype ] / [ atoms ] in {src}.itp")
# position restraints belong to the ligand moleculetype, so they go at the end of LIG.itp
out += ["", "; Ligand position restraints", "#ifdef POSRES", '#include "posre_LIG.itp"', "#endif", ""]
open("LIG.itp", "w").write("\n".join(out) + "\n")
pdb = [l for l in open(f"{src}.pdb", errors="ignore") if l.startswith(("ATOM", "HETATM"))]
pn = [l[12:16].strip() for l in pdb]
if pn != names:
    sys.exit(f"!! atom names/order differ between {src}.pdb ({len(pn)} atoms) and {src}.itp ({len(names)} atoms)")
with open("LIG.pdb", "w") as f:
    for l in pdb: f.write(l[:17] + "LIG" + l[20:])
    f.write("END\n")
print(f"  ligand atoms: {len(names)}  net charge: {q:+.3f}")
if abs(q - round(q)) > 0.01: print("  WARNING: ligand net charge is not an integer - check the ligand topology")
PY
[ $? -eq 0 ] || die "ligand preparation failed (message above)"

echo "[2] protein: pdb2gmx ($FF / $WAT)"
python3 - <<'PY'
import sys
std = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP HSD HSE HSP CYX ASH GLH LYN ACE NME NH2".split())
res, chains = {}, []
for l in open("REC.pdb", errors="ignore"):
    if l.startswith(("ATOM", "HETATM")):
        r = l[17:20].strip(); res[r] = res.get(r, 0) + 1
        c = l[21]
        if not chains or chains[-1] != c: chains.append(c)
bad = sorted(r for r in res if r not in std)
print(f"  chains in REC.pdb: {chains}")
if bad:
    sys.exit("!! REC.pdb contains non-protein residues pdb2gmx cannot handle here: " + ", ".join(bad) + " - remove or parameterise them first")
PY
[ $? -eq 0 ] || die "REC.pdb check failed (message above)"
g pdb2gmx "" pdb2gmx -f REC.pdb -o conf.gro -p topol.top -i posre.itp -ff "$FF" -water "$WAT" -ignh
cp -p topol.top setup_logs/topol_after_pdb2gmx.top

echo "[3] ligand coordinates -> complex.gro; topol.top: ligand include + [ molecules ] LIG 1"
g lig_gro "" editconf -f LIG.pdb -o LIG.gro
python3 - <<'PY'
import sys
c = open("conf.gro").read().splitlines(); l = open("LIG.gro").read().splitlines()
pa = c[2:2 + int(c[1])]; la = l[2:2 + int(l[1])]
open("complex.gro", "w").write("\n".join([c[0], str(len(pa) + len(la))] + pa + la + [c[2 + int(c[1])]]) + "\n")
print(f"  protein atoms {len(pa)} + ligand atoms {len(la)} = {len(pa) + len(la)}")
t = open("topol.top").read().splitlines()
out, done = [], False
for s in t:
    out.append(s)
    if not done and s.strip().startswith("#include") and ".ff/forcefield.itp" in s:
        out += ["", "; Include ligand topology", '#include "LIG.itp"', ""]; done = True
if not done: sys.exit("!! forcefield include not found in topol.top")
while out and not out[-1].strip(): out.pop()
out.append("LIG                 1")
open("topol.top", "w").write("\n".join(out) + "\n")
mol = False
for s in out:
    if s.strip().lower().startswith("[ molecules"): mol = True; print("  [ molecules ]:"); continue
    if mol and s.strip() and not s.startswith(";"): print("    " + s)
PY
[ $? -eq 0 ] || die "merging failed (message above)"

echo "[4] box (-d $BOXD, $BT), solvate (spc216), ions ($CONC M, neutral)"
g box "" editconf -f complex.gro -d "$BOXD" -bt "$BT" -o box.gro
g solvate "" solvate -cp box.gro -cs spc216.gro -p topol.top -o box_sol.gro
g ions_grompp "" grompp -f ions.mdp -c box_sol.gro -p topol.top -o ION.tpr -maxwarn 2
g genion "SOL\n" genion -s ION.tpr -p topol.top -conc "$CONC" -neutral -o box_sol_ion.gro
tail -6 topol.top | sed 's/^/    /'

echo "[5] ligand position restraints (heavy atoms, 1000 kJ/mol/nm2) and index (Protein | LIG)"
g ndx_lig '0 & ! a H*\nq\n' make_ndx -f LIG.gro -o index_LIG.ndx
NG=$(grep -c '^\[' index_LIG.ndx); GRP=$(( NG - 1 ))
g genrestr "$GRP\n" genrestr -f LIG.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000
g ndx_sys '"Protein" | "LIG"\nq\n' make_ndx -f box_sol_ion.gro -o index.ndx
python3 - <<'PY'
import sys, re
groups = [l.strip("[] \n") for l in open("index.ndx") if l.startswith("[")]
print("  index groups:", ", ".join(groups))
miss = []
for m in ["NVT.mdp", "NPT.mdp", "MD.mdp"]:
    for l in open(m, errors="ignore"):
        s = l.split(";")[0]
        if re.match(r"\s*tc[-_]grps\s*=", s, re.I):
            for gname in s.split("=", 1)[1].split():
                if gname not in groups: miss.append(f"{m}: {gname}")
if miss: sys.exit("!! tc-grps not found in index.ndx: " + "; ".join(miss))
print("  tc-grps of NVT/NPT/MD all present in index.ndx")
PY
[ $? -eq 0 ] || die "index check failed (message above)"

echo "[6] EM.tpr and grompp test of NVT / NPT / MD (test files go to setup_check/)"
g em_grompp "" grompp -f EM.mdp -c box_sol_ion.gro -p topol.top -o EM.tpr -maxwarn 2
for m in NVT NPT MD; do
  g "test_$m" "" grompp -f $m.mdp -c box_sol_ion.gro -r box_sol_ion.gro -p topol.top -n index.ndx -o "setup_check/${m}_test.tpr" -po "setup_check/mdout_${m}.mdp" -maxwarn 2
done
echo "  total atoms: $(sed -n 2p box_sol_ion.gro | tr -d ' ')   box: $(tail -1 box_sol_ion.gro)"
echo "=== setup OK $(date '+%F %T') - EM.tpr ready; NVT/NPT/MD inputs validated"
echo "    next: bash $(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_md_protocol.sh \"$D\""
