#!/usr/bin/env bash
# run_md_analysis.sh "<simulation folder>" <length_ns>
#
# Automated analysis of a GROMACS protein-ligand production run.
# The simulation folder is only read; all output goes to <folder>/MD_Analysis_<N>ns/.
# An existing result folder is renamed (..._old_HHMMSS), never overwritten.
set -uo pipefail
set +H
export LC_NUMERIC=C LC_ALL=C
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
die(){ echo; echo "!! STOP: $*"; exit 1; }
D_IN="${1:-}"; TNS="${2:-}"
[ -n "$D_IN" ] && [ -n "$TNS" ] || die "usage: run_md_analysis.sh <folder> <length_ns>"
[ -d "$D_IN" ] || die "folder not found: $D_IN"
command -v gmx > /dev/null || die "gmx not found - activate your GROMACS environment first"
D="$(cd "$D_IN" && pwd)"
SIM="$(basename "$D" | sed -E 's/[^A-Za-z0-9._-]+/_/g')"
TEND=$(( TNS * 1000 ))
OUT="$D/MD_Analysis_${TNS}ns"
if [ -e "$OUT" ]; then OLD="${OUT}_old_$(date +%H%M%S)"; mv "$OUT" "$OLD"; echo "existing result folder renamed -> $OLD"; fi
mkdir -p "$OUT"/{work,Raw_XVG,Excel,Figures_PNG,Data_CSV,Reports,Structures} || die "cannot create $OUT"
W="$OUT/work"; RAW="$OUT/Raw_XVG"; STR="$OUT/Structures"
exec > >(tee -a "$OUT/Reports/pipeline_log.txt") 2>&1
CMDS="$OUT/Reports/gmx_commands.txt"; : > "$CMDS"
cd "$W" || die "cannot cd $W"
echo "################################################################"
echo "SIMULATION: $SIM   ($(date '+%Y-%m-%d %H:%M:%S'))"
echo "gmx: $(command -v gmx)"
GMXVER=$(gmx --version 2>/dev/null | awk -F: '/GROMACS version/{gsub(/^ +/,"",$2); print $2; exit}')

# g <label> <stdin answers> <gmx subcommand and options...>: run, log, stop on failure
g(){ local L="$1" IN="$2"; shift 2
  { printf "printf '%s' | gmx" "$IN"; printf ' %q' "$@"; echo; } >> "$CMDS"
  gmx "$@" > "$W/log_$L.txt" 2>&1 <<< "$(printf '%b' "$IN")"
  local rc=$?
  [ $rc -eq 0 ] || { echo "!! gmx step FAILED: $L (exit $rc). Last 40 lines of $W/log_$L.txt:"; tail -40 "$W/log_$L.txt"; exit 1; }
  echo "  ok  $L"; }
# gp <label> <grompp options...>: helper tpr for analysis; warnings are shown, errors stop the run
gp(){ local L="$1"; shift
  { printf 'gmx grompp'; printf ' %q' "$@"; echo " -maxwarn 0"; } >> "$CMDS"
  gmx grompp "$@" -maxwarn 0 > "$W/log_$L.txt" 2>&1 && { echo "  ok  $L"; return 0; }
  if grep "Too many warnings" "$W/log_$L.txt" > /dev/null && ! grep -E '^ERROR [0-9]' "$W/log_$L.txt" > /dev/null; then
    echo "  grompp ($L) stopped only because of WARNINGS (shown below); helper tpr for analysis only -> repeating with -maxwarn 5"
    grep -A6 -E '^WARNING [0-9]' "$W/log_$L.txt"
    { printf 'gmx grompp'; printf ' %q' "$@"; echo " -maxwarn 5"; } >> "$CMDS"
    gmx grompp "$@" -maxwarn 5 > "$W/log_${L}_maxwarn.txt" 2>&1 && { echo "  ok  $L (warnings above)"; return 0; }
    L="${L}_maxwarn"
  fi
  echo "!! grompp FAILED ($L):"; tail -40 "$W/log_$L.txt"; exit 1; }
need(){ for f in "$@"; do [ -s "$f" ] || die "expected output missing or empty: $f"; done; }
cnt(){ awk -v g="$2" '/^\[/{gsub(/[\[\] ]/,""); f=($0==g); next} f{n+=NF} END{print n+0}' "$1"; }

echo "[1/9] detecting input files"
TPR=""; for t in "$D/MD.tpr" "$D/md.tpr"; do [ -f "$t" ] && { TPR="$t"; break; }; done
[ -n "$TPR" ] || { ls -1 "$D"; die "no MD.tpr / md.tpr found"; }
BASE="${TPR%.tpr}"; XTC="$BASE.xtc"; EDR="$BASE.edr"; MDLOG="$BASE.log"; MDMDP="$BASE.mdp"; TOP="$D/topol.top"
for f in "$XTC" "$TOP" "$EDR" "$MDMDP"; do [ -f "$f" ] || die "missing $f"; done
gmx check -f "$XTC" > "$W/check_xtc.txt" 2>&1 || { tail -20 "$W/check_xtc.txt"; die "gmx check failed"; }
LAST=$(tr '\r' '\n' < "$W/check_xtc.txt" | awk '/^Last frame/{print $5}' | tail -1)
[ -n "$LAST" ] || { tail -15 "$W/check_xtc.txt"; die "could not read the last frame time from gmx check"; }
awk -v l="$LAST" -v e="$TEND" 'BEGIN{exit !(l+0 >= e-0.001)}' || die "trajectory ends at ${LAST} ps, shorter than requested ${TEND} ps - not analysed"
echo "  xtc ends at $LAST ps (requested $TEND ps)"
if grep "Finished mdrun" "$MDLOG" > /dev/null 2>&1; then echo "  md.log: Finished mdrun found"; else echo "  WARNING: 'Finished mdrun' not found in $MDLOG"; fi
awk '/^\[ *molecules *\]/{f=1;next} /^\[/{f=0} f && $1!~/^;/ && NF>=2{print $1, $2}' "$TOP" > "$W/molecules.txt"
LIGM="${LIGAND:-}"
[ -n "$LIGM" ] || LIGM=$(awk '$1!~/^Protein/ && $1!~/^(SOL|HOH|WAT|TIP3P?|NA|CL|K|MG|CA|ZN|SOD|CLA|POT)$/ {print $1; exit}' "$W/molecules.txt")
[ -n "$LIGM" ] || { cat "$W/molecules.txt"; die "ligand molecule not detected - set LIGAND=<moleculetype name>"; }
NLIG=$(awk -v m="$LIGM" '$1==m{s+=$2} END{print s+0}' "$W/molecules.txt")
[ "$NLIG" -gt 0 ] || die "molecule $LIGM is not listed in [ molecules ]"
LITP=""
for f in "$D"/*.itp; do
  n=$(awk '/^\[ *moleculetype *\]/{f=1;next} f && $0!~/^ *;/ && NF>0 {print $1; exit}' "$f")
  if [ "$n" = "$LIGM" ]; then
    if grep -F "$(basename "$f")\"" "$TOP" > /dev/null; then LITP="$f"; break; fi
    [ -z "$LITP" ] && LITP="$f"
  fi
done
[ -n "$LITP" ] || { ls -1 "$D"/*.itp; die "itp defining $LIGM not found"; }
LNAT=$(awk '/^\[ *atoms/{a=1;next} /^\[/{a=0} a&&$1~/^[0-9]+$/{n++} END{print n+0}' "$LITP")
LRES=$(awk '/^\[ *atoms/{a=1;next} /^\[/{a=0} a&&$1~/^[0-9]+$/{print $4; exit}' "$LITP")
LTOOL=$(head -25 "$LITP" | grep -ioE 'SwissParam|cgenff|acpype|antechamber|LigParGen|ATB|MATCH' | head -1)
read -r LC LH LN LO LS LX < <(awk '/^\[ *atoms/{a=1;next} /^\[/{a=0} a&&$1~/^[0-9]+$/{m=$8+0; if(m>0.9&&m<1.2)h++; else if(m>11.5&&m<12.5)c++; else if(m>13.5&&m<14.5)n++; else if(m>15.5&&m<16.5)o++; else if(m>31&&m<33)s++; else x++} END{print c+0,h+0,n+0,o+0,s+0,x+0}' "$LITP")
LIG_NO=$(( LN + LO )); [ "$LX" -eq "$LNAT" ] && LIG_NO=-1
echo "  tpr=$TPR"; echo "  ligand molecule=$LIGM x$NLIG  resname=$LRES  itp=$(basename "$LITP") ($LNAT atoms, ${LTOOL:-tool not stated})"
echo "  ligand elements from itp masses: C=$LC H=$LH N=$LN O=$LO S=$LS other=$LX"
gmx hbond  -h > "$W/help_hbond.txt"  2>&1
gmx dssp   -h > "$W/help_dssp.txt"   2>&1
gmx gyrate -h > "$W/help_gyrate.txt" 2>&1
grep -E -- '-num' "$W/help_hbond.txt" > /dev/null || { tail -60 "$W/help_hbond.txt"; die "this gmx hbond has no -num option (needs GROMACS 2023 or newer)"; }
grep -E -- '-hmode' "$W/help_dssp.txt" > /dev/null || { tail -60 "$W/help_dssp.txt"; die "this gmx dssp has no -hmode option (needs GROMACS 2023 or newer)"; }
if grep -E -- '^ *-sel' "$W/help_gyrate.txt" > /dev/null; then GYRNEW=1; else GYRNEW=0; fi
if gmx traj -h > /dev/null 2>&1; then HAVETRAJ=1; else HAVETRAJ=0; fi

echo "[2/9] index groups"
SEL_C="\"Protein\" group \"Protein\"; \"Ligand\" resname $LRES; \"Complex\" group \"Protein\" or resname $LRES; \"Complex_Heavy\" (group \"Protein\" or resname $LRES) and mass > 1.5; \"Ligand_Heavy\" resname $LRES and mass > 1.5; \"Backbone\" group \"Backbone\"; \"C-alpha\" group \"C-alpha\"; \"Protein_Heavy\" group \"Protein\" and mass > 1.5; \"PH_ref\" group \"Protein-H\""
g index_full "" select -s "$TPR" -on "$W/idx_full.ndx" -select "$SEL_C; \"Water\" group \"Water\"; \"Water_and_ions\" group \"Water_and_ions\"; \"System\" all"
IF="$W/idx_full.ndx"
NPROT=$(cnt "$IF" Protein); NLIGAT=$(cnt "$IF" Ligand); NCX=$(cnt "$IF" Complex); NCA=$(cnt "$IF" C-alpha)
echo "  Protein=$NPROT  Ligand=$NLIGAT  Complex=$NCX  C-alpha=$NCA  Complex_Heavy=$(cnt "$IF" Complex_Heavy)  Ligand_Heavy=$(cnt "$IF" Ligand_Heavy)"
[ "$NLIGAT" -eq $(( LNAT * NLIG )) ] || die "Ligand group has $NLIGAT atoms but $(basename "$LITP") defines $LNAT x $NLIG"
[ "$NCX" -eq $(( NPROT + NLIGAT )) ] || die "Complex group size mismatch"
[ "$(cnt "$IF" Protein_Heavy)" -eq "$(cnt "$IF" PH_ref)" ] || die "heavy-atom selection (mass>1.5) differs from Protein-H"

echo "[3/9] analysis topologies (copies in work/topo; absolute includes; missing atomic numbers filled from mass)"
python3 - "$TOP" "$W" "$D" "$LIGM" <<'PY'
import sys, re, os
top, W, D, lig = sys.argv[1:]
TD = os.path.join(W, "topo"); os.makedirs(TD, exist_ok=True)
REF = [(1.008,1),(6.94,3),(10.81,5),(12.011,6),(14.007,7),(15.999,8),(18.998,9),(22.99,11),(24.305,12),
       (30.974,15),(32.06,16),(35.45,17),(39.098,19),(40.078,20),(65.38,30),(79.904,35),(126.904,53)]
def zmass(m):
    b = min(REF, key=lambda x: abs(x[0] - m))
    return b[1] if abs(b[0] - m) < 0.6 else None
def isint(x):
    try: int(x); return True
    except ValueError: return False
def fix_line(l):
    body = l.split(';')[0]; f = body.split()
    p = next((i for i, x in enumerate(f) if x in ("A", "D", "S", "V") and i >= 3), None)
    if p is None: return l, False, False
    try:
        if p == 3:                                 # name mass charge ptype ...
            z = zmass(float(f[1]))
            if not z: return l, False, True
            f = f[:1] + [str(z)] + f[1:]
        elif p == 4 and not isint(f[1]):           # name bond_type mass charge ptype ...
            z = zmass(float(f[2]))
            if not z: return l, False, True
            f = f[:2] + [str(z)] + f[2:]
        elif p == 4:                               # name at.num mass charge ptype ...
            if int(f[1]) > 0: return l, False, False
            z = zmass(float(f[2]))
            if not z: return l, False, True
            f[1] = str(z)
        elif p == 5:                               # name bond_type at.num mass charge ptype ...
            if isint(f[2]) and int(f[2]) > 0: return l, False, False
            z = zmass(float(f[3]))
            if not z: return l, False, True
            f[2] = str(z)
        else: return l, False, False
    except ValueError:
        return l, False, False
    return "  ".join(f) + "   ; atomic number from mass (analysis copy only)", True, False
done, report, samples = {}, [], []
def resolve(inc, base):
    if os.path.isabs(inc): return inc if os.path.exists(inc) else None
    for b in (base, D):
        p = os.path.normpath(os.path.join(b, inc))
        if os.path.exists(p): return p
    return None
def transform(path, molfilter=None):
    base = os.path.dirname(path); out, sec, nfix, nbad = [], None, 0, 0
    for l in open(path, errors="ignore").read().splitlines():
        s = l.strip()
        m = re.match(r'#include\s+"([^"]+)"', s)
        if m:
            r = resolve(m.group(1), base)
            if r: l = '#include "%s"' % process(r)
            out.append(l); continue
        if s.startswith('['):
            sec = s.split(';')[0].strip('[] ').lower()
        elif sec == 'atomtypes' and s and not s.startswith((';', '#')):
            if len(samples) < 6: samples.append(f"{os.path.basename(path)}: {s}")
            l2, fx, bad = fix_line(l)
            if fx: nfix += 1; l = l2
            if bad: nbad += 1
        elif molfilter is not None and sec == 'molecules' and s and not s.startswith(';'):
            if not molfilter(s.split()[0]): continue
        out.append(l)
    if nfix or nbad:
        report.append(f"{os.path.basename(path)}: atomic number filled for {nfix} atomtypes" + (f"; {nbad} could not be assigned" if nbad else ""))
    return out
def process(path):
    if path in done: return done[path]
    if ".ff/" in path: done[path] = path; return path
    done[path] = path
    out = transform(path)
    new = os.path.join(TD, f"{len(done)}_{os.path.basename(path)}")
    open(new, "w").write("\n".join(out) + "\n")
    done[path] = new
    return new
open(os.path.join(W, "complex.top"), "w").write("\n".join(transform(top, lambda n: n.lower().startswith("protein") or n == lig)) + "\n")
open(os.path.join(W, "full.top"), "w").write("\n".join(transform(top)) + "\n")
for r in sorted(set(report)): print("  " + r)
if not report: print("  no atomtypes needed an atomic number")
for s in samples: print("  atomtypes sample (non-force-field file) -> " + s)
PY
need "$W/complex.top" "$W/full.top"
g t0_whole "11\n" trjconv -s "$TPR" -f "$XTC" -n "$IF" -dump 0 -pbc whole -o "$W/t0_whole.gro"
g t0_cluster "2\n2\n" trjconv -s "$TPR" -f "$W/t0_whole.gro" -n "$IF" -pbc cluster -o "$W/t0_complex_cluster.gro"
need "$W/t0_whole.gro" "$W/t0_complex_cluster.gro"
cat > "$W/complex.mdp" <<'MDP'
integrator     = md
nsteps         = 0
dt             = 0.002
cutoff-scheme  = Verlet
coulombtype    = Reaction-Field
rcoulomb       = 1.2
rvdw           = 1.2
pbc            = xyz
tcoupl         = no
pcoupl         = no
constraints    = h-bonds
continuation   = yes
MDP
gp grompp_complex -f "$W/complex.mdp" -c "$W/t0_complex_cluster.gro" -p "$W/complex.top" -o "$W/complex.tpr" -po "$W/mdout_complex.mdp"
CT="$W/complex.tpr"
g index_complex "" select -s "$CT" -on "$W/idx_complex.ndx" -select "$SEL_C"
IC="$W/idx_complex.ndx"
[ "$(cnt "$IC" Complex)" -eq "$NCX" ] || die "complex tpr atom count differs from Complex group"
python3 - "$MDMDP" "$W/rerun.mdp" <<'PY'
import sys
src, dst = sys.argv[1:3]
keep = {'cutoff-scheme','nstlist','rlist','coulombtype','coulomb-modifier','rcoulomb','rcoulomb-switch','epsilon-r',
        'epsilon-rf','vdwtype','vdw-type','vdw-modifier','rvdw','rvdw-switch','dispcorr','fourierspacing','fourier-nx',
        'fourier-ny','fourier-nz','pme-order','ewald-rtol','ewald-geometry','constraints','constraint-algorithm',
        'lincs-order','lincs-iter','pbc','dt','verlet-buffer-tolerance'}
out = []
for l in open(src, errors="ignore"):
    s = l.split(';')[0].strip()
    if '=' not in s: continue
    k, v = [x.strip() for x in s.split('=', 1)]
    if k.lower().replace('_', '-') in keep and v: out.append(f"{k} = {v}")
out += ["integrator = md", "nsteps = 0", "tcoupl = no", "pcoupl = no", "gen-vel = no", "continuation = yes",
        "nstcalcenergy = 1", "nstenergy = 1", "energygrps = Protein Ligand"]
open(dst, 'w').write("\n".join(out) + "\n")
PY
gp grompp_full -f "$W/rerun.mdp" -c "$W/t0_whole.gro" -p "$W/full.top" -n "$IF" -o "$W/rerun.tpr" -po "$W/mdout_rerun.mdp"
FT="$W/rerun.tpr"

echo "[4/9] trajectory processing (cluster t=0 -> extract -dt 250 -> nojump -> center/compact)"
g extract "2\n" trjconv -s "$TPR" -f "$XTC" -n "$IF" -b 0 -e $TEND -dt 250 -o "$W/complex_raw.xtc"
g nojump "0\n" trjconv -s "$W/t0_complex_cluster.gro" -f "$W/complex_raw.xtc" -pbc nojump -o "$W/complex_nojump.xtc"
g center "0\n2\n" trjconv -s "$CT" -f "$W/complex_nojump.xtc" -n "$IC" -pbc mol -center -ur compact -o "$W/complex_center.xtc"
g full250 "11\n" trjconv -s "$TPR" -f "$XTC" -n "$IF" -b 0 -e $TEND -dt 250 -o "$W/full_250.xtc"
NJ="$W/complex_nojump.xtc"; CEN="$W/complex_center.xtc"
gmx check -f "$CEN" > "$W/check_center.txt" 2>&1
NFR=$(tr '\r' '\n' < "$W/check_center.txt" | awk '/^Coords/{print $2}' | tail -1)
[ "$NFR" = "$(( TEND / 250 + 1 ))" ] || die "processed trajectory has $NFR frames, expected $(( TEND / 250 + 1 ))"
echo "  processed frames: $NFR"
g pdb0 "2\n" trjconv -s "$CT" -f "$CEN" -n "$IC" -dump 0 -o "$STR/complex_0ns.pdb"
g pdbend "2\n" trjconv -s "$CT" -f "$CEN" -n "$IC" -dump $TEND -o "$STR/complex_${TNS}ns.pdb"
REF="$STR/complex_0ns.pdb"; need "$REF"
for f in "$D/$(basename "$D").pdb" "$D/REC.pdb" "$D/LIG.pdb"; do [ -f "$f" ] && cp -p "$f" "$STR/input_$(basename "$f" | sed -E 's/[^A-Za-z0-9._-]+/_/g')"; done

echo "[5/9] RMSD / RMSF / Rg / PCA"
g rmsd_ca       "6\n6\n" rms -s "$REF" -f "$NJ"  -n "$IC" -nomw -nopbc -o "$RAW/rmsd_calpha.xvg"
g rmsd_bb       "5\n5\n" rms -s "$REF" -f "$NJ"  -n "$IC" -nomw -nopbc -o "$RAW/rmsd_backbone.xvg"
g rmsd_complex  "3\n3\n" rms -s "$REF" -f "$CEN" -n "$IC" -nomw -nopbc -o "$RAW/rmsd_complex_heavy.xvg"
g rmsd_ligmove  "5\n4\n" rms -s "$REF" -f "$CEN" -n "$IC" -nomw -nopbc -o "$RAW/rmsd_ligand_after_backbone_fit.xvg"
g rmsd_ligconf  "4\n4\n" rms -s "$REF" -f "$CEN" -n "$IC" -nomw -nopbc -o "$RAW/rmsd_ligand_selffit.xvg"
g rmsf "6\n" rmsf -s "$CT" -f "$NJ" -n "$IC" -res -o "$RAW/rmsf_calpha.xvg"
if [ "$GYRNEW" = 1 ]; then
  g rg_protein "" gyrate -s "$CT" -f "$NJ"  -n "$IC" -sel 'group "Protein"' -o "$RAW/gyrate_protein.xvg"
  g rg_complex "" gyrate -s "$CT" -f "$CEN" -n "$IC" -sel 'group "Complex"' -o "$RAW/gyrate_complex.xvg"
else
  g rg_protein "0\n" gyrate -s "$CT" -f "$NJ"  -n "$IC" -o "$RAW/gyrate_protein.xvg"
  g rg_complex "2\n" gyrate -s "$CT" -f "$CEN" -n "$IC" -o "$RAW/gyrate_complex.xvg"
fi
g covar  "6\n6\n" covar  -s "$CT" -f "$NJ" -n "$IC" -o "$RAW/pca_eigenvalues.xvg" -v "$W/eigenvec.trr" -av "$W/pca_average.pdb" -l "$W/covar.log"
g anaeig "6\n6\n" anaeig -s "$CT" -f "$NJ" -n "$IC" -v "$W/eigenvec.trr" -first 1 -last 2 -proj "$RAW/pca_projection_pc1_pc2.xvg" -2d "$RAW/pca_2d_pc1_pc2.xvg"

echo "[6/9] hydrogen bonds (D-A <= 0.35 nm, angle <= 30 deg; gmx hbond defaults)"
if [ "$LIG_NO" -eq 0 ]; then
  echo "  ligand has no N or O atoms -> protein-ligand H-bonds are zero by definition (written as 0, reason in file header)"
  python3 - "$RAW/hbnum_protein_ligand.xvg" "$TEND" <<'PY'
import sys
out, tend = sys.argv[1], int(sys.argv[2])
with open(out, "w") as f:
    f.write("# Ligand itp contains no N or O atoms (element from atom masses): it has no H-bond donors or acceptors,\n")
    f.write("# so protein-ligand hydrogen bonds are 0 by definition in every frame (not computed by gmx hbond).\n")
    f.write('@    title "Protein-ligand hydrogen bonds"\n@    xaxis  label "Time (ps)"\n@    yaxis  label "Number"\n')
    for t in range(0, tend + 1, 250): f.write(f"{t:.3f} 0\n")
PY
else
  g hb_protlig "" hbond -s "$CT" -f "$CEN" -n "$IC" -r 'group "Protein"' -t 'group "Ligand"' -num "$RAW/hbnum_protein_ligand.xvg" -o "$W/hb_protlig.ndx"
  grep -E "has [0-9]+ acceptors" "$W/log_hb_protlig.txt" | sed 's/^/  /'
fi
g hb_intra   "" hbond -s "$CT" -f "$NJ"  -n "$IC" -r 'group "Protein"' -t 'group "Protein"' -num "$RAW/hbnum_intra_protein.xvg" -o "$W/hb_intra.ndx"
g hb_solute  "" hbond -s "$CT" -f "$CEN" -n "$IC" -r 'group "Complex"' -t 'group "Complex"' -num "$RAW/hbnum_solute.xvg" -o "$W/hb_solute.ndx"
g hb_sltslv  "" hbond -s "$FT" -f "$W/full_250.xtc" -n "$IF" -r 'group "Complex"' -t 'group "Water"' -num "$RAW/hbnum_solute_solvent.xvg" -o "$W/hb_sltslv.ndx"
grep -E "has [0-9]+ acceptors" "$W/log_hb_sltslv.txt" | sed 's/^/  /'

echo "[7/9] SASA, distances, contacts, secondary structure"
g sasa_protein "" sasa -s "$CT" -f "$NJ"  -n "$IC" -surface 'group "Protein"' -probe 0.14 -o "$RAW/sasa_protein.xvg"
g sasa_complex "" sasa -s "$CT" -f "$CEN" -n "$IC" -surface 'group "Complex"' -output 'group "Protein"; group "Ligand"' -probe 0.14 -o "$RAW/sasa_complex.xvg"
g sasa_ligand  "" sasa -s "$CT" -f "$CEN" -n "$IC" -surface 'group "Ligand"' -probe 0.14 -o "$RAW/sasa_ligand_alone.xvg"
g vdwsurf      "" sasa -s "$CT" -f "$CEN" -n "$IC" -surface 'group "Complex"' -probe 0 -o "$RAW/vdwsurf_complex.xvg"
g mindist "1\n0\n" mindist -s "$CT" -f "$CEN" -n "$IC" -d 0.4 -od "$RAW/mindist_ligand_protein.xvg" -on "$RAW/contacts_ligand_protein_0.4nm.xvg"
g occupancy "" select -s "$CT" -f "$CEN" -n "$IC" -select 'group "Protein" and same residue as within 0.4 of group "Ligand"' -ofpdb "$W/contact_occupancy.pdb"
g dssp "" dssp -s "$CT" -f "$NJ" -n "$IC" -sel 'group "Protein"' -hmode dssp -o "$W/dssp.dat" -num "$RAW/dssp_num.xvg"

echo "[8/9] energies, box, protein-ligand interaction energy (rerun with energy groups)"
gmx energy -f "$EDR" -o "$W/probe_energy.xvg" > "$W/energy_menu.txt" 2>&1 <<< ""
NUMS=$(python3 "$TOOLS/energy_pick.py" "$W/energy_menu.txt" md)
[ -n "$NUMS" ] || { tail -40 "$W/energy_menu.txt"; die "no energy terms recognised in $EDR"; }
g energy_md "$NUMS\n0\n" energy -f "$EDR" -o "$RAW/energy_md.xvg"
if [ "$HAVETRAJ" = 1 ]; then
  g box "1\n" traj -s "$TPR" -f "$XTC" -n "$IF" -b 0 -e $TEND -dt 250 -ob "$RAW/box.xvg"
else echo "  gmx traj not available -> box lengths will be taken from the edr Box terms"; fi
echo "gmx mdrun -s rerun.tpr -rerun full_250.xtc -deffnm rerun -nb cpu -pme cpu -bonded cpu" >> "$CMDS"
gmx mdrun -s "$FT" -rerun "$W/full_250.xtc" -deffnm "$W/rerun" -nb cpu -pme cpu -bonded cpu > "$W/log_rerun.txt" 2>&1 || { tail -40 "$W/log_rerun.txt"; die "mdrun -rerun failed"; }
echo "  ok  rerun"
gmx energy -f "$W/rerun.edr" -o "$W/probe_rerun.xvg" > "$W/energy_menu_rerun.txt" 2>&1 <<< ""
INUMS=$(python3 "$TOOLS/energy_pick.py" "$W/energy_menu_rerun.txt" inter)
[ -n "$INUMS" ] || { tail -40 "$W/energy_menu_rerun.txt"; die "Protein-Ligand energy-group terms not found in rerun.edr"; }
g energy_inter "$INUMS\n0\n" energy -f "$W/rerun.edr" -o "$RAW/interaction_energy_protein_ligand.xvg"

echo "[9/9] Excel, figures, QC"
export SIM D TNS TPR XTC EDR TOP MDLOG LIGM LRES LITP LNAT NLIG LTOOL NPROT NCA NLIGAT GMXVER GYRNEW HAVETRAJ
python3 - "$W/meta.json" <<'PY'
import os, json, sys
keys = "SIM D TNS TPR XTC EDR TOP MDLOG LIGM LRES LITP LNAT NLIG LTOOL NPROT NCA NLIGAT GMXVER GYRNEW HAVETRAJ".split()
json.dump({k: os.environ.get(k, "") for k in keys}, open(sys.argv[1], "w"), indent=1)
PY
python3 "$TOOLS/make_results.py" "$OUT"; PRC=$?
[ $PRC -eq 0 ] || [ $PRC -eq 2 ] || die "make_results.py crashed (exit $PRC) - see the Python error above"
rm -f "$W/full_250.xtc" "$W/probe_energy.xvg" "$W/probe_rerun.xvg"
echo; echo "RESULTS: $OUT"
ls -1 "$OUT/Excel"; echo "figures: $(ls -1 "$OUT/Figures_PNG" | wc -l) PNG"
[ $PRC -eq 0 ] && echo ">>> QC: 0 FAIL" || echo ">>> QC: FAIL items present - see $OUT/Reports/QC_REPORT.txt"
