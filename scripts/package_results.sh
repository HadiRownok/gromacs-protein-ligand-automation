#!/usr/bin/env bash
# package_results.sh "<root folder>" <length_ns> [output folder]
#
# Collects every <root>/<target>/<simulation>/MD_Analysis_<N>ns whose QC report says "0 FAIL"
# into one folder (results + complete analysis + simulation inputs) and zips it.
# The simulation folders are only read, never changed.
set -o pipefail; set +H; export LC_NUMERIC=C LC_ALL=C
ROOT="${1:?usage: package_results.sh <root folder> <length_ns> [output folder]}"; NS="${2:?length in ns}"
ROOT="$(cd "$ROOT" && pwd)" || exit 1
OUT="${3:-$(dirname "$ROOT")/MD_Results_${NS}ns_All}"
AN="MD_Analysis_${NS}ns"
[ -e "$OUT" ] && mv "$OUT" "${OUT}_old_$(date +%H%M%S)"
mkdir -p "$OUT"; n=0; skipped=0
for A in "$ROOT"/*/*/"$AN"; do
  [ -d "$A" ] || continue
  S="${A%/$AN}"; SIM=$(basename "$S"); T=$(basename "$(dirname "$S")")
  echo "######## $T / $SIM"
  if ! grep -qs "RESULT: 0 FAIL" "$A/Reports/QC_REPORT.txt"; then
    echo "   x QC not 0 FAIL (or missing) - skipped"; skipped=$((skipped + 1)); continue; fi
  n=$((n + 1)); D="$OUT/$(printf %02d $n)_${T}_$(echo "$SIM" | sed -E 's/[^A-Za-z0-9-]+/_/g; s/_+$//')"
  mkdir -p "$D/Protein_Ligand_Complex" "$D/Results_Excel" "$D/Figures" "$D/Simulation_Inputs"
  cp "$A"/Structures/*.pdb "$D/Protein_Ligand_Complex/"
  cp "$A"/Excel/*.xlsx "$D/Results_Excel/"
  cp "$A"/Figures_PNG/*.png "$D/Figures/"
  cp "$A/Reports/QC_REPORT.txt" "$A/Reports/Methods_parameters.txt" "$D/"
  for f in "$S"/*.mdp "$S"/*.top "$S"/*.itp "$S"/*.pdb "$S"/*.mol2 "$S"/*.gro "$S"/*.ndx "$S"/*.log \
           "$S"/*.tpr "$S"/*.edr "$S"/*.rtf "$S"/*.prm "$S"/*.par "$S"/*.psf "$S"/*.str; do
    [ -f "$f" ] && cp -p "$f" "$D/Simulation_Inputs/"; done
  cp -rp "$A" "$D/$AN"
  echo "   ok $(ls "$D/Results_Excel" | wc -l) Excel, $(ls "$D/Figures" | wc -l) figures | $(grep -h RESULT "$D/QC_REPORT.txt")"
done
cat > "$OUT/README.txt" <<README_EOF
Protein-ligand MD simulations (${NS} ns each) - results

In each folder:
  Protein_Ligand_Complex/  complex at 0 ns and ${NS} ns (protein + ligand) and input structures
  Results_Excel/           analysis_<sim>.xlsx (standard layout, every 0.25 ns + Mean/Min/Max),
                           Analysis_2_<sim>.xlsx (+ ligand columns, column definitions),
                           Protein_Ligand_MD_Analysis.xlsx (all time series), Summary_Statistics.xlsx
  Figures/                 19 figures, 600 dpi
  QC_REPORT.txt            quality-control checks (all included simulations: 0 FAIL)
  Methods_parameters.txt   simulation settings read from the input files
  Simulation_Inputs/       topology, ligand parameters, .mdp files, structures, run input, energies, logs
  ${AN}/       complete analysis output (raw .xvg, CSV, processed trajectory, commands)

Not included (size): full-system trajectories (.xtc). All values come directly from the trajectories.
Interaction energy = short-range Coulomb + Lennard-Jones (not a binding free energy).
README_EOF
echo; echo "included: $n   skipped: $skipped"
[ "$n" -gt 0 ] || { echo "!! nothing to package"; exit 1; }
cd "$(dirname "$OUT")" && rm -f "$(basename "$OUT").zip" && echo "zipping..." && zip -qr "$(basename "$OUT").zip" "$(basename "$OUT")"
ls -lh "$(dirname "$OUT")/$(basename "$OUT").zip"
