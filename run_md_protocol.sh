#!/usr/bin/env bash
# run_md_protocol.sh "<simulation folder>" [--wait] [--no-analysis]
#
# Runs EM -> NVT -> NPT -> production MD in a folder prepared by setup_md_swissparam.sh,
# then runs run_md_analysis.sh on the finished trajectory.
#   --wait         wait until no other `gmx mdrun` is running on this machine before starting
#   --no-analysis  stop after production MD
# Stages that already finished are skipped; an interrupted production run is resumed from MD.cpt.
# Progress goes to <folder>/protocol_log.txt. To keep it running after the terminal closes:
#   nohup setsid bash run_md_protocol.sh "<folder>" --wait > /dev/null 2>&1 < /dev/null &
set -o pipefail; set +H; export LC_NUMERIC=C LC_ALL=C
TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
die(){ echo "!! STOP: $*"; exit 1; }
D_IN="${1:-}"; shift || true
WAIT=0; ANALYSE=1
for a in "$@"; do case "$a" in --wait) WAIT=1 ;; --no-analysis) ANALYSE=0 ;; *) die "unknown option $a" ;; esac; done
[ -n "$D_IN" ] && [ -d "$D_IN" ] || die "usage: run_md_protocol.sh <folder> [--wait] [--no-analysis]"
command -v gmx > /dev/null || die "gmx not found - activate your GROMACS environment first"
D="$(cd "$D_IN" && pwd)"; cd "$D" || exit 1
mkdir -p setup_logs
exec >> "$D/protocol_log.txt" 2>&1
echo "=== protocol started $(date '+%F %T') in $D"
for f in topol.top index.ndx EM.tpr NVT.mdp NPT.mdp MD.mdp; do [ -f "$f" ] || die "missing $f - run setup_md_swissparam.sh first"; done

if [ "$WAIT" = 1 ]; then
  echo "--- waiting for other mdrun processes to finish"
  while pgrep -f "gmx mdrun" > /dev/null; do sleep 120; done
  echo "--- no mdrun running $(date '+%F %T')"
fi
run(){ local L="$1"; shift; echo "--- $(date '+%F %T')  gmx $*"
  gmx "$@" > "setup_logs/$L.log" 2>&1 || { echo "!! FAILED $L"; tail -40 "setup_logs/$L.log"; exit 1; }
  grep -A4 -E '^WARNING [0-9]' "setup_logs/$L.log" | head -20; }
finished(){ [ -f "$1.gro" ] && grep "Finished mdrun" "$1.log" > /dev/null 2>&1; }

if finished EM; then echo "--- EM already finished, skipped"
else run em_run mdrun -deffnm EM; fi
grep -E "converged|Maximum force|Potential Energy" EM.log | tail -3

if finished NVT; then echo "--- NVT already finished, skipped"
else
  run nvt_grompp grompp -f NVT.mdp -c EM.gro -r EM.gro -p topol.top -n index.ndx -maxwarn 2 -o NVT.tpr
  run nvt_run mdrun -deffnm NVT
fi

if finished NPT; then echo "--- NPT already finished, skipped"
else
  run npt_grompp grompp -f NPT.mdp -c NVT.gro -r NVT.gro -t NVT.cpt -p topol.top -n index.ndx -maxwarn 2 -o NPT.tpr
  run npt_run mdrun -deffnm NPT
fi

if finished MD; then echo "--- production MD already finished, skipped"
elif [ -f MD.tpr ] && [ -f MD.cpt ]; then
  echo "--- $(date '+%F %T')  resuming production MD from MD.cpt"
  run md_resume mdrun -deffnm MD -cpi MD.cpt
else
  run md_grompp grompp -f MD.mdp -c NPT.gro -t NPT.cpt -p topol.top -n index.ndx -maxwarn 2 -o MD.tpr
  echo "--- $(date '+%F %T')  production MD started"
  run md_run mdrun -deffnm MD
fi
echo "=== production MD finished $(date '+%F %T')"

if [ "$ANALYSE" = 1 ]; then
  NS=$(gmx dump -s MD.tpr 2>/dev/null | awk '/^ *nsteps/{n=$3} /^ *dt /{d=$3} END{v=n*d/1000; if (v==int(v)) printf "%d", v}')
  [ -n "$NS" ] || die "production length is not a whole number of ns - run run_md_analysis.sh by hand"
  echo "--- analysis of $NS ns"
  bash "$TOOLS/run_md_analysis.sh" "$D" "$NS" > "$D/protocol_analysis.txt" 2>&1 \
    && echo "=== analysis finished - see MD_Analysis_${NS}ns/" \
    || echo "!! analysis stopped - see protocol_analysis.txt"
fi
echo "=== protocol finished $(date '+%F %T')"
