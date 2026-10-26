#!/usr/bin/env bash
# Read-only status digest for a dockpipe run. Pure bash/awk on purpose: no
# embedded Python, so copying this around cannot break on indentation.
#
#   ./check_run.sh [run_name]        # default run name: panel500_v1
#
# Safe to run while the job is still going.

# Not -e: every section must still run when an earlier one has no data yet.
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
RUN="${1:-panel500_v1}"
R="results/$RUN"

# Columns of scores*.csv, in order. Only 1-9 are parsed with -F, -- `error`
# (10) onwards can contain commas, everything before it cannot.
C_RECEPTOR=2
C_BEST=4
C_SECONDS=6
C_STATUS=9

hr() { printf '\n===== %s =====\n' "$1"; }

# n / median / p90 / max from numbers on stdin.
stats() {
  sort -n | awk -v label="$1" '
    { v[NR] = $1 }
    END {
      if (NR == 0) { printf "  %s: no values\n", label; exit }
      p90 = int(NR * 0.9); if (p90 < 1) p90 = 1
      printf "  %s: n=%d median=%.0f p90=%.0f max=%.0f\n",
             label, NR, v[int((NR + 1) / 2)], v[p90], v[NR]
    }'
}

hr "slurm"
if command -v sacct >/dev/null 2>&1; then
  sacct -u "$USER" \
    --format=JobID%14,JobName%12,State%12,Elapsed%12,ExitCode,MaxRSS \
    -X 2>/dev/null | tail -6
else
  echo "  sacct not available (not on the cluster?)"
fi

LOG=$(ls -t logs/dockpipe_*.log 2>/dev/null | head -1)
hr "log: ${LOG:-none found}"
if [[ -n "$LOG" ]]; then
  tail -20 "$LOG"
  hr "log markers"
  grep -aE '^(env|python|commit|config|workers) |^=== |Traceback|rror|FAILED|Killed|CANCELLED|slurmstepd|srun:' \
    "$LOG" | tail -25
fi

hr "ligands"
IDX="$R/ligands/ligand_index.json"
CKL="$R/ligands/.ckpt_ligands.jsonl"
if [[ -f "$IDX" ]]; then
  N_REC=$(grep -c '"smiles"' "$IDX")
  N_ST=$(grep -o '"n_states": *[0-9]*' "$IDX" | awk '{s += $2} END {print s+0}')
  echo "  index: $N_REC records, $N_ST ligand states"
  grep -o '"status": *"[^"]*"' "$IDX" | sed 's/.*"\([^"]*\)"$/  \1/' \
    | sort | uniq -c | sort -rn
else
  echo "  no ligand_index.json yet"
fi
[[ -f "$CKL" ]] && echo "  checkpoint: $(wc -l < "$CKL") lines"

hr "redock"
RC="$R/docking/redock_control.json"
if [[ -f "$RC" ]]; then
  grep -E '"(receptor|rmsd_reported|rmsd_best|best_pose_rank|verdict|bond_orders)"' \
    "$RC" | sed 's/^ *//; s/,$//'
else
  echo "  no redock_control.json yet"
fi

hr "dock progress"
shopt -s nullglob
CKPTS=("$R"/docking/.ckpt_*.jsonl)
if ((${#CKPTS[@]})); then
  for f in "${CKPTS[@]}"; do
    printf '  %-52s %6s done\n' "$(basename "$f")" "$(wc -l < "$f")"
  done
else
  echo "  no dock checkpoints yet"
fi
ls -la "$R"/docking/scores*.csv 2>/dev/null | sed 's/^/  /'

hr "scores"
SRC=()
for f in "$R"/docking/scores*.csv; do
  case "$f" in *by_ligand*|*wide*|*merged*) continue ;; esac
  SRC+=("$f")
done

if ((${#SRC[@]} == 0)); then
  echo "  no score CSVs yet"
else
  TMP=$(mktemp) || exit 1
  trap 'rm -f "$TMP"' EXIT
  for f in "${SRC[@]}"; do
    tail -n +2 "$f" >> "$TMP"
  done
  echo "  files: ${SRC[*]##*/}"
  echo "  rows : $(wc -l < "$TMP")"

  echo "  status:"
  awk -F, -v c=$C_STATUS 'NF >= c { print "    " $c }' "$TMP" \
    | sort | uniq -c | sort -rn

  awk -F, -v c=$C_SECONDS 'NF >= c && $c != "" { print $c }' "$TMP" \
    | stats "seconds"

  echo "  per receptor (status=ok):"
  for rec in $(awk -F, -v c=$C_RECEPTOR 'NF >= c { print $c }' "$TMP" | sort -u); do
    awk -F, -v r="$rec" -v cr=$C_RECEPTOR -v cs=$C_STATUS -v cb=$C_BEST \
      'NF >= cs && $cr == r && $cs == "ok" && $cb != "" { print $cb }' "$TMP" \
      | sort -n \
      | awk -v r="$rec" '
          { v[NR] = $1 }
          END {
            if (NR == 0) { printf "    %-22s n=0\n", r; exit }
            printf "    %-22s n=%5d best=%7.2f median=%7.2f worst=%7.2f\n",
                   r, NR, v[1], v[int((NR + 1) / 2)], v[NR]
          }'
  done
fi

hr "merged outputs"
for f in scores_merged.csv scores_by_ligand.csv scores_wide.csv; do
  if [[ -f "$R/docking/$f" ]]; then
    printf '  %-22s %s rows\n' "$f" "$(( $(wc -l < "$R/docking/$f") - 1 ))"
  else
    printf '  %-22s missing\n' "$f"
  fi
done

echo
