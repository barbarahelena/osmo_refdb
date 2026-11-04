#!/usr/bin/env bash
# harvest_bacdive.sh <outdir> [sample|sweep] [max_id]
#
# Fetch BacDive strain records into <outdir>/bacdive_raw/ as 100-ID batches.
#   sample (default) -- a seeded random draw of 1000 IDs, reproducible, for
#                       coverage estimates and unbiased effect reads
#   sweep            -- every ID up to max_id (default 175000), ~1750 calls
#
# Batches already present and parseable are skipped, so an interrupted run
# resumes. Truncated files (the failure mode when a run is killed mid-write)
# are detected by a JSON parse check and re-fetched rather than silently used.
#
# No credentials: registration was dropped in Feb 2026. v2 endpoints only --
# v1 still responds but is frozen at April 2025 data.
set -uo pipefail

PV=${1:?usage: harvest_bacdive.sh <outdir> [sample|sweep] [max_id]}
MODE=${2:-sample}
MAXID=${3:-175000}
RAW="$PV/bacdive_raw"
API=https://api.bacdive.dsmz.de/v2/fetch
mkdir -p "$RAW"

case "$MODE" in
  sample)
    # seed 42 over range(1, 175001) reproduces the pilot's exact 1000 IDs
    python3 -I -c "
import random
random.seed(42)
ids = random.sample(range(1, $MAXID + 1), 1000)
for i in range(10):
    print(';'.join(str(x) for x in ids[i*100:(i+1)*100]))
" > "$PV/batches.txt"
    ;;
  sweep)
    python3 -I -c "
ids = list(range(1, $MAXID + 1))
for i in range(0, len(ids), 100):
    print(';'.join(str(x) for x in ids[i:i+100]))
" > "$PV/batches.txt"
    ;;
  *) echo "unknown mode: $MODE (expected sample or sweep)" >&2; exit 2 ;;
esac

TOTAL=$(wc -l < "$PV/batches.txt" | tr -d ' ')
echo "mode=$MODE batches=$TOTAL -> $RAW"

ok() { python3 -I -c "import json,sys; json.load(open(sys.argv[1]))['results']" "$1" 2>/dev/null; }

n=0; fetched=0; cached=0; failed=0
while read -r line; do
  n=$((n+1)); f="$RAW/b$n.json"
  if [ -s "$f" ] && ok "$f"; then cached=$((cached+1)); continue; fi
  if curl -sS -m 180 --retry 3 --retry-delay 2 "$API/$line" -o "$f" && ok "$f"; then
    fetched=$((fetched+1))
  else
    failed=$((failed+1)); echo "  batch $n FAILED" >&2; rm -f "$f"
  fi
  # pace the sweep: this is a public academic service, not our infrastructure
  [ "$MODE" = sweep ] && python3 -I -c "import time; time.sleep(0.3)"
  [ $((n % 50)) -eq 0 ] && echo "  $n/$TOTAL (fetched=$fetched cached=$cached failed=$failed)"
done < "$PV/batches.txt"

echo "done: fetched=$fetched cached=$cached failed=$failed"
python3 -I -c "
import json, glob, sys
n = sum(json.load(open(f))['count'] for f in glob.glob('$RAW/b*.json'))
print(f'strains resolved: {n}')
"
