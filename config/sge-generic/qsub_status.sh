#!/usr/bin/env bash
# qsub_status.sh <jobid>
#
# Print exactly one of: running | success | failed
#
# SGE forgets a job the moment it leaves the queue, and the accounting record
# that replaces it appears a few seconds later. So between `qstat` losing the job
# and `qacct` gaining it there is a window in which the job exists nowhere. A
# status script that reports "failed" in that window will kill a healthy run.
# We therefore report "running" until qacct can be read, and only after
# SGE_STATUS_ATTEMPTS consecutive misses conclude the job is genuinely gone.
#
# A job killed for exceeding h_rt or h_vmem shows failed=1 with exit_status=137
# (or 100+signal); a normal failure shows exit_status != 0. Both are "failed".
set -uo pipefail

jobid=$1
ATTEMPTS=${SGE_STATUS_ATTEMPTS:-20}
SLEEP=${SGE_STATUS_SLEEP:-5}
STATEFILE="${TMPDIR:-/tmp}/smk_sge_missing_${jobid}"

# still queued or running?
if qstat -j "$jobid" >/dev/null 2>&1; then
  rm -f "$STATEFILE"
  echo "running"
  exit 0
fi

# finished, deleted, or the accounting record has not landed yet
acct=$(qacct -j "$jobid" 2>/dev/null || true)
if [ -z "$acct" ]; then
  n=$(( $(cat "$STATEFILE" 2>/dev/null || echo 0) + 1 ))
  echo "$n" > "$STATEFILE"
  if [ "$n" -ge "$ATTEMPTS" ]; then
    rm -f "$STATEFILE"
    echo "failed"            # gone from qstat and never appeared in qacct
  else
    sleep "$SLEEP"
    echo "running"           # accounting lag, not a failure
  fi
  exit 0
fi
rm -f "$STATEFILE"

exit_status=$(awk '/^exit_status/ {print $2; exit}' <<< "$acct")
failed=$(awk '/^failed/ {print $2; exit}' <<< "$acct")

if [ "${exit_status:-1}" = "0" ] && [ "${failed:-1}" = "0" ]; then
  echo "success"
else
  echo "failed"
fi
