#!/usr/bin/env bash
# qsub_submit.sh <rule> <threads> <mem_mb> <runtime_min> <queue> <jobscript>
#
# Submit one Snakemake job to SGE and print ONLY the job id (qsub -terse), which
# is what the cluster-generic executor reads from stdout.
#
# Two conversions happen here, and both are the usual way to get SGE wrong.
#
#   Memory. Snakemake's mem_mb is the TOTAL the job needs. A consumable complex
#   is consumed PER SLOT under `-pe parallel N`, so requesting the total would
#   reserve N x total and the job would queue forever (or be admitted on a fat
#   node and mislead you about what it needs). We divide by the slot count,
#   rounding up. If your site's complex is consumable=JOB rather than per slot,
#   set SGE_MEM_PER_JOB=1 and the total is requested unchanged.
#
#   Wall time. runtime is minutes; h_rt is HH:MM:SS. A job whose h_rt exceeds the
#   ceiling of its queue is rejected at submission; one with no h_rt inherits the
#   queue default and is killed silently when it runs over.
#
# The memory complex name is site-specific, and only a CONSUMABLE complex
# actually reserves anything. Check with
#     qconf -sc | grep -E 'h_vmem|mem_free|m_mem_free'
# On the cluster this was written for, h_vmem is consumable and mem_free is not,
# so h_vmem is the default here. Note that h_vmem caps address space, not
# resident memory: a tool that maps far more than it touches can be killed with a
# small RSS, in which case raise the request rather than doubting the tool.
#
#   SGE_MEM_COMPLEX=none      request no memory at all (schedule purely on slots)
#   SGE_MEM_FREE_TOO=1        additionally request mem_free, which reserves
#                             nothing but keeps the job off a node whose RAM is
#                             already committed
#   SGE_MEM_PER_JOB=1         complex is consumable=JOB, not per slot
set -euo pipefail

rule=$1
threads=$2
mem_mb=$3
runtime=$4
queue=$5
jobscript=$6

PE=${SGE_PE:-parallel}
MEM_COMPLEX=${SGE_MEM_COMPLEX:-h_vmem}
LOGDIR=${SGE_LOGDIR:-.snakemake/sge_logs}
mkdir -p "$LOGDIR"

[ "$threads" -ge 1 ] 2>/dev/null || threads=1

# total MB -> per-slot MB, rounded up (unless the complex is per job)
if [ "${SGE_MEM_PER_JOB:-0}" = "1" ]; then
  mem_req=$mem_mb
else
  mem_req=$(( (mem_mb + threads - 1) / threads ))
fi

# minutes -> HH:MM:SS
h=$(( runtime / 60 )); m=$(( runtime % 60 ))
h_rt=$(printf '%02d:%02d:00' "$h" "$m")

args=(-terse -cwd -V -j y -N "smk.${rule}" -o "$LOGDIR" -q "$queue" -l "h_rt=${h_rt}")
[ "$threads" -gt 1 ] && args+=(-pe "$PE" "$threads")
if [ "$MEM_COMPLEX" != "none" ]; then
  args+=(-l "${MEM_COMPLEX}=${mem_req}M")
  # mem_free is not consumable here, so this reserves nothing; it only stops the
  # job landing on a node whose memory is already spoken for.
  [ "${SGE_MEM_FREE_TOO:-0}" = "1" ] && args+=(-l "mem_free=${mem_req}M")
fi

# -terse prints the bare job id; anything else on stdout breaks the executor
qsub "${args[@]}" "$jobscript" | tr -d '[:space:]'
