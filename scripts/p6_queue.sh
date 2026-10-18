#!/usr/bin/env bash
# Two-GPU job queue for P6 (one worker per GPU). Each line of the job file is a command template where {gpu}
# is replaced by the worker's GPU. A job is claimed atomically (mkdir); finished jobs leave claims/<n>/DONE, failed
# jobs claims/<n>/FAILED (delete the claim directory to retry). Re-running the script skips finished jobs, so it
# can be restarted after a reboot.
#
#   bash scripts/p6_queue.sh jobs/p6_e1_main.txt queues/e1_main ["0 1"]   # optional GPU list
set -u
cd "$(dirname "$0")/.."
JOBS=$1
STATE=$2
GPUS=${3:-"0 1"}
mkdir -p "$STATE/claims" "$STATE/logs"
export OMP_NUM_THREADS=2

worker() {
  local gpu=$1 pass n rc
  # stop quietly when the machine shuts down: the job resumes from its checkpoint after the reboot
  trap 'exit 0' TERM INT
  for pass in 1 2; do            # second pass picks up lines released by interrupted jobs
    n=0
    while IFS= read -r line || [ -n "$line" ]; do
      n=$((n + 1))
      [ -z "$line" ] && continue
      [[ "$line" == \#* ]] && continue
      mkdir "$STATE/claims/$n" 2>/dev/null || continue
      cmd=${line//\{gpu\}/$gpu}
      echo "[$(date '+%F %T')] gpu$gpu start $n: $cmd" >> "$STATE/queue.log"
      taskset -c 0-11,14-31 bash -c "$cmd" > "$STATE/logs/$n.log" 2>&1
      rc=$?
      if [ $rc -eq 0 ]; then
        touch "$STATE/claims/$n/DONE"
        echo "[$(date '+%F %T')] gpu$gpu done  $n" >> "$STATE/queue.log"
      elif [ $rc -ge 128 ]; then   # killed by a signal (shutdown): release the claim, do not mark FAILED
        rm -rf "$STATE/claims/$n"
        echo "[$(date '+%F %T')] gpu$gpu interrupted $n (rc $rc), released" >> "$STATE/queue.log"
        exit 0
      else
        touch "$STATE/claims/$n/FAILED"
        echo "[$(date '+%F %T')] gpu$gpu FAIL  $n (rc $rc)" >> "$STATE/queue.log"
      fi
    done < "$JOBS"
  done
}

for g in $GPUS; do
  worker "$g" &
done
wait
echo "[$(date '+%F %T')] queue finished" >> "$STATE/queue.log"
touch "$STATE/DONE"
