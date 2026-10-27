#!/bin/bash
# Keep the own Slurm jobs from crowding out other users while long runs are going.
#
# Every few minutes the script looks at the partition and decides how many of the own jobs may be
# running or waiting to run. Pending jobs beyond that number are put on hold and released again
# when there is room. Running jobs are never touched, and Nextflow simply sees held jobs as pending.
#
#   other users have jobs waiting for resources    ->  LOW   (default 24)
#   fewer than MIN_IDLE cpus idle (default 48)      ->  MID   (default 60)
#   otherwise                                       ->  HIGH  (default 150)
#
# Usage: slurm_throttle.sh <log file> [<pattern of processes to wait for>]
# The script ends when no process matching the pattern (default: run_profiles.sh) is left and the
# user has no jobs in the queue. All held jobs are released at the end.
# Environment: PARTITION (exbio-cpu), OTHER_PARTITIONS (shared-cpu,exbio-cpu,exbio-long-jobs),
#              LOW, MID, HIGH, MIN_IDLE, INTERVAL (seconds, default 300)

set -u
log=${1:?log file}
pattern=${2:-run_profiles.sh}
partition=${PARTITION:-exbio-cpu}
other_partitions=${OTHER_PARTITIONS:-shared-cpu,exbio-cpu,exbio-long-jobs}
low=${LOW:-24}; mid=${MID:-60}; high=${HIGH:-150}; min_idle=${MIN_IDLE:-48}; interval=${INTERVAL:-300}
me=$(id -un)

release_all() {
    held=$(squeue -h -u "$me" -t PD -o "%i %r" | awk '$2 == "JobHeldUser" {print $1}')
    [ -n "$held" ] && scontrol release $held 2> /dev/null
}
trap 'release_all; exit 0' TERM INT

while true; do
    running=$(squeue -h -u "$me" -t R | wc -l)
    waiting=$(squeue -h -u "$me" -t PD -o "%i %r" | awk '$2 != "JobHeldUser" {print $1}' | sort -n)
    held=$(squeue -h -u "$me" -t PD -o "%i %r" | awk '$2 == "JobHeldUser" {print $1}' | sort -n)
    others_waiting=$(squeue -h -t PD -p "$other_partitions" -o "%u %r" | awk -v me="$me" '$1 != me && ($2 == "Resources" || $2 == "Priority")' | wc -l)
    idle=$(sinfo -h -p "$partition" -o "%C" | cut -d/ -f2)

    if [ "$others_waiting" -gt 0 ]; then cap=$low; level=LOW
    elif [ "${idle:-0}" -lt "$min_idle" ]; then cap=$mid; level=MID
    else cap=$high; level=HIGH; fi

    room=$(( cap - running )); [ $room -lt 0 ] && room=0
    n_waiting=$(echo "$waiting" | grep -c .)
    n_held=$(echo "$held" | grep -c .)
    action="none"
    if [ "$n_waiting" -gt "$room" ]; then
        to_hold=$(echo "$waiting" | tail -n $(( n_waiting - room )))
        scontrol hold $to_hold 2> /dev/null
        action="held $(echo "$to_hold" | grep -c .)"
    elif [ "$n_held" -gt 0 ] && [ "$n_waiting" -lt "$room" ]; then
        to_release=$(echo "$held" | head -n $(( room - n_waiting )))
        scontrol release $to_release 2> /dev/null
        action="released $(echo "$to_release" | grep -c .)"
    fi
    echo "$(date '+%F %T') level=$level cap=$cap running=$running waiting=$n_waiting held=$n_held others_waiting=$others_waiting idle_cpus=${idle:-?} action=$action" >> "$log"

    # the pattern is also part of the command line of this script, so its own process is left out
    if [ -z "$(pgrep -u "$me" -f "$pattern" | grep -vx "$$")" ] && [ "$(squeue -h -u "$me" | wc -l)" -eq 0 ]; then
        break
    fi
    sleep "$interval"
done
release_all
echo "$(date '+%F %T') finished" >> "$log"
