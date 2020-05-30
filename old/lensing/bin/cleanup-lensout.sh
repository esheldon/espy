#!/bin/bash
# if we have to kill the jobs, we may have to clean up
# the files
#    /data/objshear/lensing/lensout/{sample}/lensout-*
#
# this runs through all machines and removes those files if
# they exist

if [[ $# -lt 1 ]]; then
    echo "cleanup-lensout run"
    exit 45
fi

run=$1

ls_command="ls /data/objshear/lensing/lensout/$run/"
rm_command="rm -r /data/objshear/lensing/lensout/$run/"

for i in $(seq -w 1 33); do

    machine="astro00$i"
    echo $machine
    ssh $machine $ls_command 2> /dev/null
    ssh $machine $rm_command 2> /dev/null

done
