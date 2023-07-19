#!/bin/bash
MOUNT=$1
SIZE=$(df -h $MOUNT | awk '{print $4}' | tail -1 | awk '/K$/{print ($0+0)*1024*1024}/M$/{print ($0+0)/1024}/G$/{print $0+0}/T$/{print ($0+0)*1024}')
MINSIZE=$2 #in GB
HOST=$(hostname)

if (( $(echo "$SIZE < $MINSIZE" |bc -l) )); then
    echo "$MOUNT free space is $SIZE GB, which is below $MINSIZE GB Alert limit. Please release some space!" | mail -s "[$HOST] $MOUNT Low Free Space Alert: $SIZE GB" lpryszcz@crg.es $3 $4
fi
