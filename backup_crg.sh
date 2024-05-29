#!/bin/bash

# move all completed bulk dumps to _toRemove
mkdir -p _toRemove
mv /ssd/data/ont/*.fast5 _toRemove
# sync everything once per day @ 6am
if [ `date +%H` -eq 6 ]; then
    for d in /ssd/data/ont/ /ssd_backup/data/ont/; do
	echo `date` $d
	rsync -av $d /users/enovoa/data/ont --exclude _toRemove --exclude "*.fast5.bam" --exclude "*.tar.gz" \
	  --exclude "*.tmp" --exclude core-dump-db --exclude intermediate --exclude queued_reads \
	  --exclude persistence --exclude pings --exclude reads --exclude user_scripts
    done
# sync recent Fast5 files every hour
else
    find /ssd/data/ont/./ -mmin -61 -type f -name "*.fast5" -o -name "*.pod5" | grep -v _toRemove > /tmp/rsync_crg.txt
    rsync -av --files-from=/tmp/rsync_crg.txt / /users/enovoa/data/ont
fi
