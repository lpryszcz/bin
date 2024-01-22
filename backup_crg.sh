#!/bin/bash

# move all completed bulk dumps to _toRemove
mkdir -p /ssd/data/ont/_toRemove
mv /ssd/data/ont/*.fast5 /ssd/data/ont/_toRemove
# sync everything once per day @ 6am
if [ `date +%H` -eq 6 ]; then    
    rsync -av /ssd/data/ont/ /users/enovoa/data/ont --exclude _toRemove \
	  --exclude "*.fast5.bam" --exclude "*.tar.gz" --exclude core-dump-db --exclude "*.tmp" \
	  --exclude "*_platform_qc*" --exclude intermediate --exclude queued_reads \
	  --exclude persistence --exclude pings --exclude reads --exclude user_scripts
# sync recent Fast5 files every hour
else
    find /ssd/data/ont/./ -mmin -61 -type f -name "*.fast5" | grep -v _toRemove > /tmp/rsync_crg.txt
    rsync -av --files-from=/tmp/rsync_crg.txt / /users/enovoa/data/ont
fi
