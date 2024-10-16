#!/bin/bash

echo `date` Starting sync...
# move all completed bulk dumps to _toRemove
cd /ssd/data/ont
mkdir -p _toRemove; mv *.fast5 _toRemove

# sync files modified in last 61 minutes
find -L . -type f -cmin -61 -regextype posix-extended -regex ".*.(pod5|html|csv|txt|md|json)" | grep -vwP "_toRemove|reads" > /tmp/rsync_crg.txt
if [ -s /tmp/rsync_crg.txt ]; then
    echo `date` Syncing `cat /tmp/rsync_crg.txt|wc -l` files...
    rsync -av --files-from=/tmp/rsync_crg.txt `pwd` /no_backup/enovoa/data/ont
fi
echo `date` Finished.
