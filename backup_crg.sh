#!/bin/bash

# move all completed bulk dumps to _toRemove
cd /ssd/data/ont
mkdir -p _toRemove; mv *.fast5 _toRemove
# consider -ctime 0 # files created in the last day (24h)
find -L . -type f -cmin -61 -regextype posix-extended -regex ".*.(pod5|html|csv|txt|md|json)" | grep -vwP "_toRemove|reads" > /tmp/rsync_crg.txt
rsync -av --files-from=/tmp/rsync_crg.txt `pwd` /users/enovoa/data/ont
