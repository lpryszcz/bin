#!/bin/bash
# Put runs older than given number of DAYS from ONTDIR onto tape.
# Run has to have a pod5 directory to be stored.
# Text files (csv, txt, json,md) are gzip compressed. 
# Runs archived before are marked with .tape file and are skipped.

DAYS=7
ONTDIR=/no_backup/enovoa/data/ont
echo `date` Archiving runs in $ONTDIR older than $DAYS days...

cd ~/_logs
header="#!/bin/sh\n[operation]\naction=archive\ncompression=no\n[objects]"
jobfn=cron.$(date +%y%m%d).sh
echo -e $header > $jobfn

# find runs in ONTDIR
for f in $(find $ONTDIR -type d -mtime +$DAYS -name pod5); do
    d=$(dirname $f)
    if [ ! -e $d/.tape ]; then
	du -sh $d
	gzip $d/*.{csv,txt,json,md} $d/other_reports/*.csv
	echo $d >> $jobfn
	touch $d/.tape
    fi
done

# start archive job only if some dirs are in the job file
if [ $(cat $jobfn |wc -l) -gt 5 ]; then
    echo `date` Submitting archiving job...
    abatch --mail-type=ALL --mail-user="leszek.pryszcz@crg.eu" $jobfn
    aqueue
fi

echo `date` Done!
