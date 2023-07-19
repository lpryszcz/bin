#!/bin/bash
# note you need to authenticate every 3 months or so..
module load treemap2

echo "Top per-user usage from /users/enovoa" > ~/du.sh.log
treemap -p /users/enovoa --by-user -o size > /tmp/treemap.out; head -n 15 /tmp/treemap.out >> ~/du.sh.log

echo -e "\nTop per-dir usage from /no_backup/enovoa" >> ~/du.sh.log
#du -sh /no_backup/enovoa/*/*/ 2> ~/du.sh.err.log | sort -hr | head -n 25 >> ~/du.sh.log
treemap -p /no_backup/enovoa -o size > /tmp/treemap.out; head -n 15 /tmp/treemap.out >> ~/du.sh.log
for d in users analysis nextflow_outputs; do
    path="/no_backup/enovoa/$d"
    echo -e "\n- $path" >> ~/du.sh.log
    treemap -p $path -o size > /tmp/treemap.out; head -n 15 /tmp/treemap.out|tail -n 6 >> ~/du.sh.log
done

echo -e "\nTop usage from /nfs/scratch01/enovoa" >> ~/du.sh.log
du -sh /nfs/scratch01/enovoa/*/ 2>> ~/du.sh.err.log | sort -hr | head -n 25 >> ~/du.sh.log
#treemap -p /nfs/scratch01/enovoa --by-user -o size

# send an email
cat ~/du.sh.log | mail -s "$(echo -e "Quota usage from "`hostname`"\nContent-Type: text/html\n<html><body><pre>")" CRGLabEvaNovoa@crg.eu #lpryszcz@crg.es

