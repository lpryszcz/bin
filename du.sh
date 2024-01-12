#!/bin/bash -l
# note you need to authenticate every 3 months or so..
#source /users/enovoa/lpryszcz/.bashrc
module load treemap2

echo "Top per-user usage from /users/enovoa" > ~/du.sh.log
treemap -p /users/enovoa --by-user -o size > /tmp/treemap.out; head -n 15 /tmp/treemap.out >> ~/du.sh.log

echo -e "\nTop per-dir usage from /no_backup/enovoa" >> ~/du.sh.log
#du -sh /no_backup/enovoa/*/*/ 2> ~/du.sh.err.log | sort -hr | head -n 25 >> ~/du.sh.log
treemap -p /no_backup/enovoa -o size > /tmp/treemap.out
head -n 15 /tmp/treemap.out >> ~/du.sh.log
for d in users analysis nextflow_outputs; do
    path="/no_backup/enovoa/$d"
    echo -e "\n- $path" >> ~/du.sh.log
    treemap -p $path -o size > /tmp/treemap.out
    head -n 15 /tmp/treemap.out|tail -n 6 >> ~/du.sh.log
done

path="/nfs/scratch01/enovoa"
echo -e "\nTop per-user usage from $path" >> ~/du.sh.log
treemap --tag scratch01 --disable-user-path -p $path --by-user -o size > /tmp/treemap.out
head -n 15 /tmp/treemap.out|tail -n 6 >> ~/du.sh.log

echo -e "\nTop per-dir usage from $path" >> ~/du.sh.log
treemap --tag scratch01 --disable-user-path -p $path -o size > /tmp/treemap.out
head -n 15 /tmp/treemap.out|tail -n 6 >> ~/du.sh.log

# send an email
cat ~/du.sh.log | mail -s "$(echo -e "Quota usage from "`hostname`"\nContent-Type: text/html\n<html><body><pre>")" lpryszcz@crg.es CRGLabEvaNovoa@crg.eu

