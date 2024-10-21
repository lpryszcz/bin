#!/bin/bash -l
# note you need to authenticate every 3 months or so..
module load treemap2

echo -e "<body><pre>\nHi!" > ~/du.sh.log

echo -e "\nBelow you can find a report of storage usage. Consider cleaning up, especially if you occupy a lot of space or generated millions of files." >> ~/du.sh.log
echo -e "You can find more details about your usage <a href='https://treemap.crg.es/treemap/dir'>here</a>." >> ~/du.sh.log
echo -e "For more general storage information, go <a href='http://www.linux.crg.es/index.php/Best_Practices'>here</a>." >> ~/du.sh.log

echo -e "\nTop per-dir usage from /users/enovoa" >> ~/du.sh.log
treemap -p /users/enovoa -o size > /tmp/treemap.out; head -n 15 /tmp/treemap.out >> ~/du.sh.log

echo -e "\nTop per-dir usage from /no_backup/enovoa" >> ~/du.sh.log
treemap -p /no_backup/enovoa -o size > /tmp/treemap.out
head -n 13 /tmp/treemap.out >> ~/du.sh.log
for d in nextflow_outputs users; do # analysis
    path="/no_backup/enovoa/$d"
    echo -e "\n- $path" >> ~/du.sh.log
    treemap -p $path -o size > /tmp/treemap.out
    head -n 15 /tmp/treemap.out|tail -n 6 >> ~/du.sh.log
done

path="/nfs/scratch01/enovoa"
echo -e "\nTop per-user usage from $path" >> ~/du.sh.log
treemap --tag scratch01 --disable-user-path -p $path --by-user -o size > /tmp/treemap.out
head -n 12 /tmp/treemap.out >> ~/du.sh.log

echo -e "\nTop per-dir usage from $path" >> ~/du.sh.log
treemap --tag scratch01 --disable-user-path -p $path -o size > /tmp/treemap.out
head -n 15 /tmp/treemap.out|tail -n 6 >> ~/du.sh.log

echo -e "\nHave a nice day!" >> ~/du.sh.log

echo -e "\n\nBTW: Leszek remember, you need to authenticate to treemap2 every 3 months or so..." >> ~/du.sh.log

# send an email
cat ~/du.sh.log | mail -s "$(echo -e "Quota usage from "`hostname`"\nContent-Type: text/html\n<html>")" lpryszcz@crg.es CRGLabEvaNovoa@crg.eu

