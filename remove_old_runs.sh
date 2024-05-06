#!/usr/bin/bash
# Remove runs older than 91 days that are not in use by any app (no *.toml file in the project directory).

# remove old runs
cd /ssd/data/ont
# get all run dirs older than 3 months
find  -mindepth 2 -maxdepth 2 -type d -mtime +91|grep -vP "persistence|pings|reads|user_scripts|queued_reads|intermediate|core-dump-db|_toRemove" > runs.txt
# sync and remove those that are not used by the app
while read d; do
    p=$(dirname $d);
    if [ $(find $p -name "*.toml"|wc -l) -eq 0 ]; then
	du -sh $d; rsync -a $d /users/enovoa/data/ont/$p; rm -r $d;
    fi;
done < runs.txt
