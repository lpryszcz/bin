#free -m | head -n2 | tail -1 | awk '{ print $7"Mb ["$7*100/$2"%]" }'
free -h | head -n2 | tail -1 | awk '{print $7"/"$2}'
