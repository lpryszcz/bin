#!/bin/bash
du -sh /no_backup_isis/enovoa/*/*/ 2> ~/du.sh.err.log | sort -hr | tee ~/du.sh.log
cat ~/du.sh.log | mail -s "Quota usage from "`hostname` lpryszcz@crg.es huanle.liu@crg.eu anna.delgado@crg.eu # eva.novoa@crg.eu
