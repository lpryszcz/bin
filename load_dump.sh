#!/bin/bash
# Logs (ps aux) to TMP processes with >10% CPU or >10% of RAM usage
# USER       PID %CPU %MEM    VSZ   RSS TTY      STAT START   TIME COMMAND

LOG=/tmp/load_dump.log

echo `date` >> $LOG
ps aux | awk '$3>10 || $4>10' >> $LOG
