#!/bin/bash
# Clean space from /
# has to be run as root

# clean /var/apt
sudo apt-get clean

# retain max 2 version of each snap - this can't be lower than 2
sudo snap set system refresh.retain=2
# clean /var/snapd
#CLOSE ALL SNAPS BEFORE RUNNING THIS
clean_snap.sh

