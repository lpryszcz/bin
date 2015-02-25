#!/bin/bash

#umount
echo "Unmounting sshfs..."
mountpoints=`mount | grep sshfs | cut -f3 -d" "`
for m in $mountpoints; do
    echo " $m"
    sudo umount $m
done

#fusermount
mountpoints=`mount | grep sshfs | cut -f3 -d" "`
if [ ! -z $mountpoints ]; then
    echo "Unmounting sshfs by 'fusermount -uz'..."
    for m in $mountpoints; do
        echo " $m"
        sudo fusermount -uz $m
    done
fi

#kill audio/video apps and sshfs
echo "Terminating sshfs"
ps aux | grep -P "sshfs" | grep -v grep | awk '{print $2}' | xargs kill -s 9

#drop disk cache if little memory
#sync && echo 3 | sudo tee /proc/sys/vm/drop_caches && 

#isilonup
sshfs -o "idmap=user,reconnect,workaround=all,compression=yes,ssh_command=ssh" ssh-server.crg.es:/nfs/users /mnt/users


