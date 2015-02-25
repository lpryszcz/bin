#!/bin/bash
##find flash files
###http://askubuntu.com/questions/289003/where-does-chromium-keep-the-youtube-video-files
ls -l $1/fd/ 2>/dev/null 3>/dev/null| grep -i 'flash' 1>/dev/null  2>/dev/null 3>/dev/null;
if [ $? -eq "0" ]; 
  then 
  echo $1/fd/;
  ls -l $1/fd/ | grep -i 'flash';
fi
