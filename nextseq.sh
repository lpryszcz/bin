#!/bin/bash
###
# Cron service for NextSeq to FastQ conversion. 
# Will send email if any WARNING encountered. 
###

me=`basename $0`
host=`hostname`

# load env variables 
## make sure line `[ -z "$PS1" ]` is commented in .bashrc
## http://stackoverflow.com/a/15574078/632242
source /home/lpryszcz/.bashrc

EMAILMESSAGE="/tmp/nextseq.log"

# Email text/message
echo "Hi!" > $EMAILMESSAGE
echo `date` "Converting NextSeq runs..." >> $EMAILMESSAGE
echo "" >> $EMAILMESSAGE

# run nextseq.py
nextseq.py -v >> $EMAILMESSAGE 2>&1

echo `date` "Done!" >> $EMAILMESSAGE
echo "" >> $EMAILMESSAGE

###
# script to send simple email
# email subject
SUBJECT="[nextseq]"
# Email To ?
EMAIL="lpryszcz@iimcb.gov.pl"
echo "Greetings from $host :)" >> $EMAILMESSAGE

# send an email using /bin/mail only if WARNING in $EMAILMESSAGE
grep=`grep WARNING $EMAILMESSAGE`
if [ ! -z "$grep" ]; then 
    mail -s "$SUBJECT" "$EMAIL" < $EMAILMESSAGE
fi
