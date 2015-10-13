#!/bin/bash

echo "Please input the theoryID you wish to upload, followed by [ENTER]:"

read theoryID
echo "Uploading Theory" $theoryID "...."

#If file exists, remove
if [ -f theory_$theoryID.tgz ] ; then
rm theory_$theoryID.tgz
fi

tar -cvzf ./theory_$theoryID.tgz ./theory_$theoryID
scp ./theory_$theoryID.tgz apfelcomb@pcteserver.mi.infn.it:~/WEB/commondatatheory/
rm ./theory_$theoryID.tgz

echo "Upload finished"
