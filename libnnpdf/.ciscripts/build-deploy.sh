#/usr/bin/bash

#Find conda
source ~/.bashrc
conda build -q .ciscripts/branch-metapackage
if [ $? != 0 ]; then
	echo failed to build metapackage
	exit 1
fi

conda build -q conda-recipe
if [ $? != 0 ]; then
	echo failed to build
	exit 1
fi

#This seems to be needed for "artifacts" to work.
cp /root/miniconda3/conda-bld/linux-64/*.tar.bz2 .

echo "Uploading package to zigzah"
KEY=$( mktemp )
#This is defined in the Gitlab variables, under the Settings Menu.
echo "$ZIGZAH_SSH_KEY" > "$KEY"

scp -i "$KEY" -o StrictHostKeyChecking=no\
    /root/miniconda3/conda-bld/linux-64/*.tar.bz2 \
    dummy@zigzah.com:~/conda-pkgs-private/linux-64 

if [ $? == 0 ]; then
	echo "Upload suceeded"
else
	echo "Upload failed"
	exit 1
fi
