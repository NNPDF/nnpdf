#/usr/bin/bash

#Find conda
source ~/.bashrc
conda build conda-recipe
if [ $? != 0 ]; then
	echo failed to build
	exit 1
fi

cp /yes/conda-bld/linux-64/*.tar.bz2 .

KEY=$( mktemp )
#This is defined in the Gitlab variables, under the Settings Menu.
echo "$ZIGZAH_SSH_KEY" > "$KEY"

scp -i "$KEY" -o StrictHostKeyChecking=no\
    /yes/conda-bld/linux-64/*.tar.bz2 \
    dummy@zigzah.com:~/conda-pkgs-private/linux-64 
