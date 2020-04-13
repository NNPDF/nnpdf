#!/bin/bash
#Find conda
source ~/.bashrc
set -e
set -o pipefail
set -u
set -v

#Set up netrc file for uploading/downloading
echo "$NETRC_FILE" | base64 --decode > ~/.netrc

# build documentation
conda config --add channels https://packages.nnpdf.science/conda;
conda config --add channels https://packages.nnpdf.science/conda-private;
conda config --set show_channel_urls true;
conda install nnpdf --yes
cd doc/sphinx
make html

echo "Uploading documentation to the NNPDF server"
KEY=$( mktemp )
#This is defined in the Travis environment variables.
echo "$NNPDF_SSH_KEY" | base64 --decode > "$KEY"

scp -r -i "$KEY" -o StrictHostKeyChecking=no\
    build/html/* \
    dummy@packages.nnpdf.science:~/sphinx-docs/

if [ $? == 0 ]; then
	echo "Documentation upload suceeded"
else
	echo "Documentation upload failed"
	exit 1
fi
