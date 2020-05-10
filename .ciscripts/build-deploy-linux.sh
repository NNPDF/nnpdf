#!/bin/bash
#Find conda
source ~/.bashrc
set -e
set -o pipefail
set -u
set -v

#Set up netrc file for uploading/downloading
echo "$NETRC_FILE" | base64 --decode > ~/.netrc

#Build package
CONDA_PY=$CONDA_PY conda build -q conda-recipe
if [ $? != 0 ]; then
	echo failed to build
	exit 1
fi

if [ "${TRAVIS_PULL_REQUEST_BRANCH:-$TRAVIS_BRANCH}" != 'master'  ] && [ "$UPLOAD_NON_MASTER" == false ];
then
  	echo "
Skiping upload because this is not master and you have not
set the UPLOAD_NON_MASTER variable."
	exit 0
fi

#This seems to be needed for "artifacts" to work.
cp /root/miniconda3/conda-bld/linux-64/*.tar.bz2 .

echo "Uploading package to the NNPDF server"
KEY=$( mktemp )
#This is defined in the Travis environment variables.
echo "$NNPDF_SSH_KEY" | base64 --decode > "$KEY"

scp -i "$KEY" -o StrictHostKeyChecking=no\
    /root/miniconda3/conda-bld/linux-64/*.tar.bz2 \
    dummy@packages.nnpdf.science:~/packages/conda-private/linux-64

if [ $? == 0 ]; then
	echo "Conda package upload suceeded"
else
	echo "Conda package upload failed"
	exit 1
fi