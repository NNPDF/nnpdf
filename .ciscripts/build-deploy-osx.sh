#/usr/bin/bash
set -u
set -v

#Set up netrc file for uploading/downloading
echo "$NETRC_FILE" | base64 --decode > ~/.netrc

# Set the version for both vp and n3fit
gitroot=$(git rev-parse --show-toplevel)
gitversion=$(git describe --long --tags)
tag=$(git describe --abbrev=0 --tags)
gitversion=${gitversion/${tag}-/${tag}.}
githash=$(git rev-parse --short HEAD)
gitversion=${gitversion/-g${githash}/+g${githash}}
echo "build_version=\"${gitversion}\"" > ${gitroot}/n3fit/src/n3fit/version.py
echo "build_version=\"${gitversion}\"" > ${gitroot}/validphys2/src/validphys/version.py

conda build -q conda-recipe
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

echo "Uploading package to the NNPDF server"
#Idiotic mac mktemp
KEY=$( mktemp  "${TMPDIR:-/tmp}/key.XXXXXXXXX" )
#This is defined in the Gitlab variables, under the Settings Menu.
echo "$NNPDF_SSH_KEY" | base64 --decode > "$KEY"

scp -i "$KEY" -o StrictHostKeyChecking=no\
    ${CONDAPATH}/conda-bld/${OUTPUT_ARCH}/*.tar.bz2 \
    dummy@packages.nnpdf.science:~/packages/${OUTPUT_CHANNEL}/${OUTPUT_ARCH}

if [ $? == 0 ]; then
	echo "Upload suceeded"
else
	echo "Upload failed"
	exit 1
fi
