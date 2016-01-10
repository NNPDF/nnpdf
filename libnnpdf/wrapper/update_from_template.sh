#!/bin/bash
# update_from_template.sh
# @author Zahari Kassabov
#
#
# This script allows to easily make a change in the template file
#    swig.template
# and then propagate the result to the existing wrappers, keeping their
# respective changes from the template.
#
# It is implemented by performing a three-way git merge between:
#
# - The wrappers generated with the current template (common parent).
# - The existing wrapper sources.
# - The wrappers generated with the edited template.
#
# The script will initialize a repository in a temporary folder and
# apply the first two, and the priompt the user to edit the template.
# It can happen that the merge cannot be done automaticallt, in which
# case the conflicts will have to be resolved with a git mergetool
# such as `meld`.
#


if [[ `git diff HEAD --name-only swig.template` ]]; then
	printf "\x1b[1m\x1b[31m"
	#Any reasonable way to break the lines?
    printf "swig.template has been modified since the last commit. You are supposed to edit it using this script."
	printf "\x1b(B\x1b[m\n"
	exit 1
fi
THIS_PATH=$PWD
TEMP_PATH=`mktemp -d -t gitmanglingXXXXXXXXX`


echo "Created temporary directory $TEMP_PATH"
cd $TEMP_PATH
python $THIS_PATH/make_base.py
cd stubs
git init .
git add *.i
git commit -am "Old template"
git tag initial
cp $THIS_PATH/src/*.i .
git commit -am "Current wrappers"
git tag current
git checkout -q initial

printf "\x1b[1m\x1b[33m"
echo "Now we will edit the template."
_EDITOR=`git config core.editor`
if [[ ! $_EDITOR ]]; then
    if [[ $EDITOR ]]; then
	    _EDITOR="$EDITOR"
    else
		_EDITOR="vi"
	fi
fi
read -p "Select editor (default: $_EDITOR) " _USER_EDITOR
printf "\x1b(B\x1b[m"
if [[ $_USER_EDITOR ]]; then
	_EDITOR="$_USER_EDITOR"
fi
$_EDITOR $THIS_PATH/swig.template
$(cd ..; python $THIS_PATH/make_base.py)

git commit -am "Updated template"
git merge current -m "Merging"
sucess=$?
if [[ "$sucess" -ne 0 ]]; then
	printf "\x1b[1m\x1b[31m"
	echo -e "Automatic merge was failed. Fix contlicts in the folder\n$TEMP_PATH/stubs\n and then copy the result."
	printf "\x1b(B\x1b[m"
else
	printf "\x1b[1m\x1b[32m"
	echo -e "Automatic merge sucessful. Updated files are in\n$TEMP_PATH/stubs\n You can update them with\n mv $TEMP_PATH/stubs/* ./src"
	printf "\x1b(B\x1b[m"
fi
