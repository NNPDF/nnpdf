#!/bin/sh

cp ./github.css ./plts/

cd plts
rm index.md
rm index.html

echo "NNPYDF Plot Index" >> index.md
echo "-----------------" >> index.md

I=1
for path in ./*; do
    [ -d "${path}" ] || continue # if not a directory, skip
    dirname="$(basename "${path}")"
    echo $dirname

    markdown "${path}"/index.md > "${path}"/index.html

    echo "<head>" >> "${path}"/index.html
    echo "<link rel=\"stylesheet\" type=\"text/css\" href=\"../github.css\">" >> "${path}"/index.html
    echo "</head>" >> "${path}"/index.html

    echo ${I}". ["${dirname}"]("${path}"/index.html)" >> index.md
    I=$((I+1))
done

markdown index.md >> index.html

echo "<head>" >> index.html
echo "<link rel=\"stylesheet\" type=\"text/css\" href=\"github.css\">" >> index.html
echo "</head>" >> index.html
cd ..
