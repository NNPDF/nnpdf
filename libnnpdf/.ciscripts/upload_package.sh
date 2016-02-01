KEY=$( mktemp )
echo "$ZIGZAH_SSH_KEY" > "$KEY"
scp -i "$KEY" -o StrictHostKeyChecking=no\
    yes/conda-bld/linux-64/*.tar.bz2 \
    dummy@zigzah.com:~/conda-pkgs-private/linux-64 
