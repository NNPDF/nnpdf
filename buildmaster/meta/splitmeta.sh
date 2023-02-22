BASE=../inc
# Split out all dataInfoRaw structs from headers, removing commas
for i in $BASE/*.h; do
    awk '/static const/ {for(i=1; i<=4; i++) {getline; print}}' $i | tr -d , >> split1.meta
done
# Remove comments
paste split1.meta | cut -f1 -d"/" > split.meta
rm split1.meta
# Separate all dataInfoRaw structs into separate files
split -l 4 split.meta pre_meta
# Attatch the yaml template
for i in pre_meta*; do
    paste template.yaml $i  > processed_${i}
done
rm pre_meta*
# Rename the files to the correct name
for i in processed_*; do
    SETNAME=$(paste $i | grep setname | awk '{print $2}' | tr -d '"')
    echo "Extracted metadata for "$SETNAME
    mv $i ${SETNAME}.yaml
done
