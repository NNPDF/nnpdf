# one-line command to remove in data point lines the comma and the percent sign
for f in $(ls *.csv); do sed -e "/^[0-9]/ s/[,%]/ /g" $f > $(basename "$f" .csv).dat; done 
