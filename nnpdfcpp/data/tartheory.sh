for i in theory_*/; do 
echo $i
tar -czf ${i:0:${#i}-1}.tgz $i; 
done
