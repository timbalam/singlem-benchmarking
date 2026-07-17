# coverage4
echo -n > coverage4.tsv
for i in $(seq 1 25);
do
    echo "Otu$i	3" >> coverage4.tsv
done
# coverage5
echo -n > coverage5.tsv
j=0
for dist in 0.1 0.5 0.6 1 2 3 10 20 30 100 200 1000;
do
    j=$((j+1))
    echo "Otu$j	$dist" >> coverage5.tsv
done