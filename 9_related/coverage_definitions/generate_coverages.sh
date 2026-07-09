echo -n > coverage4.tsv
for i in $(seq 1 60);
do
    echo "Otu$i	3.0" >> coverage4.tsv
done
echo -n > coverage5.tsv
j=0
for dist in 0.1 0.5 0.6 1.0 2.0 3.0 10.0 20.0 30.0 100.0 200.0 1000.0;
do
    j=$((j+1))
    echo "Otu$j	$dist" >> coverage5.tsv
done