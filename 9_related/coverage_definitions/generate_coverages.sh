echo -n > coverage4.tsv
for i in $(seq 1 60);
do
    echo "Otu$i	3" >> coverage4.tsv
done