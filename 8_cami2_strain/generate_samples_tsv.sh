echo -e 'names\treads1\treads2\tcoverages\ttruths' > samples.tsv
for i in $(seq 0 99); do
  samp="strain$i"
  echo -e "$samp\tlocal_reads/$samp.1.fq.gz\tlocal_reads/$samp.2.fq.gz\tcoverage_definitions/$samp.tsv\ttruths/$samp.condensed" >> samples.tsv
done