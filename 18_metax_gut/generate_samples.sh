echo -e "names\treads1\treads2\ttruths\tcoverages"
seq 0 9 | sed 's/\(.*\)/sample_\1\tsplit_reads\/sample_\1.1.fq.gz\tsplit_reads\/sample_\1.2.fq.gz\ttruths\/sample_\1.tsv\tcoverages\/sample_\1.tsv/'