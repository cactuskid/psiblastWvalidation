#! /bin/bash
for file in  ./input/*.fasta; do
	psiblast -query $file -out $file.csv -outfmt ' 10 qseqid qlen slen qcovhsp sseqid staxids bitscore score evalue pident qstart qend sstart send' -db /scratch/cluster/monthly/dmoi/metaclust/metaclust50_2017_01.fasta  -num_iterations 3 -num_threads 50
done
