samples=$(awk -F '\t' '{ print $1 }' benchmarkpairedsamples.txt)

for run in $samples

do

echo "$run"

prefetch "$run"
fastq-dump --outdir sra_fastq/ --split-files "$run"

hisat2 -x celWS288_ht/genome_tran -1 sra_fastq/"$run"_1.fastq -2 sra_fastq/"$run"_2.fastq -p 8 > temp_paired.sam

featureCounts -p -t exon -g gene_id -a c_elegans.PRJNA13758.WS288.canonical_geneset.gtf -o benchmark_counts/"$run"_FCcounts.txt temp_paired.sam

done
