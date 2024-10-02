samples=$(awk -F '\t' '{ print $1 }' benchmarksinglesamples.txt)

for run in $samples

do

echo "$run"

prefetch "$run"
fastq-dump --outdir sra_fastq/ "$run"

hisat2 -x celWS288_ht/genome_tran -U sra_fastq/"$run".fastq -p 8 > temp_single.sam

featureCounts -p -t exon -g gene_id -a c_elegans.PRJNA13758.WS288.canonical_geneset.gtf -o benchmark_counts/"$run"_FCcounts.txt temp_single.sam

done
