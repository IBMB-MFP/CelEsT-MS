samples=$(awk -F '\t' '{ print $1 }' ../../../input/benchmarksinglesamples.txt)

for run in $samples

do

echo "$run"

prefetch "$run"
fastq-dump --outdir sra_fastq/ "$run"

bowtie2 -x bowtie_indexes/CelWS288 -U sra_fastq/"$run".fastq -p 8 > temp_single.sam

featureCounts -t exon -g gene_id -a c_elegans.PRJNA13758.WS288.canonical_geneset.gtf -o ../../../input/benchmark_counts/"$run"_FCcounts.txt temp_single.sam

done
