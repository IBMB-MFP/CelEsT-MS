cd ~/Cel_GRN_manuscript

samples=$(awk -F '\t' '{ print $1 }' input/benchmarkpairedsamples.txt)

for run in $samples

do

echo "$run"

prefetch "$run"
fastq-dump --outdir sra_fastq/ --split-files "$run"

bowtie2 -x bowtie_indexes/CelWS288 -1 sra_fastq/"$run"_1.fastq -2 sra_fastq/"$run"_2.fastq -p 8 > temp_paired.sam

featureCounts -p -t exon -g gene_id -a c_elegans.PRJNA13758.WS288.canonical_geneset.gtf -o input/benchmark_counts/"$run"_FCcounts.txt temp_paired.sam

done
