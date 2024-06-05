# from ‘run info’ downloaded from SRA website extract SRR numbers
species=$(ls output/orthologue_promoter_FASTA | sed -e "s/_one2onepromoters_seq.fasta$//")

mkdir('output/FIMO_output')

for run in $species

do

mkdir(output/FIMO_output/"$run")

fasta-get-markov orthologue_promoter_FASTA/"$run"_one2onepromoters_seq.fasta > orthologue_promoter_FASTA/"$run"_Markov.txt

fimo --oc FIMO_output/"$run" --bfile orthologue_promoter_FASTA/"$run"_Markov.txt input/all_MEME_motifs.txt orthologue_promoter_FASTA/"$run"_one2onepromoters_seq.fasta

done