cd ~/Cel_GRN_manuscript

motiflist=$(ls output/uniprobe_conversions)

cd output/uniprobe_conversions

uniprobe2meme $motiflist > all_MEME_motifs.txt

mv all_MEME_motifs.txt ..

cd ..

fasta-get-markov Cel_promoters_censored.fasta > Cel_promoters_censored_Markov.txt

fimo --bfile Cel_promoters_censored_Markov.txt all_MEME_motifs.txt Cel_promoters_censored.fasta