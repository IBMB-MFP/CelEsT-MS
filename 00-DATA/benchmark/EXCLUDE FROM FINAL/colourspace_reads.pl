my $input_file = $ARGV[0]

open my $in, "<"

while(<$in>){chomp;

my @line = split("\t");

my @controlsamples = $line[9];
my @treatsamples = $line[10];

@control_array = split(', ', $controlsamples);
@treat_array = split(', ', $treatsamples);

@all_array = (@control_array, @treat_array);

print "@all_array\n";

for(@all_array){

print "$_\n";

`prefetch "$_"`;
`fastq-dump --outdir sra_fastq/"$_"`;

`bowtie 

}