FILEVEC=$(cat input/ENCODE_extraTFs_files.txt)

mkdir input/ENCODE_extraTFs

for FILE in $FILEVEC

do

curl -0 -J -L https://www.encodeproject.org/files/$FILE/@@download/$FILE.bed.gz --output input/ENCODE_extraTFs/$FILE.bed.gz

done 
