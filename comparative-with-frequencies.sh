#! /bin/bash

#Bash script to perform all-against-all BLASTP-based comparative analysis
#BLASTP and MCLBLASTLINE need to be installed and in the path

namescript="${0##*/}"

display_usage(){

cat <<EOF


Synopsis

    $namescript [-h | --help] input dbname cores

Description:

    Bash script performing all-against-all BLASTP-based comparative analysis using the NCBI+ BLASTP and mclblastine.
    As output the script returns 
    1) Blastp file
    2) Raw dump.out mcl file
    3) Mclfile_with_families.txt (dump.out file containing a family ID as first column).
    4) Frequency.txt matrix (gene frequecies per family)
    5) Mclmatrix.txt (presence/absence binary matrix)

    Additional files named locustags.txt and annotation.txt need to be provided if you want to create a frequency matrix, a presence/absence binary matrix and annotated mclmatrix.

Parameters:

    input_file
        Input file contains the concatenated multifasta aminoacidic sequences of the CDS in your genomes

    dbname
        Name of the blast database

    cores
        Number of cores

    -h, --help
        Brings up this help page

Author
    Francesca Bottacini, Munster Technological University.
    francesca.bottacini@mtu.ie


EOF
}


input="$1"
dbname="$2"
cores="$3"


#Less than two arguments supplied, display usage 
	if [  $# -le 2 ]; 
	then 
		display_usage
		exit 1
	fi 
 
# User had supplied -h or --help . If yes display usage 
	if [[ $# == "--help" || $# == "-h" ]]; 
	then 
		display_usage
		exit 0
	fi 


echo "All good to go!"


rm -fr .blast/ ;
mkdir .blast/ ;
cp $input input.fas ; mv input.fas .blast/ ;
cd .blast/ ;


#Create BLAST database
echo "Creating BLAST database..."
formatdb -i input.fas -n $dbname -p T -o T ;


#BLAST all-vs-all
echo "Now blasting..."
#blastall -i $input.fas -d $dbname -p blastp -e 0.0001 -m 8 -a $cores -o $input.blastp ;
blastp -query input.fas -out $dbname.blastp -task blastp -db $dbname -evalue 0.0001 -num_threads $cores -outfmt 6 ;


#MCLBLASTLINE
echo "Running MCL clustering"
mclblastline --blast-m9 --blast-score=b --blast-sort=a --blast-bcut=5 --mcl-I=2.5 $dbname.blastp ;

cd .. ;

cp .blast/*.blastp . ;
cp .blast/dump.out* . ;

rm -fr .blast/ ;


#MCLPARSE compute frequencies
cp dump.out.$dbname.blastp.I25 mclfile ;

sed -i 's/^/Fam\t/g' mclfile ;
awk 'BEGIN {rval=1} /Fam/ {sub("Fam", "Fam"rval++)}; 1' mclfile > mclfile_with_families.txt ;

# separate the gene ids from the genome ids
sed -e 's/_/\t/g' mclfile_with_families.txt > mclfile.tmp ;

echo "Now building frequencies.txt and mclmatrix.txt..."

#Generate an array of genome ids
(genids=`cat locustags.txt | tr -s "\n" " "; echo ""`

for genid in ${genids[@]}; do
echo -en "\t$genid";
done
echo "";

while read line ; do
fams=$(echo $line|cut -d " " -f 1 );
echo -en "$fams\t";
for genid in ${genids[@]}; do
freq=`echo $line | awk -F "$genid" '{print NF-1}'`;
echo -ne "$freq\t";
done
echo ""; done<mclfile.tmp) > frequencies.txt

#Cleanup
rm mclfile.tmp;

#MCLPARSE presence/absence matrix
awk 'NR>1 {for (i=2; i<=NF; i++) {$i = ($i > 1 ? 1 : $i)}}1' frequencies.txt > mclmatrix.txt ;
sed -i "s/ /\t/g" mclmatrix.txt ; 

#Annotate mclmatrix
#extract the first two columns from mclmatrix_with_families.txt
cut -f1,2 mclfile_with_families.txt > mclfile_reference.txt ;

#add annotation
python <<EOF

import pandas as pd

df1 = pd.read_csv('mclfile_reference.txt', delimiter="\t", header = None)
df1.columns = ["Family", "Locus"]

df2 = pd.read_csv('annotation.txt', delimiter="\t", header = None)
df2.columns = ["Locus", "Annotation"]

df3 = pd.merge(df1, df2, on='Locus', how='left')

df3.to_csv('mclfile_reference_annotated.csv', sep='\t', index=False, header = None)
df3.to_csv('mclfile_reference_annotated_head.csv', sep='\t', index=False)

#pd.exit(0)

EOF

paste mclfile_reference_annotated.csv mclfile_with_families.txt > mclfile_with_families_annot.txt ;
paste mclfile_reference_annotated_head.csv mclmatrix.txt > mclmatrix_annot.txt ;

rm mclfile_reference_annotated.csv mclfile_reference_annotated_head.csv mclfile_reference.txt

exit 0
