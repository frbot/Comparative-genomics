#! /bin/bash

#Bash script to perform all-against-all BLASTP-based comparative analysis
#BLASTP and MCLBLASTLINE need to be installed and in the path
#A file named locustags.txt must be provided in order to compute the gene frequency and presence/absence matrix


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

    An additional file named locustags.txt needs to be provided if you want to create a frequency matrix and a presence/absence binary matrix.

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
    Francesca Bottacini, APC BioIT platform, University College Cork.
    f.bottacini@umail.ucc.ie


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

exit 0
