#!/bin/sh

#To run a large number of clustal alignments.
#Inputs are a text file with a list of multifasta groups from orthofinder to align

echo "Script usage supertree_builder.sh groupfile"


for f in *.fa ; do clustalw -INFILE="$f" -ALIGN -TYPE=PROTEIN -MATRIX=BLOSUM -ITERATION=ALIGNMENT -OUTFILE="$f".tmpphy -OUTPUT=PHYLIP ; done ;

for f in *.fa ; do /usr/local/share/scripts/conversions/phylip2fasta.py "$f".tmpphy "$f".out ; done ;

for f in *.fa ; do Gblocks "$f".out -t=p -e=gbl -b4=5 ; done ;

for f in *.fa ; do /usr/local/share/scripts/conversions/fasta2phylip.py "$f".outgbl "$f".phy ; done ;

for f in *.fa ; do phyml -i "$f".phy -d 'aa' -b 1 -m JTT  ; done ;


#Cleanup
rm *.tmpphy *.outgbl ;

exit 0
