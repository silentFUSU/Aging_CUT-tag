input=$1
output=$2
echo input file: $input
echo output file: $output
awk -v OFS='\t' '{ if (NF >3) {print $4,$1,$2,$3,$6 } 
  else { print $1":"$2"-"$3,$1,$2,$3,"." } }' $1 > $2