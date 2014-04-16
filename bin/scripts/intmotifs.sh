#!/bin/bash
# Script to look at how many interesting motifs have been found
# USAGE:
# intmotifs file_to_write

temp=$(for FAM in `ls WGD2ANC0*/*motifs | xargs wc -l | awk '{sub(/.WGD2ANC[0-9]*.motifs/,"")}1 && !/total/' | awk '{if($1 > 2) print $2}'`;
do
	echo $FAM
	awk '{if($11 != "" && $11 != 0.0) print $3"\t"$11}' $FAM/$FAM.motifs
done)

echo "$temp" | head
OIFS="$IFS"
i=0
echo "Writing file $1..."
while IFS=\t read -r LINE
do	
	if [ $i == 0 ]
	then
		i=1
		echo -e "family\tdistance\tsize"
	fi

	if [[ $LINE == *WGD* ]]; then
		NAME="$LINE"
	else
		IFS=" " read -ra a <<< $LINE
		echo -e $NAME"\t"${a[0]}"\t"${a[1]}
	fi
done <<< "$temp" > $1
IFS=$OIFS
echo "Done."
