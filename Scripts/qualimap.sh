#!/bin/bash

list=$(find /home/jmeester/Internship/Qualimap/ -name genome_results.txt -print )

originalname="_genome_results.txt"
sepa="/"
bamfolder="_BamQC_20210305_qualimap"
wildcard="*"
for filename in $list; do
	sample=$(echo $filename | cut -d'/' -f6)
	location=$(echo $filename | cut -d'/' -f -7)
	oldlocation=$(echo $filename | cut -d'/' -f -8)
	oldlocationfull=$(echo $oldlocation$sepa$wildcard)
	newlocation=$(echo $location$sepa$sample$bamfolder)
        newlocationfull=$(echo $location$sepa$sample$bamfolder$sepa)
#	echo "$newlocationfull"

#	mkdir $newlocation
#	cat $filename > "$newlocation$sepa$sample$originalname"
	mv $oldlocationfull $newlocationfull
done

