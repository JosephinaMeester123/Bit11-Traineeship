#!/usr/bin/env bash

#for each sample, extract the 4th column, which contains the counts

#echo "filename	autosomal_total_count	average_count	stdev_count" > /home/jmeester/Internship/countfiles/summarycounts.txt
#for filename in /home/jmeester/Internship/countfiles/*.count; do
#	awk '$1 != "chrY" && $1 != "chrX"' $filename > $filename".autosomal"
#        awk '{print $4} ' $filename | awk '$1 != 0' | awk -v myvar="$filename" '{ total += $1; sumsq+=$1*$1 } END { print myvar, total, total/(NR-1), \
#		sqrt(sumsq/(NR-1) - (total/(NR-1))^2)}' >> /home/jmeester/Internship/countfiles/summarycounts.txt
#done

for file in /home/jmeester/Internship/countfiles/117005.binned.csv; do
        awk -F "," '{if (NR>1) {print $2, $3, $4}} ' $file > /home/jmeester/Internship/locations_500bp.txt
done

for file in /home/jmeester/Internship/countfiles/*.binned.csv; do
	awk -F "," '{if (NR>1) {print $7}} ' $file > $file".countonly"
        awk -F "," '{if (NR>1) {print $8}} ' $file > $file".gconly"
        awk -F "," '{if (NR>1) {print $9}} ' $file > $file".gr10monly"
        awk -F "," '{if (NR>1) {print $10}} ' $file > $file".gcgr10monly"
        awk -F "," '{if (NR>1) {print $11}} ' $file > $file".autogr10monly"
        awk -F "," '{if (NR>1) {print $12}} ' $file > $file".gcautogr10monly"
done

#paste together windows500bp file, which contains locations of bins, with the counts of all samples
paste /home/jmeester/Internship/locations_500bp.txt /home/jmeester/Internship/countfiles/*.countonly > /home/jmeester/Internship/Reference_test/acountonly_allsamples.txt
paste /home/jmeester/Internship/locations_500bp.txt /home/jmeester/Internship/countfiles/*.gconly > /home/jmeester/Internship/Reference_test/bgconly_allsamples.txt
paste /home/jmeester/Internship/locations_500bp.txt /home/jmeester/Internship/countfiles/*.gr10monly > /home/jmeester/Internship/Reference_test/cgr10monly_allsamples.txt
paste /home/jmeester/Internship/locations_500bp.txt /home/jmeester/Internship/countfiles/*.gcgr10monly > /home/jmeester/Internship/Reference_test/dgcgr10monly_allsamples.txt
paste /home/jmeester/Internship/locations_500bp.txt /home/jmeester/Internship/countfiles/*.autogr10monly > /home/jmeester/Internship/Reference_test/eautogr10monly_allsamples.txt
paste /home/jmeester/Internship/locations_500bp.txt /home/jmeester/Internship/countfiles/*.gcautogr10monly > /home/jmeester/Internship/Reference_test/fgcautogr10monly_allsamples.txt

#calculate average and standard deviation per row to create a reference set
awk 'BEGIN{OFS=","} {T=0; V=0;
	for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
	for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
	stdev = sqrt(V/(NF-3));
	print $1, $2, $3, T, stdev }' /home/jmeester/Internship/Reference_test/acountonly_allsamples.txt > /home/jmeester/Internship/Reference_test/acountonly_allsamples.csv

awk 'BEGIN{OFS=","} {T=0; V=0;
        for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
        for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
        stdev = sqrt(V/(NF-3));
        print T, stdev }' /home/jmeester/Internship/Reference_test//bgconly_allsamples.txt > /home/jmeester/Internship/Reference_test//bgconly_allsamples.csv

awk 'BEGIN{OFS=","} {T=0; V=0;
        for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
        for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
        stdev = sqrt(V/(NF-3));
        print T, stdev }' /home/jmeester/Internship/Reference_test/cgr10monly_allsamples.txt > /home/jmeester/Internship/Reference_test/cgr10monly_allsamples.csv

awk 'BEGIN{OFS=","} {T=0; V=0;
        for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
        for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
        stdev = sqrt(V/(NF-3));
        print T, stdev }' /home/jmeester/Internship/Reference_test/dgcgr10monly_allsamples.txt > /home/jmeester/Internship/Reference_test/dgcgr10monly_allsamples.csv

awk 'BEGIN{OFS=","} {T=0; V=0;
        for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
        for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
        stdev = sqrt(V/(NF-3));
        print T, stdev }' /home/jmeester/Internship/Reference_test/eautogr10monly_allsamples.txt > /home/jmeester/Internship/Reference_test/eautogr10monly_allsamples.csv

awk 'BEGIN{OFS=","} {T=0; V=0;
        for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
        for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
        stdev = sqrt(V/(NF-3));
        print T, stdev }' /home/jmeester/Internship/Reference_test/fgcautogr10monly_allsamples.txt > /home/jmeester/Internship/Reference_test/fgcautogr10monly_allsamples.csv

paste -d "," /home/jmeester/Internship/Reference_test/*_allsamples.csv > /home/jmeester/Internship/Reference_test/referenceNormals.csv
# rm /home/jmeester/Internship/countfiles/*.countonly


