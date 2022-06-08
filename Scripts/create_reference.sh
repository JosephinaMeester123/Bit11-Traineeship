#!/usr/bin/env bash

#for each sample, extract the 4th column, which contains the counts

for filename in /home/jmeester/Internship/countfiles/101557.count; do
       awk '$1 != "chrY" && $1 != "chrX"' $filename > $filename".autosomal"
       awk '{print $1, $2, $3} ' $filename".autosomal" > /home/jmeester/Internship/countfiles/locations.txt
done

echo "filename	autosomal_total_count	average_count	stdev_count" > /home/jmeester/Internship/countfiles/summarycounts.txt
for filename in /home/jmeester/Internship/countfiles/*.count; do
	awk '$1 != "chrY" && $1 != "chrX"' $filename > $filename".autosomal"
	awk '{print $1, $2, $3} ' $filename".autosomal" > /home/jmeester/Internship/countfiles/locations.txt
        awk '{print $4} ' $filename | awk '$1 != 0' | awk -v myvar="$filename" '{ total += $1; sumsq+=$1*$1 } END { print myvar, total, total/(NR-1), \
		sqrt(sumsq/(NR-1) - (total/(NR-1))^2)}' >> /home/jmeester/Internship/countfiles/summarycounts.txt
done

#perform normalization
totalavg="$(awk '{ total += $2} END { print total/(NR-1)}' /home/jmeester/Internship/countfiles/summarycounts.txt)"
echo $totalavg 

for file in /home/jmeester/Internship/countfiles/*.count.autosomal; do
	totalsample="$(awk '{total += $4} END {print total}' $file)"
	awk -v totalavg="$totalavg" -v totalsample="$totalsample" '{print $1, $2, $3, $4, totalsample, totalavg, $4/totalsample*totalavg}' $file > $file".norm"
	awk '{print $7} ' $file".norm" > $file".norm.countonly"
done

#paste together windows500bp file, which contains locations of bins, with the counts of all samples
paste /home/jmeester/Internship/windows500bp.bed /home/jmeester/Internship/countfiles/*.norm.countonly > /home/jmeester/Internship/Reference_test/count_allsamples_norm.txt

#calculate average and standard deviation per row to create a reference set
awk 'BEGIN{OFS=","} {T=0; V=0;
	for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
	for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
	stdev = sqrt(V/(NF-3));
	print $1, $2, $3, T, stdev }' /home/jmeester/Internship/Reference_test/count_allsamples_norm.txt > /home/jmeester/Internship/Reference_test/referenceNormals_norm.csv

# rm /home/jmeester/Internship/countfiles/*.countonly


