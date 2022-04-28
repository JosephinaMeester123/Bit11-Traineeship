#!/usr/bin/env bash

#for each sample, extract the 4th column, which contains the counts
for filename in /home/jmeester/WGS_Jeannette/Testing/*.count.stdout.*; do
	awk '{print $4} ' $filename > $filename".count"
done

#paste together windows1kb file, which contains locations of bins, with the counts of all samples
paste /home/jmeester/Internship/windows1kb.bed /home/jmeester/WGS_Jeannette/Testing/*.count > /home/jmeester/Internship/Reference_test/count_allsamples.txt

awk 'BEGIN{OFS=","} {T=0; V=0;
	for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
	for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
	stdev = sqrt(V/(NF-3));
	print $1, $2, $3, T, stdev }' /home/jmeester/Internship/Reference_test/count_allsamples.txt > /home/jmeester/Internship/Reference_test/referenceNormals.stats.binned.csv
