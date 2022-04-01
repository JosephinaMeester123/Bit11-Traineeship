#!/usr/bin/env bash

#make adjustments next time: loop over files in a folder, which is easier if you don't know how many samples are in the folder
paste /home/jmeester/WGS_Jeannette/Testing/*.pbs.biomina.be | awk '{print $1, $2, $3, $4, $8, $12, $16, $20, $24, $28, $32, $36}' > /home/jmeester/Internship/Reference_test/count_allsamples.txt

awk '{T=0; V=0;
	for(N=4; N<=NF; N++) T+=$N; T/=(NF-3);
	for(N=4; N<=NF; N++) V+=(($N-T)*($N-T));
	stdev = sqrt(V/(NF-3));
	print $1, $2, $3, T, stdev }' /home/jmeester/Internship/Reference_test/count_allsamples.txt > /home/jmeester/Internship/Reference_test/average-stdev.txt
