 results/WGD2ANC0*/*comp | xargs wc -l | awk '{if(!/total/ && $1 > 0) print $2}' | xargs awk FNR-1 > results/15may14consensuscomp.txt
