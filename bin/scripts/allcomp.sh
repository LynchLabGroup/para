# Allcomp.sh
# Produces a list with all motifs in .comp files
# USAGE:
# bin/allcomp.sh output
output=$1
results/WGD2ANC0*/*comp | xargs wc -l | awk '{if(!/total/ && $1 > 0) print $2}' | xargs awk FNR-1 > $output
