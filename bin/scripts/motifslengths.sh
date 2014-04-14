ls WGD2ANC0*/*motifs | xargs wc -l | awk '{sub(/.WGD2ANC[0-9]*.motifs/,"")}1 && !/total/' >> motifs.lengths

# Command to look at how many families had no motifs found
# awk '{if($1 == 2) print $0}' results/lengthsmotifsfiles.txt | wc -l
