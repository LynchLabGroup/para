ls WGD2ANC0*/*motifs | xargs wc -l | awk '{sub(/.WGD2ANC[0-9]*.motifs/,"")}1 && !/total/' >> motifs.lengths
