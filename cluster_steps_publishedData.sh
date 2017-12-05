#! /bin/bash -l

#contains the steps to extract information from cd-hit-est-2d output to R workspace
# note: you need to set the path to CD-HIT

filter="python findAbove1PerSampleUnmapFasta.py ./*/aln_hg38/*hg38.sam -s 2 -n .sam"
echo $filter
eval $filter

cdhit=/path_to_cdhit/cd-hit-est-2d #absolutePaths
CDHIT_BASE=seqsOI.fa
CDHIT_IN=gt1psUnmappedSequences.fa
CDHIT_OUT=gt1psSequences.unmapped.100.clustered.14.seqIO
CDHIT_REPORT=cdhit2D.unmapped.100.14.log

time $cdhit -i ${CDHIT_BASE} -i2 ${CDHIT_IN} -o ${CDHIT_OUT} -c 1 -n 8 -G 0 -A 14 -S2 40 -g 1 -r 0 -M 40960 -d 0 -T 12 -bak 1 > ${CDHIT_REPORT}

module load lang/R/3.0.2-ictce-5.3.0

Rscript evaluate_seqOICounts.R





