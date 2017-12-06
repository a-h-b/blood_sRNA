#! /bin/bash -l

filter="python findAboveMinSeqsUnmapFasta.py ../analysis/*/aln_hg38/*hg38.sam -s 30 -n .sam"
echo $filter
eval $filter

cdhit=/path_to_CDHIT/cd-hit-est #absolutePaths
CDHIT_IN=validUnmappedSequences.fa
CDHIT_OUT=validUnmappedSequences.100.100.clustered
CDHIT_REPORT=cdhit.100.100.validUnmappedSequences.log

time $cdhit -i ${CDHIT_IN} -o ${CDHIT_OUT} -c 1 -n 10 -G 0 -aL 1 -aS 1 -uL 0 -uS 0 -S 0 -g 1 -r 0 -M 40960 -d 0 -T 12 -bak 1 > ${CDHIT_REPORT}

#module load lang/R/3.0.2-ictce-5.3.0 #for university of Luxembourg's GAIA cluster
#module load module load lang/R/3.4.0-foss-2017a-X11-20170314-bare #new OS

Rscript firstClustering.unmapped.R

CDHIT_IN2=validUnmappedSequences.100.100.clustered
CDHIT_OUT2=validUnmappedSequences.100.100.clustered.100.oh
CDHIT_REPORT2=cdhit.100.oh.validUnmappedSequences.log

time $cdhit -i ${CDHIT_IN2} -o ${CDHIT_OUT2} -c 1 -n 10 -G 0 -aL 0.3 -aS 1 -uL 0 -uS 0 -S 34 -g 1 -r 0 -M 40960 -d 0 -T 12 -bak 1 > ${CDHIT_REPORT2}

Rscript bigClusters.unmapped.R

