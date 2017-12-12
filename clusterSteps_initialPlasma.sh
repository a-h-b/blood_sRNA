#! /bin/bash -l

filter="python /home/users/aheintzbuschart/myScripts/findAboveMinSeqsUnmapFasta.py ./*/aln_hg38/*hg38.sam -s 10 -n .sam"
echo $filter
eval $filter

cdhit=/work/projects/ecosystem_biology/local_tools/cd-hit-v4.6.1-2012-08-27_OpenMP/cd-hit-est
CDHIT_IN=validUnmappedSequences.fa
CDHIT_OUT=validUnmappedSequences.100.100.clustered
CDHIT_REPORT=cdhit.100.100.validUnmappedSequences.log

time $cdhit -i ${CDHIT_IN} -o ${CDHIT_OUT} -c 1 -n 10 -G 0 -aL 1 -aS 1 -uL 0 -uS 0 -S 0 -g 1 -r 0 -M 40960 -d 0 -T 12 -bak 1 > ${CDHIT_REPORT}

#module load lang/R/3.0.2-ictce-5.3.0

module load lang/R/3.4.0-foss-2017a-X11-20170314-bare
#170811_clusteredCounts.100.100.unmapped.R
Rscript 171108_clusteredCounts.100.100.unmapped.R

#CDHIT_IN=validUnmappedSequences.100.100.clustered
#CDHIT_OUT=validUnmappedSequences.100.100.clustered.100.oh
#CDHIT_REPORT=cdhit.100.oh.validUnmappedSequences.log

#time $cdhit -i ${CDHIT_IN} -o ${CDHIT_OUT} -c 1 -n 10 -G 0 -aL 0.3 -aS 1 -uL 0 -uS 0 -S 34 -g 1 -r 0 -M 40960 -d 0 -T 12 -bak 1 > ${CDHIT_REPORT}

#module load module load lang/R/3.4.0-foss-2017a-X11-20170314-bare

#Rscript 171006_clusteredCounts.100.100.unmappedR.bigClusters.R

#/work/projects/ecosystem_biology/Qiagen/plasmaCombined$ for i in $(cat 171018_potentialExogenousReads.txt); do grep ^$i -w validUnmappedSequences.100.100.clustered.bak.clstr | grep -v %$ | cut -f 2 -d ">" | sed 's/[.][.][.]//g' | cut -f 1 -d " " | grep -f - -w -A 1 validUnmappedSequences.fa >> 171018_potentialExogenousReads.fa; done

#for i in $(cat  171018_potentialExogenousReads.noHuman.noSalter.contignames.txt); do grep $i validUnmappedSequences.100.100.clustered.bak.clstr | cut -f 1 >> 171018_potentialExogenousReads.noHuman.noSalter.clusters; done

#for lib in $(cat ../analysis/ids); do grep "^@" -v ../analysis/$lib/aln_salmonella/${lib}_trim.fa_salmonella.sam >> $lib.small.sam; grep "^@" -v ../analysis/$lib/aln_hg38/${lib}_trim.fa_hg38.sam >> $lib.small.sam; done
#stepthree="python /home/users/aheintzbuschart/myScripts/161130_qiagen_clusterSam.py *.small.sam -c validSequences.clustered.clstr -n .small.sam" 
#to quick-check the clusters of interest:
#for clus in `cut -f 2 -d " " clustersOfInterest`; do grep ^${clus} mappersPerCluster -w | cut -f 1,4,5,10,11,16,17; done


