## This repository documents the analysis steps to identify potential or confirmed contaminants in small RNA sequencing data 
The results of these analyses form part of a manuscript, which you can find [here](https://www.biorxiv.org/). The links in the sub-headings below lead to short descriptions of the workflows. Below each subheading are lists of the scripts which are described in more detail in the linked description and the actual data.

### [Analysis steps of the central dataset](main_steps.md) (control samples and the plasma titration experiment):
This is the central part of the project. The steps are related to per-sample read processing, comparison of samples by forming cross-sample read clusters and visualization of results.
  * [runAnalysis_smallRNA_noSalmonella.sh](runAnalysis_smallRNA_noSalmonella.sh)
  * [runAnalysis_smallRNA_withSalmonella.sh](runAnalysis_smallRNA_withSalmonella.sh)
  * [clean_lib.py](clean_lib.py)
  * [main_cluster_steps.sh](main_cluster_steps.sh)
  * [findAboveMinSeqsUnmapFasta.py](findAboveMinSeqsUnmapFasta.py)
  * [firstClustering.unmapped.R](firstClustering.unmapped.R)
  * [bigClusters.unmapped.R](bigClusters.unmapped.R)
  * [clusteringToSeqOI.unmapped.R](clusteringToSeqOI.unmapped.R)
  * [clusteringToSeqOI.byMapping.R](clusteringToSeqOI.byMapping.R)
  * [titration4plot.wostop.txt](titration4plot.wostop.txt)
  * [readCountsPerSeqOIPerSample.unmapped.14.100forplot]()
  * [Figure5_titration.R](Figure5_titration.R)
  * [unmappedReadCountsPerCluster100.100PerSample_withBiggerCluster.RDS](unmappedReadCountsPerCluster100.100PerSample_withBiggerCluster.RDS)
  * [stops.sorted.txt](stops.sorted.txt)
  * [algs.sorted.txt](algs.sorted.txt)
  * [8182s.sorted.txt](8182s.sorted.txt)
  * [3869s.sorted.txt](3869s.sorted.txt)
  * [3378s.sorted.txt](3378s.sorted.txt)
  * [8030s.sorted.txt](8030s.sorted.txt)
  * [6848s.sorted.txt](6848s.sorted.txt)
  * [anna3tab2venn.R](anna3tab2venn.R)
  * [figure4B-D_potentialContaminants_Exogenous.R](figure4B-D_potentialContaminants_Exogenous.R)

### [Analysis of published](published.md) sRNA sequencing datasets:
This part documents the steps taken to detect the confirmed contaminants in published datasets and the visualization of the results.
  * [runAnalysis_smallRNA_noSalmonella.sh](runAnalysis_smallRNA_noSalmonella.sh)
  * [cluster_steps_publishedData.sh](cluster_steps_publishedData.sh)
  * [findAbove1PerSampleUnmapFasta.py](findAbove1PerSampleUnmapFasta.py)
  * [evaluate_seqOICounts.R](evaluate_seqOICounts.R)
  * [allLibsMeta.txt](allLibsMeta.txt)
  * [allReads.tsv](allReads.tsv)
  * [public.readCountsPerSeqOIPerSample.unmapped.14.t.tsv](public.readCountsPerSeqOIPerSample.unmapped.14.t.tsv)
  * [heatmap.2m.R](heatmap.2m.R)
  * [figure3_publishedData.R](figure3_publishedData.R)

### [Miscellaneous](miscellaneous.md):
qPCR data and the scripts used to plot the figures.
  * [figure2_WS.Rdata](figure2_WS.Rdata)
  * [figure2_qPCR-results.R](figure2_qPCR-results.R)
  * [161013_qiaSummary_withDnase_extractUCTvsOld.txt](161013_qiaSummary_withDnase_extractUCTvsOld.txt)
  * [figure4A_qPCR-results.R](figure4A_qPCR-results.R)
  * [matchcountSeqOIPerSample4plot.txt](matchcountSeqOIPerSample4plot.txt)
  * [supplementaryFigure3_reducedContaminants.R](supplementaryFigure3_reducedContaminants.R)
  
