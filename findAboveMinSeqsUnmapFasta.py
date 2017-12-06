#!/usr/bin/env python

#script to extract reads that are present and not mapped to a target genome in a minimum number of times over a group of samples from sam files
#Anna Heintz-Buschart, May 2017

import os
import sys
import argparse
import re

parser = argparse.ArgumentParser(description='Process samfiles to find reads with minimum sum of reads.')
parser.add_argument('inputFiles', help='samfiles to filter',nargs='+')
parser.add_argument('-s','--minSumReads', default=1,type=int,help='minimum number of sum of reads in all samples')
parser.add_argument('-o','--outputFile', default="validUnmappedSequences.fa",type=str,help='name of output fasta file')
parser.add_argument('-n','--commonName', default="",type=str,help='part of name of all input files to be removed from sample names for fasta file')

args = parser.parse_args()
minsum = int(args.minSumReads)
sams = args.inputFiles
outfastafile = args.outputFile
removeSample = args.commonName

dictval = {}
dictseen = {}

out_fasta_file =  open(outfastafile,"w")

for sample in sams:
    basesample = os.path.basename(sample)
    shortsample = basesample.replace(removeSample,"",1)
    outfile = basesample[:-4] + ".validCounts.tsv"
    dictsam = {}
    in_file =  open(sample,"r")
    out_file =  open(outfile,"w")
    while 1:
        lines = in_file.readline()
        if lines == "":
            break
        lines = lines.rstrip()
        if re.match("@",lines) is None:
            tabs = lines.split("\t")
            if int(tabs[1]) == 4:
                name = tabs[0]
                if name not in dictsam:
                    dictsam[name] = 1
                    seq = tabs[9]
                    countc = name.split("-")[1]
                    if seq in dictval:
                        out_file.write(name + "\t" + seq + "\t" + countc + "\n")
                        out_fasta_file.write(">" + shortsample + "." + name + "\n")
                        out_fasta_file.write(seq + "\n")
                    else:
                        counti = int(countc)
                        if counti >= minsum:
                            dictval[seq] = 1
                            out_file.write(name + "\t" + seq + "\t" + countc + "\n")
                            out_fasta_file.write(">" + shortsample + "." + name + "\n")
                            out_fasta_file.write(seq + "\n")
                        else:
                            if seq in dictseen:
                                newtot = dictseen[seq]["totalcount"] + counti
                                if newtot >= minsum:
                                    out_file.write(name + "\t" + seq + "\t" + countc + "\n")
                                    out_fasta_file.write(">" + shortsample + "." + name + "\n")
                                    out_fasta_file.write(seq + "\n")
                                    for i in range(len(dictseen[seq]["values"])):
                                        curr_sample = dictseen[seq]["values"][i][0]
                                        curr_outfile = curr_sample[:-4] + ".validCounts.tsv"
                                        curr_out_file =  open(curr_outfile,"a")
                                        curr_out_file.write(dictseen[seq]["values"][i][1] + "\t" + seq + "\t" + dictseen[seq]["values"][i][2] + "\n")
                                        curr_out_file.close()
                                        curr_shortsample = curr_sample.replace(removeSample,"",1)
                                        out_fasta_file.write(">" + curr_shortsample + "." + dictseen[seq]["values"][i][1] + "\n")
                                        out_fasta_file.write(seq + "\n")
                                    dictval[seq] = 1
                                    del dictseen[seq]
                                else:
                                    dictseen[seq]["totalcount"] = newtot
                                    dictseen[seq]["values"].append([basesample,name,countc])
                            else:
                                dictseen[seq] = {}
                                dictseen[seq]["totalcount"] = counti
                                dictseen[seq]["values"] = [[basesample,name,countc]]
    in_file.close()
    out_file.close()
    

