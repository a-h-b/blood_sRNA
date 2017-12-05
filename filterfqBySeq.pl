#!/usr/bin/perl
# script to filter fastq by sequence (exact match)
# Anna Heintz-Buschart

use strict;
use Getopt::Long;

my ($sequenceFile,$fastqFile,$exclude,$output);
GetOptions('s=s' => \$sequenceFile,
	   'f=s' => \$fastqFile,
        'e=s' => \$exclude,
        'o=s' => \$output);
$exclude = $exclude ? 0 : 1;

#read fasta with sequences to be kept (generic so it could be normal text too)
my %reads = ();
open (IN, $sequenceFile);
while (my $line = <IN>){
    unless ($line =~ /^>/) { 
        chomp $line;
		$reads{$line} = 1;
		}
}
close(IN);

#read fastq
#keep as l1, l2, l3, l4 (count lines)
#compare l2 with the read sequences
#- if present, write l1-l4 to file
my $linecount = 0;
my @ls = ("","","","","");
open (IN, $fastqFile);
open(FASTQout, ">", $output) or die "cannot open $output \n";
while (my $lines = <IN>){
        chomp $lines;
		if($linecount == 4){
			$linecount = 1;
		}else{
			$linecount += 1;
		}
		$ls[$linecount-1] = $lines;
		if($linecount==4){
			if($exclude == 1){
				if(exists($reads{$ls[1]})) {
					print FASTQout $ls[0], "\n";
					print FASTQout $ls[1], "\n";
					print FASTQout $ls[2], "\n";
					print FASTQout $ls[3], "\n";
				}
			}else{
				unless(exists($reads{$ls[1]})) {
					print FASTQout $ls[0], "\n";
					print FASTQout $ls[1], "\n";
					print FASTQout $ls[2], "\n";
					print FASTQout $ls[3], "\n";
				}
			}
		}
}
close(IN);
close(FASTQout);
