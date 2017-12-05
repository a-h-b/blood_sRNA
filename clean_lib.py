#! /usr/bin/env python
"""
script for cleaning  adapter/primer seqeunces from short RNA reads
"""
#contributed by Dilmurat Yusuf

def usage():
    print
    print 'script for cleaning the adapter/primer seqeunces from short RNA reads'
    print 'Version: 0.1'
    print
    print 'Command example:'
    print 'clean adaptor: clean_lib.py -a adapter fastq_dir/*fastq.gz'
    print 'run fastqc: clean_lib.py -q fastq_dir/*fastq'
    print 'generate table  of primer contamination: clean_lib.py -t fastqc_dir/*/fastqc_data.txt'
    print 'clean  primers: clean_lib.py -p primer_contamination.csv  fastq_dir/*fastq'
    print 'quality_filter/collapse identical reads: clean_lib.py -c fastq_dir/*fastq'
    print 
    
    sys.exit()
    
def check_cmd_opt():
    try:
        cmd_opt, input_list = getopt.getopt(sys.argv[1:],'a:qtp:c')
    except getopt.GetoptError as error:
        print '!', error
        usage()
    
    if not input_list:
        print '\n!please specify input file(s)\n'
        usage()
    
    #cmd_opt, input_list = getopt.getopt(sys.argv[1:],'c:')
        
    return cmd_opt, input_list

def omit_tag(file_name):
        file_tag = file_name.split('/')[-1].split('.')[0]
        return file_tag

def run_cmd(cmd):    
    print 'CMD: %s' % cmd
    process = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    print ''.join(process.stderr.readlines())
    print ''.join(process.stdout.readlines())

def clean_adapter(adapter, fastqs):
    # trim adapter
    # quality filter
    
    #clipper = commands.getoutput("which fastx_clipper")
    clipper = "~/bin/fastx_clipper"
    
    for fastq in fastqs:
	file_tag = omit_tag(fastq)
	ta_fq =   file_tag + '.ta.fastq'
	
	cmd = 'zcat %s | %s -v -a %s -l 11 -Q33 -c > %s' % (fastq, clipper,  adapter, ta_fq)
	
	run_cmd(cmd)

def run_qc(fastqs):
    qc = '~/bin/fastqc'
    for fastq in fastqs:
	cmd = '%s %s' % (qc, fastq)
	run_cmd(cmd)

def primer_info(fastqc):
    primers = {}
    
    for qc in fastqc:
        on = False
        for line in open(qc):
            if line[:2] == ">>":
                line = line.split()
                if line[:-1] == ['>>Overrepresented', 'sequences']:
                    on = True
                if line ==  ['>>END_MODULE']:
                    on = False
            if on:
                if type(line) != list and line[0] != "#":
                    line = line.strip()
                    line = line.split("\t")
                    
                    read = line[0]
                    percent = float(line[2])
                    primer_type = line[-1]
                    
                    if primer_type == "No Hit":
                        continue
                    
                    primer_len = primer_type.split()[-1].replace("bp)","")
                    primer_len = int(primer_len)
                    primer_seq = read[-primer_len:]
                    
                    primer_type = primer_type.replace(",","")
                    
                    if primer_type not in primers:
                        primers[primer_type] = {}
                        primers[primer_type]["seq"] = primer_seq
                    if qc not in primers[primer_type]:
                        primers[primer_type][qc] = percent
                    else:
                        primers[primer_type][qc] += percent
    return  primers
    
def primer_table(fastqc, primers):
    primer_types = sorted(primers.keys())
    primers_len =  len(primer_types)
    
    primer_contam = open("primer_contamination.csv", "w")
    
    primer_line = "primer\t"+"%s\t"*primers_len % tuple(primer_types)
    print >> primer_contam, primer_line.rstrip()
    seq_line = "primer_seq\t" + "%s\t"*primers_len % tuple([primers[primer_type]["seq"] for primer_type in primer_types])
    print >> primer_contam, seq_line.rstrip()
    
    for qc in sorted(fastqc):
        qc_tag = qc.split('/')[-2].split('.')[0]
        qc_line = "%s\t" % qc_tag
        for primer_type in primer_types:
            if qc in primers[primer_type].keys():
                percent = primers[primer_type][qc]
                percent = "%.2f" % percent
                
                primer_seq = primers[primer_type]["seq"]
                
            else:
                percent = 0
            qc_line += "%s\t" % percent
            
        print >>primer_contam, qc_line.rstrip()    
    
    print "Contamination table: primer_contamination.csv" 
	
def lib_primer(primer_csv):
        primer_seqs = {}
        
	reader = csv.reader(open(primer_csv), delimiter="\t")
	for row in reader:
		if row[0] == "primer_seq":
			primer = row
			break
	for row in reader:
		qc_tag = row[0]
		primer_indexs = [row.index(percent) for percent in row[1:] if float(percent) > 0]
                for primer_index in primer_indexs:
                        primer_seq = primer[primer_index]
                        #clean_lib(qc_tag, fastq_dir, primer_seq)
                        primer_seqs[qc_tag] = primer_seq
        
        return primer_seqs
    
def clean_primer(primer_seqs, fastqs):
        #clipper_path = commands.getoutput("which fastx_clipper")
        clipper_path = '~/bin/fastx_clipper' #_t    
        
        for fastq in fastqs:
            qc_tag = fastq.split('/')[-1].split('.')[0]
	    if qc_tag in primer_seqs:
		primer_seq = primer_seqs[qc_tag] 
		cmd = "%s -v -a %s -l 11 -Q33 -i %s -o %s.clean.fastq" % (clipper_path, primer_seq, fastq, qc_tag)
		run_cmd(cmd)
	    

def collapse(fastqs):
    #qa = commands.getoutput("which fastq_quality_filter")
    qa = "~/bin/fastq_quality_filter"
    #clap = commands.getoutput("which fastx_collapser") 
    clap = "~/bin/fastx_collapser"
    
    for fastq in fastqs:
	file_tag = omit_tag(fastq)
        pros_fa = file_tag + '.pros.fa'
	cmd = '%s -v -q 25 -p 100 -Q33  -i %s |\
	%s -v -Q33  > %s' % (qa, fastq, clap, pros_fa)
        run_cmd(cmd)
	
if __name__ == '__main__':
    import sys, getopt
    import subprocess, commands, time
    import csv
    
    cmd_opt, infile = check_cmd_opt()
    if cmd_opt[0][0] == '-a':
	fastqs = infile
	adapter = cmd_opt[0][1]
	clean_adapter(adapter, fastqs)
    
    elif cmd_opt[0][0] == '-q':
	fastqs = infile
	run_qc(fastqs)
	
    elif cmd_opt[0][0]  == '-t':
        fastqc_dir = infile
        primers =  primer_info(infile)
        primer_table(infile, primers)
        
    elif cmd_opt[0][0] == '-p':
        primer_csv = cmd_opt[0][-1]
        fastqs = infile
        primer_seqs = lib_primer(primer_csv)
        clean_primer(primer_seqs, fastqs)
    
    elif cmd_opt[0][0] == '-c':
	fastqs = infile
	collapse(fastqs)
        
        
        
        
