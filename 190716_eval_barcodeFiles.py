#!/usr/bin/python

import os
import argparse
import re

parser = argparse.ArgumentParser(description='Accession file and eval path')
parser.add_argument('-d','--barcodedirectory', type=str,help='directory where all results of searching sequences in front of primers are kept')
parser.add_argument('-o','--out', type=str,help='output')

args = parser.parse_args()
PDIR = args.barcodedirectory
OUT = args.out

fwd="GTGYCAGCMGCCGCGGTAA"
rev="GGACTACNVGGGTWTCTAAT"
fwd_rc="TTACCGCGGCKGCTGRCAC"
rev_rc="ATTAGAWACCCBNGTAGTCC"
primers = [fwd, rev, fwd_rc, rev_rc]

print "Reading fastqs from ", PDIR
fastqs = [f for f in os.listdir(PDIR) if re.match(r'.*.fastq', f)]

all = {}
for f in fastqs:
	lcnt = 0
	tmpfile = open(PDIR+"/"+f,"r")
	barcodes = []
	while 1:
		line = tmpfile.readline()
		#print line
		if line == "":
			break
		if lcnt % 4 == 1:
			line = line.rstrip()
			barcodes.append(line)
		lcnt += 1
	tmpfile.close()
	run = re.sub('_.+fastq', '',f)
	if len(barcodes) > 0:
		if run in all:
			all[run] +=barcodes 
		else:
			all[run] = barcodes

out_file = open(OUT, "w")	
out_file.write("\t".join(["run","barcoded fragments","different barcodes"]) + "\n")
for key in all:
	tmplen = len(all[key])
	tmpuni = len(set(all[key]))
	out_file.write("\t".join([key,str(tmplen),str(tmpuni)]) + "\n")
	
out_file.close()

