#!/usr/bin/python

import os
import argparse
import re

parser = argparse.ArgumentParser(description='Accession file and eval path')
parser.add_argument('-d','--primerdirectory', type=str,help='directory where all results of primers are kept')
parser.add_argument('-o','--out', type=str,help='output')

args = parser.parse_args()
PDIR = args.primerdirectory
OUT = args.out

fwd="GTGYCAGCMGCCGCGGTAA"
rev="GGACTACNVGGGTWTCTAAT"
fwd_rc="TTACCGCGGCKGCTGRCAC"
rev_rc="ATTAGAWACCCBNGTAGTCC"
primers = [fwd, rev, fwd_rc, rev_rc]

print "Reading runfiles from ", PDIR
primf = os.listdir(PDIR)
primf1 = [f for f in primf if "_1_" in f]
primf1runs = list(set([re.sub('_1_.+txt', '', f) for f in primf1]))
#primf2 = [f for f in primf if re.match(r'_2_', f)]
#primf2runs = list(set([f.replace('_2_[A-Z]+.fastq') for f in primf2]))
primfs = [f for f in primf if "_1_" not in f]
primfsruns = list(set([re.sub('_.+txt', '',f) for f in primfs]))

#all = dict.fromkeys(list(set(primf1runs + primfsruns)))

out_file = open(OUT, "w")
out_file.write("\t".join(["run","both","fwd","rev","missing"]) + "\n")
for f in primf1runs:
	rundict = {"_1_": {p: [] for p in primers}, "_2_": {p: [] for p in primers}}
	for r in ["_1_", "_2_"]:
		for p in range(len(primers)):
			tmpfile = PDIR+"/"+f+r+primers[p]+".txt"
			tmpfile = open(tmpfile,"r")
			while 1:
				line = tmpfile.readline()
				#print line
				if line == "":
					break
				line = line.rstrip()
				tabs = line.split("\t")
				#print len(tabs)
				if len(tabs) <= 4:
					rundict[r][primers[p]].append(0)
				else:
					if len(tabs[4]) < len(tabs[6]):
						rundict[r][primers[p]].append(1)
					else:
						rundict[r][primers[p]].append(0)
			tmpfile.close()
	tmpbothfwd = [x + y == 2 for x,y in zip (rundict["_1_"][primers[0]], rundict["_2_"][primers[1]])]
	tmpbothrev = [x + y== 2 for x,y in zip (rundict["_1_"][primers[2]], rundict["_2_"][primers[3]])]
	tmpboth = sum([x or y for x,y in zip (tmpbothfwd, tmpbothrev)])
	tmp1 = [sum(x) == 1 for x in zip (*[rundict["_1_"][primers[0]], rundict["_1_"][primers[1]], rundict["_1_"][primers[2]], rundict["_1_"][primers[3]],
											rundict["_1_"][primers[0]], rundict["_1_"][primers[1]], rundict["_1_"][primers[2]], rundict["_1_"][primers[3]]])]
	tmpfwd1 = [x > 0 and y for x,y in zip (rundict["_1_"][primers[0]], tmp1)]
	tmpfwd2 = [x > 0 and y for x,y in zip (rundict["_1_"][primers[2]], tmp1)]
	tmpfwd3 = [x > 0 and y for x,y in zip (rundict["_2_"][primers[0]], tmp1)]
	tmpfwd4 = [x > 0 and y for x,y in zip (rundict["_2_"][primers[2]], tmp1)]
	tmprev1 = [x > 0 and y for x,y in zip (rundict["_1_"][primers[1]], tmp1)]
	tmprev2 = [x > 0 and y for x,y in zip (rundict["_1_"][primers[3]], tmp1)]
	tmprev3 = [x > 0 and y for x,y in zip (rundict["_2_"][primers[0]], tmp1)]
	tmprev4 = [x > 0 and y for x,y in zip (rundict["_2_"][primers[2]], tmp1)]
	tmpfwd = sum([w or x or y or z for w,x,y,z in zip (tmpfwd1, tmpfwd2,tmpfwd3,tmpfwd4)])
	tmprev = sum([w or x or y or z for w,x,y,z in zip (tmprev1, tmprev2,tmprev3,tmprev4)])
	tmp0 = sum([sum(x) == 0 for x in zip (*[rundict["_1_"][primers[0]], rundict["_1_"][primers[1]], rundict["_1_"][primers[2]], rundict["_1_"][primers[3]],
											rundict["_1_"][primers[0]], rundict["_1_"][primers[1]], rundict["_1_"][primers[2]], rundict["_1_"][primers[3]]])])
	out_file.write("\t".join([f,str(tmpboth),str(tmpfwd),str(tmprev),str(tmp0)]) + "\n")
	
for f in primfsruns:
	if f not in primf1runs:
		rundict = {p: [] for p in primers}
		for p in range(len(primers)):
			tmpfile = PDIR+"/"+f+"_"+primers[p]+".txt"
			tmpfile = open(tmpfile,"r")
			while 1:
				line = tmpfile.readline()
				if line == "":
					break
				line = line.rstrip()
				tabs = line.split("\t")
				if len(tabs) <= 4:
					rundict[primers[p]].append(0)
				else:
					if len(tabs[4]) < len(tabs[6]):
						rundict[primers[p]].append(1)
					else:
						rundict[primers[p]].append(0)
			tmpfile.close()
		tmpboth1= [x + y == 2 for x,y in zip (rundict[primers[0]], rundict[primers[1]])]
		tmpboth2= [x + y == 2 for x,y in zip (rundict[primers[2]], rundict[primers[3]])]
		tmpboth = sum([x or y for x,y in zip (tmpboth1, tmpboth2)])
		tmp1 = [sum(x) == 1 for x in zip (*[rundict[primers[0]], rundict[primers[1]], rundict[primers[2]], rundict[primers[3]]])]
		tmpfwd = [x + y > 0 for x,y in zip (rundict[primers[0]], rundict[primers[2]])]
		tmpfwd1 = sum([x and y for x,y in zip (tmp1, tmpfwd)])
		tmprev = [x + y > 0 for x,y in zip (rundict[primers[1]], rundict[primers[3]])]
		tmprev1 = sum([x and y for x,y in zip (tmp1, tmprev)])
		tmp0 = sum([sum(x) == 0 for x in zip (*[rundict[primers[0]], rundict[primers[1]], rundict[primers[2]], rundict[primers[3]]])])
		out_file.write("\t".join([f,str(tmpboth),str(tmpfwd1),str(tmprev1),str(tmp0)]) + "\n")


#for f in all.keys():
	
out_file.close()

