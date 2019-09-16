#!/usr/bin/env python

import os
import sys
sys.path.append('/home/users/aheintzbuschart/lib/')
import argparse
import re
import numpy
import collections
from pymongo import MongoClient

topdirectory = "/scratch/users/snarayanasamy/LAO_TS_IMP-v1.3/"
keggdirectory = "/scratch/users/snarayanasamy/LAO_TS_IMP-v1.3/Databases/Annotations/KEGG/"
drepdirectory ="/scratch/users/snarayanasamy/LAO_TS_IMP-v1.3/TemporalBinning/RepresentativeBins/data_tables/"
popdirectory = "/scratch/users/snarayanasamy/LAO_TS_IMP-v1.3/PopulationAnalysis/Calculations/GeneLevel/"
idfile = "/scratch/users/snarayanasamy/LAO_TS_IMP-v1.3/mongoDB/ids"

# get arguments from call
parser = argparse.ArgumentParser(description='Feed mongo DB with contig and gene annotations')
parser.add_argument('-l','--lib', type=str,help='sample ID EG A01')

args = parser.parse_args()
LIB = args.lib # combined assembly ID of this sample

# read file with higher order cluster information
drepdict = {}
drepFileAll = drepdirectory + "Cdb.csv"
header = 1
print "Reading top-level cluster information from ", drepFileAll
drep_file = open(drepFileAll, "r")
while 1:
    line = drep_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        if not line.startswith("Isolate"):
            tab = line.split(",") # 0:genome, 1: secondary_cluster, 2: threshold, 3: cluster_method, 4: comparison_algorithm, 5: primary_cluster
            drepdict[tab[0]] = {'primary_cluster': tab[5], 'secondary_cluster': tab[1]}
drep_file.close()

# read file with names of representative genomes
drepFileRepres = drepdirectory + "Widb.csv"
header = 1
print "Reading top-level cluster information from ", drepFileRepres
drep_file = open(drepFileRepres, "r")
while 1:
    line = drep_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split(",") # 0:N50, 1: closest_cluster_member, 2: cluster, 3: cluster_members, 4: completeness, 5: completeness_metric, 6: contamination,
        # 7: contamination_metric, 8: furthest_cluster_member, 9: genome, 10: score, 11: size, 12: strain_heterogeneity,
        # 13: tax_confidence, 14: taxonomy
        if tab[9] in drepdict:
            drepdict[tab[9]]['representsGenomes'] = tab[3]
drep_file.close()

# read file of contig information and make dictionary with all contigs
contig_dict = {}
contigFile = topdirectory + LIB + "/Analysis/mgmt.assembly.length.txt"
print "Reading contig information from ", contigFile
contig_file = open(contigFile, "r")
while 1:
    line = contig_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:sequenceID, 1:length
    contig_dict[tab[0]] = {'sample': LIB,'length' : int(tab[1]), 'GCperc' : 0, 'aveCov' : 0, 'cluster' : "S"}
contig_file.close()

# read file of contig GC content and update dictionary 
header = 1
contigFile = topdirectory + LIB + "/Analysis/mgmt.assembly.gc_content.txt"
print "Reading contig information from ", contigFile
contig_file = open(contigFile, "r")
while 1:
    line = contig_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:sequenceID, 1:GC
        contig_dict[tab[0]]['GCperc'] = float(tab[1])
contig_file.close()

# read file of DNA coverage of the same sample for contigs and update contig dictionary
covDNAcontigFile = topdirectory + LIB + "/Analysis/mg.assembly.contig_depth.txt"
print "Reading coverage data for DNA on contigs from ", covDNAcontigFile
covDNAcontig_file = open(covDNAcontigFile, "r")
while 1:
    line = covDNAcontig_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:SequenceID,	1:Average.coverage 
    contig_dict[tab[0]]['aveCov'] = float(tab[1])
covDNAcontig_file.close()

#read cluster membership
header = 1
clusterFile = topdirectory + LIB + "/Binning/contigs2clusters.10.4.tsv"
print "Reading contig cluster membership from ", clusterFile
cluster_file = open(clusterFile, "r")
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['cluster'] = tab[1]
        if tab[1] in drepdict:
            contig_dict[tab[0]]['primary_cluster'] = drepdict[tab[1]]['primary_cluster']
            contig_dict[tab[0]]['secondary_cluster'] = drepdict[tab[1]]['secondary_cluster']
            if "representsGenomes" in drepdict[tab[1]]:
                contig_dict[tab[0]]['representing'] = drepdict[tab[1]]['representsGenomes']
cluster_file.close()
contigNum = len(contig_dict)
print "gathered information on ", contigNum, "contigs"

# read file of prokka predictions and make dictionary with all predicted genes
protein_dict = {}
lct = 0
protFile = topdirectory + LIB + "/Analysis/annotation/annotation.filt.gff"
print "Reading protein information from ", protFile
prot_file = open(protFile, "r")
while 1:
    line = prot_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        lct += 1
        tab = line.split("\t") # 0:contig, 1:source, 2:kind, 3:start, 4:end, 5:.n, 6:sense, 7:0, 8:attributes
        if len(tab) >= 9:
            atts = tab[8].split(";")
            for att in atts:
                if att.startswith("ID="):
                    gene = att.replace("ID=","")
                    break
            if tab[3] != "" and tab[4]!= "":
                protein_dict[gene] = {'contig' : tab[0],'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]), 'kind' : tab[2], 'aveCovDNA' : 0, 'aveCovRNA' : 0}
                products = ""
                symbols =  ""
                ecs = ""
                inferences = ""
                for att in atts:
                    if att.startswith("product="):
                        products += att.replace("product=","") + ";"
                    if att.startswith("gene="):
                        symbols += att.replace("gene=","") + ";"
                    if att.startswith("eC_number="):
                        ecs += ecs + att.replace("eC_number=","") + ";"
                    if att.startswith("inference="):
                        infer = att.replace("inference=","")
                        infers = infer.split(":")
                        if "UniProtKB" in infers or "CLUSTERS" in infers or "Pfam" in infers or "HAMAP" in infers:
                            inferences += infers[4]
                #genes without these attributes are not annotated
                prods = products.split(";")
                del prods[-1]
                symbs = symbols.split(";")
                del symbs[-1]
                ecsl = ecs.split(";")
                del ecsl[-1]
                infs = inferences.split(";")
                if len(prods) > 0:
                    protein_dict[gene]['product'] = prods
                if len(symbs) > 0:
                    protein_dict[gene]['symbol'] = symbs
                if len(ecsl) > 0:
                    protein_dict[gene]['EC'] = ecsl
                if len(infs) > 0:
                    protein_dict[gene]['annoID'] = infs
            else:
                print "Missing genomic coordinates for ", gene, " in sample ", LIB
        else:
            print "missing attributes for the gene in line ", lct, "of ",  protFile
prot_file.close()

# read file of DNA coverage for genes and update gene dictionary
covDNAgeneFile = topdirectory + LIB + "/Analysis/mg.gene_depth.avg"
print "Reading coverage data for DNA on genes from ", covDNAgeneFile
covDNAgene_file = open(covDNAgeneFile, "r")
while 1:
    line = covDNAgene_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:geneID,	1:average depth
    protein_dict[tab[0]]['aveCovDNA'] = float(tab[1])
covDNAgene_file.close()

# read file of RNA coverage for genes and update gene dictionary
covRNAgeneFile = topdirectory + LIB + "/Analysis/mt.gene_depth.avg"
print "Reading forward coverage data for RNA on genes from ", covRNAgeneFile
covRNAgene_file = open(covRNAgeneFile, "r")
while 1:
    line = covRNAgene_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:geneID,	1:average depth
    protein_dict[tab[0]]['aveCovRNA'] = float(tab[1])
covRNAgene_file.close()

# read file with all IDs
print "Reading IDs of all samples from ", idfile
ids = []
id_file = open(idfile, "r")
while 1:
    line = id_file.readline()
    if line == "":
        break
    line = line.rstrip()
    ids.append(line)
covDNAgene_file.close()

# read files of DNA reads per gene by all samples for genes in representative clusters and update gene dictionary
for currlib in ids:
    readDNAgeneFile = popdirectory + currlib + ".mg.annotation.featureCounts.txt"
    print "Reading coverage data for DNA on genes from ", readDNAgeneFile
    readDNAgene_file = open(readDNAgeneFile, "r")
    while 1:
        line = readDNAgene_file.readline()
        if line == "":
            break
        if line.startswith(LIB):
            line = line.rstrip()
            tab = line.split("\t") # 0: Geneid, 1: Chr, 2: Start, 3: End, 4: Strand, 5: Length, 6: number of reads
            currgene = re.sub(".+_PROKKA","PROKKA",tab[0])
            if 'crossSampleDNAreads' not in protein_dict[currgene]:
                protein_dict[currgene]['crossSampleDNAreads'] = dict()
            protein_dict[currgene]['crossSampleDNAreads'][currlib] = int(tab[6])
    readDNAgene_file.close()

# read files of RNA reads per gene by all samples for genes in representative clusters and update gene dictionary
for currlib in ids:
    readRNAgeneFile = popdirectory + currlib + ".mt.annotation.featureCounts.txt"
    print "Reading coverage data for RNA on genes from ", readRNAgeneFile
    readRNAgene_file = open(readRNAgeneFile, "r")
    while 1:
        line = readRNAgene_file.readline()
        if line == "":
            break
        if line.startswith(LIB):
            line = line.rstrip()
            tab = line.split("\t") # 0: Geneid, 1: Chr, 2: Start, 3: End, 4: Strand, 5: Length, 6: number of reads
            currgene = re.sub(".+_PROKKA","PROKKA",tab[0])
            if 'crossSampleRNAreads' not in protein_dict[currgene]:
                protein_dict[currgene]['crossSampleRNAreads'] = dict()
            protein_dict[currgene]['crossSampleRNAreads'][currlib] = int(tab[6])
    readRNAgene_file.close()

# read file of essential genes and keep only best hits
header = 3
ess_dict = {}
essFile = topdirectory + LIB + "/Binning/ORFS.hmm.orfs.essential.hits"
print "Reading essential gene information from ", essFile
ess_file = open(essFile, "r")
while 1:
    line = ess_file.readline()
    if line == "":
        break
    if header > 0:
        header -= 1
    else:
        line = line.rstrip()
        tab = line.split() # 0:target name, 1:accession, 2:query name, 3: accession, 4: E-value, 5: score, 6: bias, 7: E-value, 8:score, 9: bias, 10: exp, 11: reg, 12: clu, 13: ov, 14: env, 15: dom, 16: rep, 17: inc, 18: description of target
        if "#" not in tab[0]:
            if tab[0] not in ess_dict:
                ess_dict[tab[0]] = [tab[2],float(tab[4])]
            elif float(tab[4]) < ess_dict[tab[0]][1]:
                ess_dict[tab[0]] = [tab[2],float(tab[4])]
ess_file.close()

# update gene dictionary with essentiality information
for ess in ess_dict:
    if ess in protein_dict:
            protein_dict[ess]['essentialGene'] = ess_dict[ess][0] # I do not make an entry for unessential genes
ess_dict = {}

# read kegg annotation file for genes and update gene dictionary
keggFile = keggdirectory + LIB + ".KOs_RXNs_CMPs.tsv"
print "Reading KEGG annotations of genes from ", keggFile
kegg_file = open(keggFile, "r")
kct = 0
while 1:
    line = kegg_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:Gene with sample ID etc, 1:KO with KEGG, 2:maxScore, 3:hitNumber, 4: KO again, 5:reaction IDs, 6: compound IDs, 7: compound names
    geneID = re.sub(".+_PROKKA","PROKKA",tab[0]) 
    if geneID in protein_dict:
        kofield = tab[1].replace("KEGG:","")
        kos = kofield.split(";")
        protein_dict[geneID]['KO'] = kos # I do not make an entry for unannotated genes
        compounds = tab[6].split(";")
        protein_dict[geneID]['compound'] = compounds 
kegg_file.close()


#insert protein information into the contig dictionary
print "inserting information on ", len(protein_dict), "genes into contig data"
contcount=0
for prot in protein_dict:
    protein_dict[prot]['gene'] = prot
    contig = protein_dict[prot]['contig']
    if contig not in contig_dict:
        print "no contig information available for protein", prot, " on ", contig
    else:
        if "genes" not in contig_dict[contig]:
            contig_dict[contig]['genes'] = [protein_dict[prot]]
            contcount += 1
        else:
            contig_dict[contig]['genes'].append(protein_dict[prot])
protein_dict = {}
            
print "gene information for ", contcount, " contigs"

#make clusterInfo-like table
clusterStatFile =  LIB + ".clusterInfoFromPy.tsv"
cluster_dict ={}
for cont in contig_dict:
    #XXX
    tmpCont = contig_dict[cont]
    complCnt = 0
    expCnt = 0
    expRead = 0
    proUCnt = 0
    proPCnt = 0
    ess = []
    ann = 0
    annKOL = 0
    annProkL = 0
    annECL = 0
    annSPL = 0
    annPfL = 0
    annTIL = 0
    trna = 0
    annKO = []
    annProk = []
    annEC = []
    annSP = []
    annPf = []
    annTI = []
    amPhy = []
    amGen = []
    mpPhy = []
    mpGen = []
    mp = []
    mbPhy = []
    mbGen = []
    mb = []
    geneno = 0
    if "genes" in tmpCont:
        tmpGenes = tmpCont['genes']
        geneno = len(tmpGenes)
        for g in tmpGenes:
            if g['aveCovRNA'] >= 0.01:
                expCnt += 1
            if g['kind'] == "tRNA":
                trna +=1
            if 'essentialGene' in g:
                ess.append(g['essentialGene'])
            if 'KO' in g or 'annoID' in g:
                ann += 1
            if 'KO' in g:
                annKO += g['KO']
                annKOL += 1
            if 'annoID' in g:
                annProk += g['annoID']
                annProkL += 1
            if 'EC' in g:
                annEC += g['EC']
                annECL += 1
    if tmpCont['cluster'] != "S":
        if tmpCont['cluster'] not in cluster_dict:
            #XXX
            cluster_dict[tmpCont['cluster']] = {'length' : tmpCont['length'],'contigs' : 1,'aveCov' : tmpCont['length']*tmpCont['aveCov'],
                                               'genes' : geneno, 'tRNA' : trna, 
                                               'expressedGenes' : expCnt, 'ess' : ess,'annotated' : ann,
                                                'annotatedKO' : annKOL, 'annotatedPROKKA' : annProkL, 'annotatedEC' : annECL, 
                                                'catsKO' : annKO, 'catsPROKKA' : annProk, 'catsEC' : annEC }
        else:
            cluster_dict[tmpCont['cluster']]['length'] += tmpCont['length']
            cluster_dict[tmpCont['cluster']]['contigs'] +=  1
            cluster_dict[tmpCont['cluster']]['aveCov'] += tmpCont['length']*tmpCont['aveCov']
            cluster_dict[tmpCont['cluster']]['genes'] += geneno
            cluster_dict[tmpCont['cluster']]['tRNA'] += trna
            cluster_dict[tmpCont['cluster']]['expressedGenes'] += expCnt
            cluster_dict[tmpCont['cluster']]['ess'] += ess
            cluster_dict[tmpCont['cluster']]['annotated'] += ann
            cluster_dict[tmpCont['cluster']]['annotatedKO'] += annKOL
            cluster_dict[tmpCont['cluster']]['annotatedPROKKA'] += annProkL
            cluster_dict[tmpCont['cluster']]['annotatedEC'] += annECL
            cluster_dict[tmpCont['cluster']]['catsKO'] += annKO
            cluster_dict[tmpCont['cluster']]['catsPROKKA'] += annProk
            cluster_dict[tmpCont['cluster']]['catsEC'] += annEC
    else:
        if "S" not in cluster_dict:
            cluster_dict['S'] = {'length' : tmpCont['length'],'contigs' : 1,'aveCov' : tmpCont['length']*tmpCont['aveCov'],
                                  'genes' : geneno, 'tRNA' : trna, 
                                  'expressedGenes' : expCnt,
                                  'ess' : ess,'annotated' : ann,
                                  'annotatedKO' : annKOL, 'annotatedPROKKA' : annProkL, 'annotatedEC' : annECL, 
                                  'catsKO' : annKO,'catsPROKKA' : annProk, 'catsEC' : annEC }
        else:
            cluster_dict['S']['length'] += tmpCont['length']
            cluster_dict['S']['contigs'] +=  1
            cluster_dict['S']['aveCov'] += tmpCont['length']*tmpCont['aveCov']
            cluster_dict['S']['genes'] += geneno
            cluster_dict['S']['tRNA'] += trna
            cluster_dict['S']['expressedGenes'] += expCnt
            cluster_dict['S']['ess'] += ess
            cluster_dict['S']['annotatedKO'] += annKOL
            cluster_dict['S']['annotatedPROKKA'] += annProkL
            cluster_dict['S']['annotatedEC'] += annECL
            cluster_dict['S']['catsKO'] += annKO
            cluster_dict['S']['catsPROKKA'] += annProk
print "collected data on ", len(cluster_dict), " clusters."
print "Writing cluster data to ", clusterStatFile
clusterStat_file = open(clusterStatFile, "w")
clusterStat_file.write("sample\tcluster\tlength\tcontigs\taveCov\tgenes\texpressedGenes\ttRNAloci\tuniqueEss\tnumEss\t"+
                       "catsKO\tcatsPROKKA\tcatsEC\tannotated\tannotatedKO\tannotatedPROKKA\tannotatedEC\n")
for clus in cluster_dict:
    length = str(cluster_dict[clus]['length'])
    contigs = str(cluster_dict[clus]['contigs'])
    aveCov = str(cluster_dict[clus]['aveCov']/cluster_dict[clus]['length'])
    genes = str(cluster_dict[clus]['genes'])
    exprGenes = str(cluster_dict[clus]['expressedGenes'])
    trnaU = str(cluster_dict[clus]['tRNA'])
    uniqueEss = str(len(set(cluster_dict[clus]['ess'])))
    numEss = str(len(cluster_dict[clus]['ess']))
    annotated = str(cluster_dict[clus]['annotated'])
    catsKO = str(len(set(cluster_dict[clus]['catsKO'])))
    annoKO = str(cluster_dict[clus]['annotatedKO'])
    catsPROKKA = str(len(set(cluster_dict[clus]['catsPROKKA'])))
    annoPROKKA = str(cluster_dict[clus]['annotatedPROKKA'])
    catsEC = str(len(set(cluster_dict[clus]['catsEC'])))
    annoEC = str(cluster_dict[clus]['annotatedEC'])
    writeClusterList = [LIB,clus,length,contigs,aveCov,genes,exprGenes,trnaU,uniqueEss,numEss,
                        catsKO,catsPROKKA,catsEC,annotated,annoKO,annoPROKKA,annoEC]
    writeCluster = "\t".join(writeClusterList)
    clusterStat_file.write(writeCluster+"\n")
    cluster_dict[clus] = []
clusterStat_file.close()
cluster_dict = {}

#insert contig dictionaries into MongoDB
client = MongoClient()
db = client['laodb']
oldCollSize = db.laots.count()
for cont in contig_dict:
    contig_dict[cont]['contig'] = cont
    db.laots.insert_one(contig_dict[cont])
newCollSize = db.laots.count()
print "contigs inserted into database:", newCollSize - oldCollSize
print "there are now", newCollSize, "documents in the collection."
