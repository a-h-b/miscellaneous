#!/usr/bin/python

import os
import sys
import argparse
import re
import numpy
import collections
from pymongo import MongoClient

# get arguments from call
parser = argparse.ArgumentParser(description='Feed mongo DB with contig and gene annotations')
parser.add_argument('-d','--topdirectory', type=str,help='directory where all samples are kept')
parser.add_argument('-l','--lib', type=str,help='sample ID ie W2I02')

args = parser.parse_args()
LIB = args.lib # combined assembly ID of this sample
DDIR = args.topdirectory
SDIR = DDIR + "/" + LIB 

#function to be used to simplify lists of taxa
def enumerateTaxa(taxList):
    kp = collections.Counter(taxList).most_common()
    init = 1
    taxString = ""
    for tax,occ in kp:
        if init == 1:
            init = 0
            taxString = tax + "(" + str(occ) + ");"
        else:
            taxString += tax + "(" + str(occ) + ");"
    taxString = taxString.rstrip(";")
    return taxString


# read file with selected MAG information
## here we want to have the MAG name, taxonomy, quality and GRiD result + sample of origin and method;
## all MAG will become one collection 
magdict = {}

dasFile = SDIR + "/Binning/" + "selected_DASTool_summary.txt"
print "Reading MAG information from ", dasFile
header = 1
das_file = open(dasFile, "r")
while 1:
    line = das_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        origin = "MetaBAT"
        line = line.rstrip()
        tab = line.split("\t") # 0: bin 1:uniqueBacSCGs 2:redundantBacSCGs 3:uniqueArcSCGs 4:redundantArcSCGs 5:bacRatio 6:arcRatio 7:size 8:contigs N50 9:binScore 10:SCG_completeness  11:SCG_redundancy
        if re.match(r"maxbin",tab[0]):
            origin = "MaxBin"
        elif re.match(r"\D",tab[0]):
            origin = "binny"
        magdict[tab[0]] = {'sample': LIB, 'MAG': tab[0], 'binScore': float(tab[9]), 'binningMethod': origin}
das_file.close()
print "Read MAG summary information"

for bin in magdict:
    gridFile = SDIR + "/Binning/" + "selected_DASTool_bins/" + bin + "/grid/reads.GRiD.txt"
    print "Reading MAG growth information from ", gridFile
    header = 1
    grid_file = open(gridFile, "r")
    while 1:
        line = grid_file.readline()
        if line == "":
            break
        if header == 1:
            header = 0
        else:
            origin = "MetaBAT"
            line = line.rstrip()
            tab = line.split("\t") # 0: Sample 1: GRiD 2:95% CI  3: GRiD unrefined  4: Species heterogeneity  5:Coverage  6:dnaA/ori ratio  7:ter/dif ratio
            magdict[bin].update([ ('GRiD', float(tab[1])) , ('depth', float(tab[5])) , ('dnaAToOri', float(tab[6])) , ('terToDif' , float(tab[7]))])
    grid_file.close()
print "Read MAG growth information"


gtdbBFile = SDIR + "/Binning/selected_GTDBtk_classifications/gtdbtk.bac120.summary.tsv"
gtdbAFile = SDIR + "/Binning/selected_GTDBtk_classifications/gtdbtk.ar122.summary.tsv"
if os.path.isfile(gtdbBFile):
    print "Reading MAG taxonomy from ", gtdbBFile
    header = 1
    gtdb_file = open(gtdbBFile, "r")
    while 1:
        line = gtdb_file.readline()
        if line == "":
            break
        if header == 1:
            header = 0
        else:
            line = line.rstrip()
            tab = line.split("\t") # 0: user_genome 1:classification 2:fastani_reference  3:fastani_taxonomy 4:fastani_ani 5:fastani_af 6:closest_placement_reference  7:closest_placement_taxonomy 8:closest_placement_ani 9:closest_placement_af  10:classification_method  11:note  12:other_related_references(genome_id,species_name,ANI,AF)
            magdict[tab[0]].update([('taxstring', tab[1]), ('taxMethod', tab[10])])
    gtdb_file.close()
    print "Read MAG bacterial taxonomy information"
    
if os.path.isfile(gtdbAFile):
    print "Reading MAG taxonomy from ", gtdbAFile
    header = 1
    gtdb_file = open(gtdbAFile, "r")
    while 1:
        line = gtdb_file.readline()
        if line == "":
            break
        if header == 1:
            header = 0
        else:
            line = line.rstrip()
            tab = line.split("\t") # 0: user_genome 1:classification 2:fastani_reference  3:fastani_taxonomy 4:fastani_ani 5:fastani_af 6:closest_placement_reference  7:closest_placement_taxonomy 8:closest_placement_ani 9:closest_placement_af  10:classification_method  11:note  12:other_related_references(genome_id,species_name,ANI,AF)
            magdict[tab[0]].update([('taxstring', tab[1]), ('taxMethod', tab[10])])
    gtdb_file.close()
    print "Read MAG archaean taxonomy information"


# read file of contig information and make dictionary with all contigs
contig_dict = {}
contigFile = SDIR + "/Analysis/mgmt.assembly.length.txt"
print "Reading contig information from ", contigFile
contig_file = open(contigFile, "r")
while 1:
    line = contig_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:sequenceID, 1:length
    contig_dict[tab[0]] = {'sample': LIB,'length' : int(tab[1]), 'GCperc' : 0, 'aveCov' : 0, 'binnyCluster' : "S"}
contig_file.close()

# read file of contig GC content and update dictionary 
header = 1
contigFile = SDIR + "/Analysis/mgmt.assembly.gc_content.txt"
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
covDNAcontigFile = SDIR + "/Analysis/mg.assembly.contig_depth.txt"
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
print "gathered basic information on ", len(contig_dict), "contigs"

# read coordinates
vizCnt = 0
coordFile = SDIR + "/Analysis/mgmt.vizbin.with-contig-names.points"
print "Reading contig BHSNE coordinates from ", coordFile
coord_file = open(coordFile, "r")
while 1:
    line = coord_file.readline()
    if line == "":
        break
    else:
        vizCnt += 1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig 1:x, 2:y
        contig_dict[tab[0]]['coords'] = [float(tab[1]),float(tab[2])] 
coord_file.close()
print "gathered vizbin coordinates of ", vizCnt, "contigs"

#read cluster membership
header = 1
binnyCnt = 0
clusterFile = SDIR + "/Binning/binny/contigs2clusters.10.4.tsv"
print "Reading contig cluster membership from ", clusterFile
cluster_file = open(clusterFile, "r")
while 1:
    line = cluster_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        binnyCnt +=  1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['binnyCluster'] = tab[1]
cluster_file.close()
print "gathered binny-based information on ", binnyCnt, "contigs"

#read maxbin
header = 1
maxbinCnt = 0
maxbinFile = SDIR + "/Binning/MaxBin/scaffold2bin.tsv"
print "Reading contig bin membership from ", maxbinFile
maxbin_file = open(maxbinFile, "r")
while 1:
    line = maxbin_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        maxbinCnt +=  1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['maxbin'] = tab[1]
maxbin_file.close()
print "gathered MaxBin-based information on ", maxbinCnt, "contigs"

#read metabat
header = 1
metabatCnt = 0
metabatFile = SDIR + "/Binning/MetaBAT/scaffold2bin.tsv"
print "Reading contig bin membership from ", metabatFile
metabat_file = open(metabatFile, "r")
while 1:
    line = metabat_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        metabatCnt +=  1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['metabatbin'] = tab[1]
metabat_file.close()
print "gathered MetaBAT-based information on ", metabatCnt, "contigs"

#read DAStool
header = 1
dastoolCnt = 0
dastoolFile = SDIR + "/Binning/selected_DASTool_scaffolds2bin.txt"
print "Reading contig bin membership from ", metabatFile
dastool_file = open(dastoolFile, "r")
while 1:
    line = dastool_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        dastoolCnt +=  1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:class (contig membership)
        contig_dict[tab[0]]['dastoolbin'] = tab[1]
dastool_file.close()
print "read DAStool information on ", dastoolCnt, "contigs"


# read Kraken data and update contig dictionary
header = 1
krakenFile = SDIR + "/Analysis/mgmt.GTDBkraken.parsed.tsv"
print "Reading contig Kraken information from ", krakenFile
kraken_file = open(krakenFile, "r")
while 1:
    line = kraken_file.readline()
    if line == "":
        break
    if header == 1:
        header = 0
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:length, 2:divSpec, 3:divGen, 4:divFam, 5:divOrd, 6:divClass, 7:divPhylum, 8:divKingdom,
        ### 9:annotationLevel, 10:species, 11:genus, 12:family, 13:order, 14:class, 15:phylum, 16:kingdom
        contig_dict[tab[0]]['krakenAnnotationLevel'] = tab[9]
        if tab[10] != "unknown":
            contig_dict[tab[0]]['krakenSpecies'] = tab[10]
        if tab[10] != "unknown":
            contig_dict[tab[0]]['krakenGenus'] = tab[11]
        if tab[11] != "unknown":
            contig_dict[tab[0]]['krakenFamily'] = tab[12]
        if tab[12] != "unknown":
            contig_dict[tab[0]]['krakenOrder'] = tab[13]
        if tab[13] != "unknown":
            contig_dict[tab[0]]['krakenClass'] = tab[14]
        if tab[14] != "unknown":
            contig_dict[tab[0]]['krakenPhylum'] = tab[15]
        if tab[15] != "unknown":
            contig_dict[tab[0]]['krakenKingdom'] = tab[16]
        # I only use annotated contigs; the names of the fields are different from R, because I cannot have "." in the keys
kraken_file.close()




# read file of prokka and barrnap etc predictions and make dictionary with all predicted genes
protein_dict = {}
rrna_dict = {}
trna_dict = {}
pct = 0
rct = 0
tct = 0
lct = 0
protFile = SDIR + "/Analysis/annotation/annotation.CDS.RNA.hmms.gff"
print "Reading protein information from ", protFile
prot_file = open(protFile, "r")
while 1:
    line = prot_file.readline()
    if line == "":
        break
    else:
        lct += 1
        line = line.rstrip()
        tab = line.split("\t") # 0:contig, 1:source, 2:kind, 3:start, 4:end, 5:.n, 6:sense, 7:0, 8:attributes
        if len(tab) >= 9:
            atts = tab[8].split(";")
            for att in atts:
                if att.startswith("ID="):
                    gene = att.replace("ID=","")
                    break
            if tab[3] != "" and tab[4]!= "":
                if tab[2] == "CDS":
                    pct += 1
                    protein_dict[gene] = {'contig' : tab[0],'prokkaID' : "",'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]), 'completeness': "complete",'kind' : tab[2], 'aveCovDNA' : 0, 'aveCovRNA' : 0, 'readsDNA' : 0, 'readsRNA': 0}
                    if int(tab[3]) == 1 or int(tab[4]) == contig_dict[tab[0]]["length"]:
                        protein_dict[gene]["completeness"] = "incomplete"
                    KEGG = ""
                    Pfam_A =  ""
                    essential = ""
                    dbCAN = ""
                    Resfam = ""
                    for att in atts:
                        if att.startswith("KEGG="):
                            protein_dict[gene]['KEGG'] = att.replace("KEGG=","")
                        if att.startswith("Pfam_A="):
                            protein_dict[gene]['Pfam'] = att.replace("Pfam_A=","")
                        if att.startswith("essential="):
                            protein_dict[gene]['essential'] = att.replace("essential=","")
                        if att.startswith("dbCAN="):
                            protein_dict[gene]['dbCAN'] = att.replace("dbCAN=","")
                        if att.startswith("Resfams="):
                            protein_dict[gene]['Resfam'] = att.replace("Resfams=","")
                #genes without these attributes are not annotated
                elif tab[2] == "rRNA":
                    rct += 1
                    rrna_dict[gene] = {'contig' : tab[0],'progigalID' : "",'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]), 'kind' : tab[2], 'aveCovDNA' : 0, 'aveCovRNA' : 0}
                    for att in atts:
                        if att.startswith("product="):
                            att.replace("product=","")
                            if att.endswith(" (partial)"):
                                rrna_dict[gene]['completeness'] = "incomplete"
                                rrna_dict[gene]['rRNA'] = att.replace(" ribosomal RNA (partial)","")
                            else:
                                rrna_dict[gene]['completeness'] = "complete"
                                rrna_dict[gene]['rRNA'] = att.replace(" ribosomal RNA","")
                elif tab[2] == "tRNA":    
                    tct += 1
                    trna_dict[gene] = {'contig' : tab[0],'prodigalID' : "",'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]), 'kind' : tab[2], 'aveCovDNA' : 0, 'aveCovRNA' : 0}
                    for att in atts:
                        if att.startswith("product="):
                            trna_dict[gene]['tRNA'] = att.replace("product=","")
                elif tab[2] == "tmRNA":    
                    tct += 1
                    trna_dict[gene] = {'contig' : tab[0],'prodigalID' : "",'sense' : tab[6], 'start' : int(tab[3]), 'end' : int(tab[4]), 'length' : 1+int(tab[4])-int(tab[3]), 'kind' : tab[2], 'aveCovDNA' : 0, 'aveCovRNA' : 0}
                    for att in atts:
                        if att.startswith("product="):
                            trna_dict[gene]['tRNA'] = att.replace("product=","")
            else:
                print "Missing genomic coordinates for ", gene, " in sample ", LIB
        else:
            print "missing attributes for the gene in line ", lct, "of ",  protFile
prot_file.close()

print "read information on ", lct, " regions: ", pct, " proteins, ", rct, " rRNAs, ", tct, " tRNAs."

#read prodigal-style name for genes:
prodIDFile = SDIR + "/Analysis/annotation/annotation.filt.contig2ID.tsv"
print "Reading alternative names from ", prodIDFile
prodID_file = open(prodIDFile, "r")
while 1:
    line = prodID_file.readline()
    if line == "":
        break
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:contig_name, 1:prokka_name
        if rrna_dict.has_key(tab[1]):
            rrna_dict[tab[1]]['prodigalID'] = tab[0]
        elif trna_dict.has_key(tab[1]):
            trna_dict[tab[1]]['prodigalID'] = tab[0]
        elif protein_dict.has_key(tab[1]):
            protein_dict[tab[1]]['prodigalID'] = tab[0]
prodID_file.close()
print "Read alternative names."

# read file of DNA coverage for genes and update gene dictionary
covDNAgeneFile = SDIR + "/Analysis/mg.gene_depth.avg"
print "Reading coverage data for DNA on genes from ", covDNAgeneFile
covDNAgene_file = open(covDNAgeneFile, "r")
while 1:
    line = covDNAgene_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:geneID,	1:average depth
    if protein_dict.has_key(tab[0]):
        protein_dict[tab[0]]['aveCovDNA'] = float(tab[1])
    elif rrna_dict.has_key(tab[0]):
        rrna_dict[tab[0]]['aveCovDNA'] = float(tab[1])
    elif trna_dict.has_key(tab[0]):
        trna_dict[tab[0]]['aveCovDNA'] = float(tab[1])
covDNAgene_file.close()

# read file of RNA coverage for genes and update gene dictionary
covRNAgeneFile = SDIR + "/Analysis/mt.gene_depth.avg"
print "Reading forward coverage data for RNA on genes from ", covRNAgeneFile
covRNAgene_file = open(covRNAgeneFile, "r")
while 1:
    line = covRNAgene_file.readline()
    if line == "":
        break
    line = line.rstrip()
    tab = line.split("\t") # 0:geneID,	1:average depth
    if protein_dict.has_key(tab[0]):
        protein_dict[tab[0]]['aveCovRNA'] = float(tab[1])
    elif rrna_dict.has_key(tab[0]):
        rrna_dict[tab[0]]['aveCovRNA'] = float(tab[1])
    elif trna_dict.has_key(tab[0]):
        trna_dict[tab[0]]['aveCovRNA'] = float(tab[1])
covRNAgene_file.close()

# read files of RNA reads per gene and update protein dictionary
readsDNAgeneFile = SDIR + "/Analysis/mg.CDS_counts.tsv"
print "Reading forward coverage data for DNA on genes from ", covDNAgeneFile
header = 2
readsDNAgene_file = open(readsDNAgeneFile, "r")
while 1:
    line = readsDNAgene_file.readline()
    if line == "":
        break
    if header > 0:
        header -= 1
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Geneid  1:Chr 2:Start  3:End  4:Strand  5:Length 6:Assembly/mg.reads.sorted.ba
        if protein_dict.has_key(tab[0]):
            protein_dict[tab[0]]['readsDNA'] = float(tab[6])
readsDNAgene_file.close()

# read files of RNA reads per gene and update protein dictionary
readsRNAgeneFile = SDIR + "/Analysis/mt.CDS_counts.tsv"
print "Reading forward coverage data for RNA on genes from ", covRNAgeneFile
header = 2
readsRNAgene_file = open(readsRNAgeneFile, "r")
while 1:
    line = readsRNAgene_file.readline()
    if line == "":
        break
    if header > 0:
        header -= 1
    else:
        line = line.rstrip()
        tab = line.split("\t") # 0:Geneid  1:Chr 2:Start  3:End  4:Strand  5:Length 6:Assembly/mt.reads.sorted.bam
        if protein_dict.has_key(tab[0]):
            protein_dict[tab[0]]['readsRNA'] = float(tab[6])
readsRNAgene_file.close()


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

print "inserting information on ", len(rrna_dict), "rRNA genes into contig data"
for rna in rrna_dict:
    rrna_dict[rna]['gene'] = rna
    contig = rrna_dict[rna]['contig']
    if contig not in contig_dict:
        print "no contig information available for rRNA", rna, " on ", contig
    else:
        if "genes" not in contig_dict[contig]:
            contig_dict[contig]['genes'] = [rrna_dict[rna]]
            contcount += 1
        else:
            contig_dict[contig]['genes'].append(rrna_dict[rna])
rrna_dict = {}

print "inserting information on ", len(trna_dict), "rRNA genes into contig data"
for rna in trna_dict:
    trna_dict[rna]['gene'] = rna
    contig = trna_dict[rna]['contig']
    if contig not in contig_dict:
        print "no contig information available for tRNA", rna, " on ", contig
    else:
        if "genes" not in contig_dict[contig]:
            contig_dict[contig]['genes'] = [trna_dict[rna]]
            contcount += 1
        else:
            contig_dict[contig]['genes'].append(trna_dict[rna])
trna_dict = {}
print "gene information for ", contcount, " contigs"


#make clusterInfo-like table
# we have MAG name, taxonomy, quality and GRiD result + sample of origin and method;
# we want to add number of genes, number of rRNAs, tRNAs, number of essential genes, duplication of essential genes,
#      number of annotated and expressed genes, summary of kraken-result
clusterStatFile =  SDIR + "/clusterInfoFromPy.tsv"
bin_dict ={}
for cont in contig_dict:
    #XXX
    tmpCont = contig_dict[cont]
    complCnt = 0
    expCnt = 0
    expRead = 0
    ess = []
    ann = 0
    annKOL = 0
    annPfamL = 0
    annRESL = 0
    annCANL = 0
    trna = 0
    rrna = 0
    cds = 0
    amPhy = []
    amGen = []
    mpPhy = []
    mpGen = []
    mp = []
    mbPhy = []
    mbGen = []
    mb = []
    kP = []
    kC = []
    kO = []
    kF = []
    kG = []
    x = []
    y = []
    geneno = 0
    if tmpCont.has_key('dastoolbin'):
        if "genes" in tmpCont:
            tmpGenes = tmpCont['genes']
            geneno = len(tmpGenes)
            for g in tmpGenes:
                if g['aveCovRNA'] >= 0.01:
                    expCnt += 1
                if g['kind'] == "CDS":
                    cds +=1
                    if g['completeness'] == "complete":
                        complCnt += 1
                if g['kind'] == "tRNA":
                    trna +=1
                if g['kind'] == "rRNA":
                    rrna +=1
                if 'essential' in g:
                    ess.append(g['essential'])
                if 'KEGG' in g or 'Pfam' in g or 'dbCAN' in g or 'Resfam' or 'essential' in g:
                    ann += 1
                if 'KO' in g:
                    annKOL += 1
                if 'Pfam' in g:
                    annPfamL += 1
                if 'dbCAN' in g:
                    annCANL += 1
                if 'Resfam' in g:
                    annRESL += 1
        if "krakenPhylum" in tmpCont:
            kP = [tmpCont['krakenPhylum']]
        if "krakenClass" in tmpCont:
            kC = [tmpCont['krakenClass']]
        if "krakenOrder" in tmpCont:
            kO = [tmpCont['krakenOrder']]
        if "krakenFamily" in tmpCont:
            kF = [tmpCont['krakenFamily']]
        if "krakenGenus" in tmpCont:
            kG = [tmpCont['krakenGenus']]
        if "coords" in tmpCont:
            x = [tmpCont['coords'][0]]
            y = [tmpCont['coords'][1]]
        if tmpCont['dastoolbin'] not in bin_dict:
            bin_dict[tmpCont['dastoolbin']] = {'length' : tmpCont['length'],'contigs' : 1,'aveCov' : tmpCont['length']*tmpCont['aveCov'],
                                               'CDS' : cds, 'completeCDS' : complCnt, 'tRNA' : trna, 'rRNA' : rrna, 
                                               'expressedGenes' : expCnt, 'ess' : ess,'annotated' : ann,
                                                'annotatedKO' : annKOL, 'annotatedPfam' : annPfamL, 'annotatedDbCAN' : annCANL,
                                                'annotatedResfam' : annRESL, 'krakenPhylum': kP, 'krakenClass': kC, 'krakenOrder': kO,
                                                'krakenFamily': kF, 'krakenGenus': kG,
                                                'coordX' : x, 'coordY' : y}
        else:
            bin_dict[tmpCont['dastoolbin']]['length'] += tmpCont['length']
            bin_dict[tmpCont['dastoolbin']]['contigs'] +=  1
            bin_dict[tmpCont['dastoolbin']]['aveCov'] += tmpCont['length']*tmpCont['aveCov']
            bin_dict[tmpCont['dastoolbin']]['CDS'] += cds
            bin_dict[tmpCont['dastoolbin']]['completeCDS'] += complCnt
            bin_dict[tmpCont['dastoolbin']]['rRNA'] += rrna
            bin_dict[tmpCont['dastoolbin']]['tRNA'] += trna
            bin_dict[tmpCont['dastoolbin']]['expressedGenes'] += expCnt
            bin_dict[tmpCont['dastoolbin']]['ess'] += ess
            bin_dict[tmpCont['dastoolbin']]['annotated'] += ann
            bin_dict[tmpCont['dastoolbin']]['annotatedKO'] += annKOL
            bin_dict[tmpCont['dastoolbin']]['annotatedPfam'] += annPfamL
            bin_dict[tmpCont['dastoolbin']]['annotatedDbCAN'] += annCANL
            bin_dict[tmpCont['dastoolbin']]['annotatedResfam'] += annRESL
            bin_dict[tmpCont['dastoolbin']]['krakenPhylum'] += kP
            bin_dict[tmpCont['dastoolbin']]['krakenClass'] += kC
            bin_dict[tmpCont['dastoolbin']]['krakenOrder'] += kO
            bin_dict[tmpCont['dastoolbin']]['krakenFamily'] += kF
            bin_dict[tmpCont['dastoolbin']]['krakenGenus'] += kG
            bin_dict[tmpCont['dastoolbin']]['coordX'].append(x)
            bin_dict[tmpCont['dastoolbin']]['coordY'].append(y)

print "collected data on ", len(bin_dict), " bins."
print "Writing cluster data to ", clusterStatFile
clusterStat_file = open(clusterStatFile, "w")
clusterStat_file.write("sample\tbin\tlength\tcontigs\taveCov\tCDS\tcompleteCDS\texpressedGenes\trRNAloci\ttRNAloci\tuniqueEss\tnumEss\t"+
                       "annotatedGenes\tannotatedKEGG\tannotatedPfam\tannotatedDbCAN\tannotatedResfam\tphylaKraken\tclassesKraken\t"+
                       "ordersKraken\tfamiliesKraken\tgeneraKraken\tx\ty\n")
for bin in bin_dict:
    length = str(bin_dict[bin]['length'])
    contigs = str(bin_dict[bin]['contigs'])
    aveCov = str(bin_dict[bin]['aveCov']/bin_dict[bin]['length'])
    cdsU = str(bin_dict[bin]['CDS'])
    complC = str(bin_dict[bin]['completeCDS'])
    exprGenes = str(bin_dict[bin]['expressedGenes'])
    rrnaU = str(bin_dict[bin]['rRNA'])
    trnaU = str(bin_dict[bin]['tRNA'])
    uniqueEss = str(len(set(bin_dict[bin]['ess'])))
    numEss = str(len(bin_dict[bin]['ess']))
    annotated = str(bin_dict[bin]['annotated'])
    annoKO = str(bin_dict[bin]['annotatedKO'])
    annoPfam = str(bin_dict[bin]['annotatedPfam'])
    annoCAN = str(bin_dict[bin]['annotatedDbCAN'])
    annoResfam = str(bin_dict[bin]['annotatedResfam'])
    phylaKraken = enumerateTaxa(bin_dict[bin]['krakenPhylum'])
    classKraken = enumerateTaxa(bin_dict[bin]['krakenClass'])
    orderKraken = enumerateTaxa(bin_dict[bin]['krakenOrder'])
    familyKraken = enumerateTaxa(bin_dict[bin]['krakenFamily'])
    generaKraken = enumerateTaxa(bin_dict[bin]['krakenGenus'])
    coordX = str(numpy.median(numpy.array(bin_dict[bin]['coordX'])))
    coordY = str(numpy.median(numpy.array(bin_dict[bin]['coordY'])))
    magdict[bin].update([('length',length),('contigs',contigs),('aveCov',aveCov),('CDS',cdsU),('completeCDS',complC),
                         ('expressedGenes',exprGenes),('rRNA',rrnaU),('tRNA',trnaU),('uniqueEss',uniqueEss),
                         ('ess',numEss),('annotated',annotated),('annotatedKO',annoKO),
                         ('annotatedPfam',annoPfam),('annotatedDbCAN',annoCAN),('annotatedResfam',annoResfam),
                         ('krakenPhylum',phylaKraken),('krakenClass',classKraken),('krakenOrder',orderKraken),
                         ('krakenFamily',familyKraken),('krakenGenus',generaKraken),('coordX',coordX),
                         ('coordY',coordY)])
    writeClusterList = [LIB,bin,length,contigs,aveCov,cdsU,complC,exprGenes,rrnaU,trnaU,uniqueEss,numEss,
                        annotated,annoKO,annoPfam,annoCAN,annoResfam,phylaKraken,classKraken,orderKraken,
                        familyKraken,generaKraken,coordX,coordY]
    writeCluster = "\t".join(writeClusterList)
    clusterStat_file.write(writeCluster+"\n")
    bin_dict[bin] = []
clusterStat_file.close()
bin_dict = {}

#insert contig dictionaries into MongoDB
client = MongoClient()
db = client['litterDB']
print "Filling contig information into the database."
oldCollSize = db.W2Icontigs.count()
for cont in contig_dict:
    contig_dict[cont]['contig'] = cont
    db.W2Icontigs.insert_one(contig_dict[cont])
newCollSize = db.W2Icontigs.count()
print "contigs inserted into database:", newCollSize - oldCollSize
print "there are now", newCollSize, "documents in the contig collection."
print "Filling bin information into the database."
oldCollSize = db.W2Ibins.count()
for bin in magdict:
    magdict[bin]['bin'] = bin
    db.W2Ibins.insert_one(magdict[bin])
newCollSize = db.W2Ibins.count()
print "contigs inserted into database:", newCollSize - oldCollSize
print "there are now", newCollSize, "documents in the collection."
