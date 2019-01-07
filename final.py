#---------------------------------
# Plasmidannotation
#
#---------------------------------
# Bioinformatik -  Projekt 2
#
#---------------------------------
# Date:           07.01.2019
#
#---------------------------------
# Authors:
#
# Auer Daphne
# Galli Christopher
# Tosoni Deniz
#
#---------------------------------
# Description:
#
# This program combines the following
# approaches:
#
# 1.
# Based on a number of annotated plasmid
# sequences in .gb format, common
# features BLABLA WILL DO MORE, THIS RUNS!!
#
#---------------------------------
# Imports:
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from collections import Counter
from collections import defaultdict
from itertools import chain
from random import *
import pickle
#---------------------------------
# File to read
inputFile = "/Users/deniztosoni/FHNW/bioinformatik/daten plasmid annotation 2.0/vectors.gb"

# Initialise dictionary
featureDict = dict()
testDict = dict()

# Initialise Containers
counterStorage = []

# Define tags
interested_only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb',
                        'LTR', 'enhancer']
interested_note_and_gen = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
interested_gene_and_product = ['CDS']
interested_note_and_bound_moiety = ['protein_bind', 'misc_binding']
interested_note_and_mobile_element_type = ['mobile_element']
interested_gene = ['mRNA']
interested_product = ['tRNA', 'rRNA']
other_features = ['5\'UTR', 'RBS', 'exon', 'intron', '3\'UTR']

# Read File
record = SeqIO.parse(inputFile, "genbank")

# Get feature.types and feature.qualifiers with the defined tags
for search in record:
    if len(search.seq) >= 1500:
        for feature in search.features:
            if feature.type in interested_only_note:
                try:  # Try because some of the entry don't match the requirements!
                    test = []
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_note_and_gen:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['gene'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_gene_and_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'],
                                                                         feature.qualifiers['gene'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_note_and_bound_moiety:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['bound_moiety'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_note_and_mobile_element_type:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['mobile_element_type'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_gene:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['gene'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in other_features:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['gene'],
                                                                         feature.qualifiers['note'])
                    counterStorage.append(feature.location.extract(search).seq)
                except:
                    break

# Manage the Dictionary
CounterInformation = Counter(counterStorage)  # Creates a Dictionary with all the counts

for k in list(CounterInformation):
    if CounterInformation[k] <= 2:
        del featureDict[k]
        del CounterInformation[k]

testDict.update(CounterInformation)

# finalize the dictionary
finalDictionary = defaultdict(list)
for a, b in chain(featureDict.items(), testDict.items()):
    finalDictionary[a].append(b)

for k in list(finalDictionary.keys()):
    if len(k) <= 10:
        del finalDictionary[k]
        del featureDict[k]

resultDict = dict()

for (x, y), z in zip(finalDictionary.values(), finalDictionary.keys()):
    insertKey = ''.join(str(vals) + ', ' for vals in x)
    insertKey = insertKey[:-2]
    if resultDict.get(insertKey) is None:
        resultDict[insertKey] = y, z
    newVal = resultDict.get(insertKey)[0]
    if y > newVal:
        resultDict[insertKey] = y, z

with open('testDump.pickle', 'wb') as handle:
    pickle.dump(resultDict, handle, protocol=pickle.HIGHEST_PROTOCOL)

datapath = "/Users/deniztosoni/FHNW/bioinformatik/daten plasmid annotation 2.0/"

# Take 50 random plasmids from the file
fiftyRandomPlasmids = []
i = 0
for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
    x = randint(0, 1)
    if x == 1 and i < 50:
        fiftyRandomPlasmids.append(record)
        i = i + 1
        print("Following records will be checked:")
        print(record.id + "\t " + str(i))

if i < 50:
    for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
        fiftyRandomPlasmids.append(record)
        i += 1
        if i == 49:
            break

# Create primerBindingSites list from common_primer.mfasta
primerBindingSites = []
for record in SeqIO.parse(datapath + "common_primer.mfasta", "fasta"):
    recordEnding = record[-15:-1]+record[-1]
    primerBindingSites.append([record,recordEnding])

# Create specialTranslatedFeatures from tags_epitropes.mfasta
specialTranslatedFeatures = []
for record in SeqIO.parse(datapath + "tags_epitopes.mfasta", "fasta"):
    specialTranslatedFeatures.append(record)

# Task 1
for plasmid in fiftyRandomPlasmids:

    plasmidlen = len(plasmid)
    searchsequence = plasmid + plasmid
    searchsequencelen = len(searchsequence)

    longseqStart = 0
    longseqEnd = 0
    longseqStartRev = 0
    longseqEndRev = 0

    for prim in primerBindingSites:

        shortseq = prim[1].seq
        shortseqRev = prim[1].reverse_complement().seq
        shortseqlen = len(shortseq)

        longseq = prim[0].seq
        longseqRev = prim[0].reverse_complement().seq
        longseqlen = len(longseq)

        i = 0
        while i < searchsequencelen:
            if i - plasmidlen < 0:
                # Search on strand 1
                longseqStart = i
                longseqEnd = i+longseqlen

                shortseqStart = i + (longseqlen-shortseqlen)
                foundPrim = bool(0)

                if searchsequence[longseqStart:longseqEnd].seq == longseq:
                    foundPrim = bool(1)
                    newPrimer = SeqFeature(FeatureLocation(start=longseqStart, end=longseqEnd), type="primer_bind", strand=1, qualifiers={"note":[prim[0].name, "perfect match"]})
                elif searchsequence[shortseqStart:longseqEnd].seq == shortseq:
                    foundPrim = bool(1)
                    newPrimer = SeqFeature(FeatureLocation(start=longseqStart, end=longseqEnd), type="primer_bind", strand=1, qualifiers={"note":[prim[0].name, "partial match"]})

                if foundPrim == bool(1):
                    plasmid.features.append(newPrimer)
                    print("Primer found on strand 1, at position: " + str(i))
                    print(longseqStart)
                    print(longseqEnd)

                # Reverse complement to find positions on strand -1
                longseqRevStart = plasmidlen - (i + longseqlen)
                longseqRevEnd = plasmidlen - i

                shortseqRevStart = plasmidlen - (i + shortseqlen)
                foundRevPrim = bool(0)

                if searchsequence[longseqRevStart:longseqRevEnd].seq == longseqRev:
                    foundRevPrim = bool(1)
                    newPrimer = SeqFeature(FeatureLocation(start=longseqRevStart, end=longseqRevEnd), type="primer_bind", strand=-1, qualifiers={"note":[prim[0].name, "perfect match"]})
                elif searchsequence[shortseqRevStart:longseqEnd].seq == shortseqRev:
                    foundRevPrim = bool(1)
                    newPrimer = SeqFeature(FeatureLocation(start=longseqRevStart, end=longseqRevEnd), type="primer_bind", strand=-1, qualifiers={"note":[prim[0].name, "partial match"]})

                if foundRevPrim == bool(1):
                    plasmid.features.append(newPrimer)
                    print("Primer found on strand -1, at position:" + str(longseqRevStart))
                    print(shortseqRevStart)
                    print(longseqRevEnd)

            i += 1

    print("ready")

# Testoutput file after first task
    # newFile = open("fileAfterPart1.gb", "w")
    # SeqIO.write(fiftyRandomPlasmids, newFile, "genbank")

    vorderstrang = plasmid + plasmid
    vorderstrang.name = "strand 1"
    gegenstrang = vorderstrang.reverse_complement()
    gegenstrang.name = "strand -1"
    straenge = [vorderstrang, gegenstrang]
    positions = [0, 1, 2]

    openReadingFrames = []
    existingAnnotations = []
    for feat in plasmid.features:
        if feat.type == "CDS" or feat.type == "rep_origin":
            existingAnnotations.append(feat)

    for strang in straenge:
        print("working on " + strang.name)

        for i in positions:
            print("position: " + str(i))

            pos = strang[i:plasmidlen+i]
            transPlas = pos.translate()
            print("frame: ")
            print(transPlas.seq)
            transPlaslen = len(transPlas)

            # Task 2
            for specFeature in specialTranslatedFeatures:
                specFeaturelen = len(specFeature)
                j = 0

                while j < transPlaslen:
                    pos = transPlas[j:(j+specFeaturelen)]
                    if specFeature.seq == pos.seq:
                        featureLocation = i + j*3
                        if featureLocation < plasmidlen:
                            print("special feature on " + strang.name + " at position " + str(featureLocation))
                            if strang.name == "strand 1":
                                strand=1
                            if strang.name == "strand -1":
                                strand=-1
                            newSpec = SeqFeature(FeatureLocation(start=featureLocation,
                                                                 end=featureLocation+specFeaturelen*3),
                                                 type="misc_feature", strand=strand,
                                                 qualifiers={"note":[specFeature.name]})
                            plasmid.features.append(newSpec)
                        # print(specFeature.seq)
                        # print(pos.seq)
                    j += 1

# Testoutput file after second task
            # newFile2 = open("fileAfterPart2.gb", "w")
            # SeqIO.write(fiftyRandomPlasmids, newFile2, "genbank")

            # Task 3
            openedFrame = bool(0)
            begin = 0
            end = 0
            takeOnlyCareOfLongSequences = 0

            k = 0
            while k < transPlaslen:

                if transPlas[k] == "M":
                    begin = k
                    openedFrame = bool(1)
                    takeOnlyCareOfLongSequences = 0

                if transPlas[k] == "*":

                    # only if "M" was found -> open reading frame -> closed by "*"
                    if openedFrame == bool(1):
                        end = k
                        # "*" doesnt count
                        openReadingFrame = strang[begin*3+i:end*3+i]
                        search = bool(0)

                        if (takeOnlyCareOfLongSequences > 49):
                            # check if already annotated
                            for anno in existingAnnotations:
                                if (begin*3+i) != anno.location.start and (end*3+i) != anno.location.end:
                                    search = bool(1)

                        # Take open reading frame if long enough & not yet annotated
                        if search == bool(1):
                            if strang.name == "strand 1":
                                strand = 1
                            if strang.name == "strand -1":
                                strand = -1
                            openReadingFrames.append([openReadingFrame, begin*3+1, end*3+i, strand])
                            print("long open reading frame at: " + str(begin*3+1))
                            openedFrame = bool(0)
                k += 1
                takeOnlyCareOfLongSequences += 1

            print("Finished open reading frames.")

    # BLAST
    print("BLAST search starting.")
    for opFr in openReadingFrames:

        # Frame is SeqRecord of sequence
        # opFr[1], opFr[2] contain position of start/end of sequence
        frame = opFr[0]

        result_handle = NCBIWWW.qblast("blastx", "refseq_protein", frame.seq)

        blastRecord = NCBIXML.read(result_handle)
        print("blastRecord")
        print(blastRecord)

        hit_title = blastRecord.alignments[0].title
        print(hit_title)
        hit_seq = blastRecord.alignments[0].hsps[0].query
        print(hit_seq)

        i = 0
        sameletters = 0
        for letter in frame.translate().seq:
            if i == len(hit_seq):
                break
            if hit_seq[i] == letter:
                try:
                    sameletters += 1
                except:
                    print("out of bound")

            i += 1

        similarity = sameletters / len(frame.seq)

        newHit = SeqFeature(FeatureLocation(start=opFr[1], end=opFr[2]), type="protein", strand=opFr[3], qualifiers={"note":[hit_title, "similarity = "+str(similarity)]})
        plasmid.features.append(newHit)

# Testoutput file after third task
    # newFile3 = open("fileAfterPart3.gb", "w")
    # SeqIO.write(fiftyRandomPlasmids, newFile3, "genbank")

with open('testDump.pickle', 'rb') as handle:
    unserializedData = pickle.load(handle)

# print(unserializedData.items())

for plasmid in fiftyRandomPlasmids:
    plasmidLen = len(plasmid)
    searchSequence = plasmid + plasmid
    for (drop, toAnnotate), additionalData in zip(unserializedData.values(), unserializedData.keys()):

        seq = toAnnotate
        seqRev = toAnnotate.reverse_complement()
        seqLen = len(toAnnotate)

        counter = 0
        while counter < len(searchSequence):
            if searchSequence[counter:counter + seqLen].seq == seq:
                if counter - plasmidLen < 0:
                    print("Found at position: " + str(counter) + " in main strand")
                    if additionalData[:3] == 'LTR':
                        note = additionalData[7:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="LTR",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:4] == 'oriT':
                        note = additionalData[8:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="oriT",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:8] == 'enhancer':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="enhancer",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:8] == 'promoter':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="promoter",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:10] == 'rep_origin':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="rep_origin",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:11] == 'primer_bind':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="primer_bind",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:10] == 'terminator':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="terminator",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:11] == 'misc_signal':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="misc_signal",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:11] == 'misc_recomb':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="misc_recomb",
                                               strand=1, qualifiers={"note": note})
                    if additionalData[:10] == '-35_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="-35_signal",
                                               strand=1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:10] == '-10_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="-10_signal",
                                               strand=1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:3] == 'RBS':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="RBS",
                                               strand=1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:12] == 'polyA_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="polyA_signal",
                                               strand=1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:11] == 'sig_peptide':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="sig_peptide",
                                               strand=1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:3] == 'CDS':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="CDS",
                                               strand=1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:12] == 'protein_bind':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="protein_bind",
                                               strand=1, qualifiers={"note": note, "bound_moiety": boundM})
                    if additionalData[:12] == 'misc_binding':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="misc_binding",
                                               strand=1, qualifiers={"note": note, "bound_moiety": boundM})
                    if additionalData[:14] == 'mobile_element':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        mobEmTyp = additionalData.split('||')[2]
                        mobEmTyp = mobEmTyp[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="mobile_element", strand=1,
                                               qualifiers={"note": note, "mobile element type": mobEmTyp})
                    if additionalData[:4] == 'mRNA':
                        gene = additionalData[8:]
                        gene = gene[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="mRNA",
                                               strand=1, qualifiers={"gene": gene})
                    if additionalData[:4] == 'tRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="tRNA",
                                               strand=1, qualifiers={"product": product})
                    if additionalData[:4] == 'rRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="rRNA",
                                               strand=1, qualifiers={"product": product})

                    plasmid.features.append(newPrimer)

            if searchSequence[counter:counter + seqLen].seq == seqRev:
                if counter - plasmidLen < 0:
                    print("Fount at position: " + str(plasmidLen - (counter + seqLen)) + " in complement strand")
                    if additionalData[:3] == 'LTR':
                        note = additionalData[7:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="LTR",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:4] == 'oriT':
                        note = additionalData[8:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="oriT",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:8] == 'enhancer':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="enhancer",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:8] == 'promoter':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="promoter",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:10] == 'rep_origin':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="rep_origin",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:11] == 'primer_bind':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="primer_bind",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:10] == 'terminator':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="terminator",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:11] == 'misc_signal':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="misc_signal",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:11] == 'misc_recomb':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="misc_recomb",
                                               strand=-1, qualifiers={"note": note})
                    if additionalData[:10] == '-35_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="-35_signal",
                                               strand=-1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:10] == '-10_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="-10_signal",
                                               strand=-1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:3] == 'RBS':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="RBS",
                                               strand=-1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:12] == 'polyA_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="polyA_signal",
                                               strand=-1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:11] == 'sig_peptide':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="sig_peptide",
                                               strand=-1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:3] == 'CDS':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="CDS",
                                               strand=-1, qualifiers={"gene": gene, "product": product})
                    if additionalData[:12] == 'protein_bind':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="protein_bind",
                                               strand=-1, qualifiers={"note": note, "bound_moiety": boundM})
                    if additionalData[:12] == 'misc_binding':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="misc_binding",
                                               strand=-1, qualifiers={"note": note, "bound_moiety": boundM})
                    if additionalData[:14] == 'mobile_element':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        mobEmTyp = additionalData.split('||')[2]
                        mobEmTyp = mobEmTyp[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen),
                                               type="mobile_element",
                                               strand=-1, qualifiers={"note": note, "mobile element type": mobEmTyp})
                    if additionalData[:4] == 'mRNA':
                        gene = additionalData[8:]
                        gene = gene[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="mRNA",
                                               strand=-1, qualifiers={"gene": gene})
                    if additionalData[:4] == 'tRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="tRNA",
                                               strand=-1, qualifiers={"product": product})
                    if additionalData[:4] == 'rRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=counter, end=counter + seqLen), type="rRNA",
                                               strand=-1, qualifiers={"product": product})

                    plasmid.features.append(newPrimer)
            counter += 1

o_file = open("finalOutput.gb", "w")
SeqIO.write(fiftyRandomPlasmids, o_file, "genbank")

print("Done.")
