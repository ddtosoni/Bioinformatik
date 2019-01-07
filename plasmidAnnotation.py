"""
Plasmidannotation

Bioinformatik - Projekt 2

@author:
Auer Daphne
Galli Christopher
Tosoni Deniz
"""

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from random import *
from collections import Counter
from itertools import chain
from collections import defaultdict
import pickle

# Create the dictionary from vectors
'''File to read'''
inputFile = "vectors.gb"

'''Initialise dictionary'''
featureDict = dict()
testDict = dict()

'''Initialise Containers'''
storage = []

# Define interested tags
interested_only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb',
                        'LTR', 'enhancer']
interested_note_and_gen = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']
interested_gene_and_product = ['CDS']
interested_note_and_bound_moiety = ['protein_bind', 'misc_binding']
interested_note_and_mobile_element_type = ['mobile_element']
interested_gene = ['mRNA']
interested_product = ['tRNA', 'rRNA']
other_features = ['5\'UTR', 'RBS', 'exon', 'intron', '3\'UTR']

# test_dict = SeqIO.to_dict(SeqIO.parse("vectors-100.gb", "genbank"))

record = SeqIO.parse(inputFile, "genbank")

for search in record:
    if len(search.seq) >= 1500:
        for feature in search.features:
            if feature.type in interested_only_note:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'])
                    storage.append(feature.location.extract(search).seq)
                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_note_and_gen:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['gene'])
                    storage.append(feature.location.extract(search).seq)
                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_gene_and_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'],
                                                                         feature.qualifiers['gene'])
                    storage.append(feature.location.extract(search).seq)
                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_note_and_bound_moiety:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['bound_moiety'])
                    storage.append(feature.location.extract(search).seq)
                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_note_and_mobile_element_type:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['mobile_element_type'])
                    storage.append(feature.location.extract(search).seq)
                except:
                    break

            if feature.type in interested_gene:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['gene'])
                    storage.append(feature.location.extract(search).seq)
                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'])
                    storage.append(feature.location.extract(search).seq)
                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in  other_features:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['gene'],
                                                                         feature.qualifiers['note'])
                    storage.append(feature.location.extract(search).seq)
                except:
                    break

'''Manage the Dictionary'''
CounterInformation = Counter(storage)

for k in list(CounterInformation):
    if CounterInformation[k] <= 2:
        del featureDict[k]
        del CounterInformation[k]

testDict.update(CounterInformation)

'''Create correct Dictionary'''
finalDictionary = defaultdict(list)
for a, b in chain(featureDict.items(), testDict.items()):
    finalDictionary[a].append(b)

for k in list(finalDictionary.keys()):
    if len(k) <= 10:
        del finalDictionary[k]
        del featureDict[k]

testIt = dict()

# keys:     TYPE||QUALIFIER(||QUALIFIER)   <- if second qualifier is asked for
# values:   Count, SeqObject
for (x, y), z in zip(finalDictionary.values(), finalDictionary.keys()):
    insertKey = ''.join(str(vals) + '||' for vals in x)
    insertKey = insertKey[:-2]
    if testIt.get(insertKey) is None:
        testIt[insertKey] = y, z
    newVal = testIt.get(insertKey)[0]
    if y > newVal:
        testIt[insertKey] = y, z

for i in testIt.items():
    print(i)

# save data as pickle file
with open('testDump.pickle', 'wb') as handle:
    pickle.dump(testIt, handle, protocol=pickle.HIGHEST_PROTOCOL)

# read the pickle file
with open('testDump.pickle', 'rb') as handle:
    unserialized_data = pickle.load(handle)

# is it really the same? -> yes!
# print(testIt == unserialized_data)

datapath = ""

fiftyRandomPlasmids = []
i=0
for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
    x = randint(0, 1)
    if x == 1 and i < 50:
        fiftyRandomPlasmids.append(record)
        i = i+1
        print("untersuche folgende Records in der angegebenen Reihenfolge:")
        print(record.id + "\t " + str(i))

if i < 50:
    for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
        fiftyRandomPlasmids.append(record)
        i += 1
        if i == 49:
            break

# common_primer.mfasta -> list primerBindingSites
primerBindingSites = []
for record in SeqIO.parse(datapath + "common_primer.mfasta", "fasta"):
    recordEnding = record[-15:-1]+record[-1]
    primerBindingSites.append([record, recordEnding])

# tags_epitropes.mfasta -> list specialTranslatedFeatures
specialTranslatedFeatures = []
for record in SeqIO.parse(datapath + "tags_epitopes.mfasta", "fasta"):
    specialTranslatedFeatures.append(record)

# task 1
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
                    print("Primer auf dem Vorderstrang gefunden. Position: " + str(i))
                    print(longseqStart)
                    print(longseqEnd)

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
                    print("Primer auf dem Gegenstrang gefunden. Position: " + str(longseqRevStart))
                    print(shortseqRevStart)
                    print(longseqRevEnd)

            i += 1

    print("ready")

    newFile = open("fileAfterPart1.gb", "w")
    SeqIO.write(fiftyRandomPlasmids, newFile, "genbank")

    vorderstrang = plasmid + plasmid
    vorderstrang.name = "vorderstrang"
    gegenstrang = vorderstrang.reverse_complement()
    gegenstrang.name = "gegenstrang"
    straenge = [vorderstrang, gegenstrang]
    positions = [0, 1, 2]

    openReadingFrames = []
    existingAnnotations = []
    for feat in plasmid.features:
        if feat.type == "CDS" or feat.type == "rep_origin":
            existingAnnotations.append(feat)

    for strang in straenge:
        print("behandle den " + strang.name)

        for i in positions:
            print("position: " + str(i))

            pos = strang[i:plasmidlen+i]
            transPlas = pos.translate()
            print("frame sieht folgendermaßen aus: ")
            print(transPlas.seq)
            transPlaslen = len(transPlas)

            # task 2
            for specFeature in specialTranslatedFeatures:
                specFeaturelen = len(specFeature)
                j = 0

                while j < transPlaslen:
                    pos = transPlas[j:(j+specFeaturelen)]
                    if specFeature.seq == pos.seq:
                        featureLocation = i + j*3
                        if featureLocation < plasmidlen:
                            print("special feature auf dem " + strang.name + " an Position " + str(featureLocation) + " gefunden.")
                            if strang.name == "vorderstrang":   strand=1
                            if strang.name == "gegenstrang":    strand=-1
                            newSpec = SeqFeature(FeatureLocation(start=featureLocation, end=featureLocation+specFeaturelen*3), type="misc_feature", strand=strand, qualifiers={"note":[specFeature.name]})
                            plasmid.features.append(newSpec)
                    j += 1

            newFile2 = open("fileAfterPart2.gb", "w")
            SeqIO.write(fiftyRandomPlasmids, newFile2, "genbank")

            # task 3

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

                    # if an "M" was found, search until "*"
                    if openedFrame == bool(1):
                        end = k

                        # "*" doesn't count
                        openReadingFrame = strang[begin*3+i:end*3+i]
                        search = bool(0)

                        if (takeOnlyCareOfLongSequences > 49):
                            # is already annotated?
                            for anno in existingAnnotations:
                                if (begin*3+i) != anno.location.start and (end*3+i) != anno.location.end:
                                    search = bool(1)

                        if search == bool(1):
                            if strang.name == "vorderstrang":   strand=1
                            if strang.name == "gegenstrang":    strand=-1

                            openReadingFrames.append([openReadingFrame, begin*3+1, end*3+i, strand])
                            print("langes, offenes Leseraster gefunden an Stelle " + str(begin*3+1))

                            openedFrame = bool(0)
                k += 1
                takeOnlyCareOfLongSequences += 1

            print("Habe alle offenen Leseraster des Plasmids aufgenommen.")

    # -> BLAST it!
    print("BLAST-Suche startet.")
    for opFr in openReadingFrames:

        # frame is SeqRecord
        frame = opFr[0]

        # search in protein database
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
        print(sameletters + " von " + len(frame.seq) + " stimmen überein.")

        newHit = SeqFeature(FeatureLocation(start=opFr[1], end=opFr[2]), type="protein", strand=opFr[3], qualifiers={"note":[hit_title, "similarity = "+str(similarity)]})
        plasmid.features.append(newHit)

    newFile3 = open("fileAfterPart3.gb", "w")
    SeqIO.write(fiftyRandomPlasmids, newFile3, "genbank")

o_file = open("fileWithAnnotations.gb", "w")
SeqIO.write(fiftyRandomPlasmids, o_file, "genbank")