from collections import Counter
from itertools import chain
from collections import defaultdict
import pickle
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
import os

"""
Part One Member 1 and 2
@author Deniz Tosoni, Christopher Galli
"""


'''File to read'''
inputFile = "vectors.gb"

'''Initialise dictionary'''
featureDict = dict()
testDict = dict()

'''Initialise Containers'''
counterStorage = []

# Define interested tags!
interested_only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb',
                        'LTR', 'enhancer']  # note

interested_note_and_gen = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']

interested_gene_and_product = ['CDS']

interested_note_and_bound_moiety = ['protein_bind', 'misc_binding']

interested_note_and_mobile_element_type = ['mobile_element']

interested_gene = ['mRNA']

interested_product = ['tRNA', 'rRNA']

other_features = ['5\'UTR', 'RBS', 'exon', 'intron', '3\'UTR']

'''Read in the File'''
record = SeqIO.parse(inputFile, "genbank")

'''Get feature.types and feature.qualifiers with the defined tags'''
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


'''Manage the Dictionary'''
CounterInformation = Counter(counterStorage)  # Creates a Dictionary with all the counts

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

for (x, y), z in zip(finalDictionary.values(), finalDictionary.keys()):
        insertKey = ''.join(str(vals) + ', ' for vals in x)
        insertKey = insertKey[:-2]
        if testIt.get(insertKey) is None:
            testIt[insertKey] = y, z
        newVal = testIt.get(insertKey)[0]
        if y > newVal:
            testIt[insertKey] = y, z

print(len(testIt))

with open('testDump.pickle', 'wb') as handle:
    pickle.dump(testIt, handle, protocol=pickle.HIGHEST_PROTOCOL)


"""
Plasmidannotation

-- part of member 3 --

@author: daphne
"""

dataPath = ""

fiftyRandomPlasmids = []
counter = 0

for record in SeqIO.parse(dataPath + "vectors.gb", "genbank"):
    if counter < 10:
        fiftyRandomPlasmids.append(record)
        counter += 1
        print(record.id + "\t " + str(counter))
        if counter == 2:
            break


o_file = open("fileWithAnnotations.gb", "w")

SeqIO.write(fiftyRandomPlasmids, o_file, "genbank")

primerBindingSites = []
primerBindingSiteEndings = []

for record in SeqIO.parse(dataPath + "common_primer.mfasta", "fasta"):
    primerBindingSites.append(record)
    recordEnding = record[-15:-1]+record[-1]
    primerBindingSiteEndings.append(recordEnding)

specialTranslatedFeatures = []

for record in SeqIO.parse(dataPath + "tags_epitopes.mfasta", "fasta"):
    specialTranslatedFeatures.append(record)


'''Sub-Solution 1'''
for plasmid in fiftyRandomPlasmids:

    plasmidLen = len(plasmid)
    searchSequence = plasmid + plasmid

    for prim in primerBindingSiteEndings:

        seq = prim.seq
        seqRev = prim.reverse_complement().seq
        seqLen = len(prim)

        counter = 0

        while counter < len(searchSequence):
            if searchSequence[counter:counter + seqLen].seq == seq:
                if counter - plasmidLen < 0:

                    print("Found at position: " + str(counter) + " in main string")

                    newPrimer = SeqFeature(FeatureLocation(start=counter,
                                           end=counter + seqLen),
                                           type="primer_bind",
                                           strand=1,
                                           qualifiers={"note": [prim.name]})

                    plasmid.features.append(newPrimer)

            if searchSequence[counter:counter + seqLen].seq == seqRev:
                if counter - plasmidLen < 0:

                    print("Found at position: " + str(plasmidLen - (counter + seqLen)) + " complement string")

                    newPrimer = SeqFeature(FeatureLocation(start=plasmidLen - (counter + seqLen),
                                           end=plasmidLen - counter),
                                           type="primer_bind",
                                           strand=-1,
                                           qualifiers={"note": [prim.name]})

                    plasmid.features.append(newPrimer)

            counter += 1

print("ready")

newFile = open("fileAfterPart12.gb", "w")

SeqIO.write(fiftyRandomPlasmids, newFile, "genbank")

for plasmid in fiftyRandomPlasmids:
    mainString = plasmid + plasmid
    mainString.name = "main string"
    complementString = mainString.reverse_complement()
    complementString.name = "complement string"
    strings = [mainString, complementString]
    positions = [0, 1, 2]
    plasmidLen = len(plasmid)
    openReadingFrames = []
    existingAnnotations = []

    for feat in plasmid.features:
        if feat.type == "CDS" or feat.type == "rep_origin":
            existingAnnotations.append(feat)

        print("feat type start")
        print(feat.type)
        print(feat.location)
        print(feat.location.start)
        print("feat type end")

    for string in strings:

        print(strings[1])

        for counter in positions:
            print(counter)
            pos = string[counter:]
            transPlasmid = pos.translate()
            print(transPlasmid.seq)
            transPlasmidLen = len(transPlasmid)

            '''Sub-Solution 2'''
            for specFeature in specialTranslatedFeatures:
                specFeatureLen = len(specFeature)

                j = 0

                while j < transPlasmidLen:
                    pos = transPlasmid[j:(j + specFeatureLen)]
                    if specFeature.seq == pos.seq:
                        print(len(plasmid))
                        featureLocation = counter + j * 3
                        if featureLocation < plasmidLen:
                            print("Found it! String " + string.description + " +position " + str(featureLocation))
                            if string.name == "main string":
                                strand = 1
                            if string.name == "complement string":
                                strand = -1
                            newSpec = SeqFeature(FeatureLocation(start=featureLocation,
                                                 end=featureLocation + specFeatureLen * 3),
                                                 type="misc_feature",
                                                 strand=strand,
                                                 qualifiers={"note": [specFeature.name]})

                            plasmid.features.append(newSpec)

                        else:
                            print("Repeat at position: " + str(featureLocation - plasmidLen))

                        print(specFeature.seq)
                        print(pos.seq)
                    j += 1

            newFile2 = open("fileAfterPart2.gb", "w")
            SeqIO.write(fiftyRandomPlasmids, newFile2, "genbank")

            '''Sub-Solution 3'''

            openedFrame = bool(0)
            begin = 0
            end = 0
            takeOnlyCareOfLongSequences = 0

            k = 0

            while k < transPlasmidLen:
                if transPlasmid[k] == "M":
                    begin = k
                    openedFrame = bool(1)
                    takeOnlyCareOfLongSequences = 0
                if transPlasmid[k] == "*":
                    if openedFrame == bool(1):
                        end = k
                        openReadingFrame = string[begin * 3 + counter:end * 3 + counter]
                        search = bool(0)

                        if takeOnlyCareOfLongSequences > 49:
                            for anno in existingAnnotations:
                                if (begin * 3 + counter) != anno.location.start and (end * 3 + counter)\
                                        != anno.location.end:
                                    search = bool(1)

                        if search == bool(1):
                            openReadingFrames.append(openReadingFrame)

                            print("Langes, offenes Leseraster gefunden! Warte auf Rückmeldung der online-BLAST-Suche.")
                            print("Suche Sequenz: " + openReadingFrame)
                            result_handle = NCBIWWW.qblast("blastx", "refseq_protein", openReadingFrame.seq)

                            blast_result = open("my_blast.xml", "w")
                            blast_result.write(result_handle.read())
                            blast_result.close()
                            result_handle.close()

                            hit = NCBIXML.read(open("my_blast.xml"))
                            hit_title = hit.alignments[0].title
                            print(hit_title)
                            newHit = SeqFeature(FeatureLocation(start=begin * 3 + counter,
                                                                end=end * 3 + counter),
                                                type="protein",
                                                strand=1,
                                                qualifiers={"note": [hit_title]})

                            plasmid.features.append(newHit)

                            os.remove("my_blast.xml")
                        else:
                            print("Sequence is to short")
                        openedFrame = bool(0)
                    else:
                        print("counldn't find a new point to start for opening the reading frame")
                k += 1
                takeOnlyCareOfLongSequences += 1




newFile3 = open("fileAfterPart3.gb", "w")
SeqIO.write(fiftyRandomPlasmids, newFile3, "genbank")


with open('testDump.pickle', 'rb') as handle:
    unserializedData = pickle.load(handle)

print(unserializedData.items())

for plasmid in fiftyRandomPlasmids:
    plasmidLen = len(plasmid)
    searchSequence = plasmid + plasmid
    for (drop, toAnnotate), additionalData in zip(unserializedData.values(), unserializedData.keys()):

        print("NEXT")

        print(additionalData.split('||')[0])
        print(additionalData.split('||')[1])

        seq = toAnnotate
        seqRev = toAnnotate.reverse_complement()
        seqLen = len(toAnnotate)

        counter = 0
        while counter < len(searchSequence):
            if searchSequence[counter:counter + seqLen].seq == seq:
                if counter - plasmidLen < 0:
                    print("Found at position: " + str(counter) + " in main string")
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
                    print("Fount at position: " + str(plasmidLen - (counter + seqLen)) + " in complement string")
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

print("ready")

newFile = open("CommonFeatures_TEST.gb", "w")
SeqIO.write(fiftyRandomPlasmids, newFile, "genbank")
