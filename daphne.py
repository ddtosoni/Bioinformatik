#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmidannotation

-- part of member 3 --

@author: daphne
"""

#TODO Stellen sind jetzt immer ab 0! Evtl. Änderungsbedarf
#TODO Längen durch 3 teilbar machen

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from random import *
import os
import pickle

#nur für Testzwecke
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

datapath = "/Users/deniztosoni/FHNW/bioinformatik/daten plasmid annotation 2.0/"

fiftyRandomPlasmids = []
i=0
#(TODO): aktuell keine Sicherheit, dass es wirklich 50 werden.. aber das ist gerade nicht so wichtig
for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
    x = 1
    #x = randint(0, 1)
    #todo 10 auf 50 setzen
    if (x == 1 and i < 10):
        fiftyRandomPlasmids.append(record)
        i = i+1
        print(record.id + "\t " + str(i))
        if i == 2:  break

#print("fiftyRandomPlasmids")
#print(fiftyRandomPlasmids)

"""
for record in fiftyRandomPlasmids:
    output_file = open("fileWithoutAnnotations.gb", "w")
    SeqIO.write(record, output_file, "genbank")
    feature = SeqFeature(FeatureLocation(start=3, end=12), type="primer_bind", qualifiers={"note":["Hello"]})
    record.features.append(feature)
    for feat in record.features:
        if feat.type == "primer_bind":
            print(feat.type)
            print(feat.qualifiers)
            print(feat.qualifiers['note'][0])
            feat.qualifiers["note"][0] = "ABC"
            feat.qualifiers["note"].append("DEF")
            feat.qualifiers["note"][1] = "GHI"
            print(feat.qualifiers['note'][0])
            print(feat.qualifiers)
        for i in feat.qualifiers:
            print(i)
    print(record.annotations["source"])
"""
o_file = open("fileWithAnnotations.gb", "w")
SeqIO.write(fiftyRandomPlasmids, o_file, "genbank")

#aus common_primer.mfasta Liste primerBindingSites erstellen
primerBindingSites = []
# da die Sequenz nur auf den Inhalt der letzten 15 Basen überprüft werden muss,
# werden die Enden separat in der Liste primerBindingEndings gespeichert
primerBindingSiteEndings = []
for record in SeqIO.parse(datapath + "common_primer.mfasta", "fasta"):
    primerBindingSites.append(record)
    recordEnding = record[-15:-1]+record[-1]
    primerBindingSiteEndings.append(recordEnding)

# aus tags_epitropes.mfasta Liste specialTranslatedFeatures erstellen
specialTranslatedFeatures = []
for record in SeqIO.parse(datapath + "tags_epitopes.mfasta", "fasta"):
    specialTranslatedFeatures.append(record)

"""für Testzwecke:
newfeature = SeqRecord(Seq("ILLVCF"), id="Testfeature")
print(newfeature.seq)
specialTranslatedFeatures.append(newfeature)
"""

#Teilaufgabe 1
#todo später Funktion daraus machen
for plasmid in fiftyRandomPlasmids:
    plasmidlen = len(plasmid)
    searchsequence = plasmid + plasmid
    for prim in primerBindingSiteEndings:

        seq = prim.seq
        seqRev = prim.reverse_complement().seq
        # Achtung, Länge des Motivs wird als Ausgangspunkt gewählt, d.h. bitte schon richtig übergeben
        seqlen = len(prim)

        #todo brauche für Verallgemeinerung noch motifRevComp, das überprüft wird
        i = 0
        while i < len(searchsequence):
            if searchsequence[i:i+seqlen].seq == seq:
                if i - plasmidlen < 0:
                    print("gefunden an Stelle: " + str(i) + " im Vorderstrang")
                    newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="primer_bind", strand=1, qualifiers={"note":[prim.name]})
                    plasmid.features.append(newPrimer)
                    #TODO  The feature.qualifier note shall also indicate in the whether there is only a partial (but perfect match to 15 bases of the 3’ end of the primer).
            if searchsequence[i:i+seqlen].seq == seqRev:
                if i - plasmidlen < 0:
                    print("gefunden an Stelle: " + str(plasmidlen - (i + seqlen)) + " im Gegenstrang")
                    newPrimer = SeqFeature(FeatureLocation(start=plasmidlen-(i+seqlen), end=plasmidlen-i), type="primer_bind", strand=-1, qualifiers={"note":[prim.name]})
                    plasmid.features.append(newPrimer)
                    #TODO  The feature.qualifier note shall also indicate in the whether there is only a partial (but perfect match to 15 bases of the 3’ end of the primer).
            i += 1

print("ready")

newFile = open("fileAfterPart12.gb", "w")
SeqIO.write(fiftyRandomPlasmids, newFile, "genbank")

for plasmid in fiftyRandomPlasmids:
    vorderstrang = plasmid + plasmid
    vorderstrang.name = "vorderstrang"
    gegenstrang = vorderstrang.reverse_complement()
    gegenstrang.name = "gegenstrang"
    straenge = [vorderstrang, gegenstrang]
    positions = [0,1,2]
    plasmidlen = len(plasmid)

    print("hello")
    openReadingFrames = []
    existingAnnotations = []
    for feat in plasmid.features:
        if feat.type == "CDS" or feat.type == "rep_origin":
            existingAnnotations.append(feat)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print("feat type start")
        print(feat.type)
        print(feat.location)
        print(feat.location.start)
        print("feat type end")

    for strang in straenge:
        print(straenge[1])
        for i in positions:
            print(i)
            pos = strang[i:]
            transPlas = pos.translate()
            print(transPlas.seq)
            transPlaslen = len(transPlas)

            #Teilaufgabe 2
            for specFeature in specialTranslatedFeatures:
                specFeaturelen = len(specFeature)
                j = 0

                while j < transPlaslen:
                    pos = transPlas[j:(j+specFeaturelen)]
                    if specFeature.seq == pos.seq:
                        print(len(plasmid))
                        featureLocation = i + j*3
                        if featureLocation < plasmidlen:
                            print("Found it! Strang " + strang.description + " +Position " + str(featureLocation))
                            if strang.name == "vorderstrang":
                                strand=1
                            if strang.name == "gegenstrang":
                                strand=-1
                            newSpec = SeqFeature(FeatureLocation(start=featureLocation, end=featureLocation+specFeaturelen*3), type="misc_feature", strand=strand, qualifiers={"note":[specFeature.name]})
                            plasmid.features.append(newSpec)
                        else:
                            print("Wiederholungstäter der Stelle " + str(featureLocation-plasmidlen))
                        print(specFeature.seq)
                        print(pos.seq)
                    j += 1

            newFile2 = open("fileAfterPart2.gb", "w")
            SeqIO.write(fiftyRandomPlasmids, newFile2, "genbank")

            #Teilaufgabe 3

            openedFrame = bool(0)
            print(openedFrame)
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
                    if openedFrame == bool(1):
                        end = k

                        # das * wird nicht mehr dazugezählt --> Auslesen bis zur Stelle des Endes
                        openReadingFrame = strang[begin*3+i:end*3+i]
                        search = bool(0)

                        #todo: später 199 auf 49 setzen
                        if (takeOnlyCareOfLongSequences > 49):
                            #TODO: vorher überprüfen, ob bereits im gb-file annotiert
                            for anno in existingAnnotations:
                                if (begin*3+i) != anno.location.start and (end*3+i) != anno.location.end:
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
                            hit_title = hit.alignments[0].title  #[:60]  # evtl. Titel kürzen
                            print(hit_title)
                            newHit = SeqFeature(FeatureLocation(start=begin*3+i, end=end*3+i), type="protein", strand=1, qualifiers={"note":[hit_title]})
                            plasmid.features.append(newHit)

                            os.remove("my_blast.xml")
                        else:
                            print("sequence is to short")
                        openedFrame = bool(0)
                    else:
                        print("counldn't find a new point to start for opening the reading frame")
                k += 1
                takeOnlyCareOfLongSequences += 1

newFile3 = open("fileAfterPart3.gb", "w")
SeqIO.write(fiftyRandomPlasmids, newFile3, "genbank")



## TEST THIS SHIT YOOO


from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from random import *
import os
import pickle

datapath = "/Users/deniztosoni/Downloads/Bioinformatik3/"

fiftyRandomPlasmids = []
i=0
#(TODO): aktuell keine Sicherheit, dass es wirklich 50 werden.. aber das ist gerade nicht so wichtig
for record in SeqIO.parse(datapath + "fileAfterPart3.gb", "genbank"):
    x = 1
    #x = randint(0, 1)
    #todo 10 auf 50 setzen
    if (x == 1 and i < 10):
        fiftyRandomPlasmids.append(record)
        i = i+1
        print(record.id + "\t " + str(i))
        if i == 1:  break

with open('testDump.pickle', 'rb') as handle:
    unserialized_data = pickle.load(handle)

print(unserialized_data.items())

for plasmid in fiftyRandomPlasmids:
    plasmidlen = len(plasmid)
    searchsequence = plasmid + plasmid
    for (drop, toAnnotate), additionalData in zip(unserialized_data.values(), unserialized_data.keys()):
#        print(toAnnotate)
#        print(additionalData)
        print("NEXT")

        print(additionalData.split('||')[0])
        print(additionalData.split('||')[1])

        seq = toAnnotate
        seqRev = toAnnotate.reverse_complement()
        # Achtung, Länge des Motivs wird als Ausgangspunkt gewählt, d.h. bitte schon richtig übergeben
        seqlen = len(toAnnotate)

        #todo brauche für Verallgemeinerung noch motifRevComp, das überprüft wird
        i = 0
        while i < len(searchsequence):
            if searchsequence[i:i+seqlen].seq == seq:
                if i - plasmidlen < 0:
                    print("gefunden an Stelle: " + str(i) + " im Vorderstrang")
                    if additionalData[:3] == 'LTR':
                        note = additionalData[7:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="LTR", strand=1, qualifiers={"note":note})
                    if additionalData[:4] == 'oriT':
                        note = additionalData[8:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="oriT", strand=1, qualifiers={"note":note})
                    if additionalData[:8] == 'enhancer':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="enhancer", strand=1, qualifiers={"note":note})
                    if additionalData[:8] == 'promoter':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="promoter", strand=1, qualifiers={"note":note})
                    if additionalData[:10] == 'rep_origin':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="rep_origin", strand=1, qualifiers={"note":note})
                    if additionalData[:11] == 'primer_bind':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="primer_bind", strand=1, qualifiers={"note":note})
                    if additionalData[:10] == 'terminator':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="terminator", strand=1, qualifiers={"note":note})
                    if additionalData[:11] == 'misc_signal':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="misc_signal", strand=1, qualifiers={"note":note})
                    if additionalData[:11] == 'misc_recomb':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="misc_recomb", strand=1, qualifiers={"note":note})
                    if additionalData[:10] == '-35_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="-35_signal", strand=1, qualifiers={"gene":gene,
                                                                                                                                 "product":product})
                    if additionalData[:10] == '-10_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="-10_signal", strand=1, qualifiers={"gene":gene,
                                                                                                                                 "product":product})
                    if additionalData[:3] == 'RBS':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="RBS", strand=1, qualifiers={"gene":gene,
                                                                                                                          "product":product})
                    if additionalData[:12] == 'polyA_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="polyA_signal", strand=1, qualifiers={"gene":gene,
                                                                                                                                   "product":product})
                    if additionalData[:11] == 'sig_peptide':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="sig_peptide", strand=1, qualifiers={"gene":gene,
                                                                                                                                  "product":product})
                    if additionalData[:3] == 'CDS': #gene and prod
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="CDS", strand=1, qualifiers={"gene":gene,
                                                                                                                          "product":product})
                    if additionalData[:12] == 'protein_bind':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="protein_bind", strand=1, qualifiers={"note":note,
                                                                                                                                   "bound_moiety":boundM})
                    if additionalData[:12] == 'misc_binding':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="misc_binding", strand=1, qualifiers={"note":note,
                                                                                                                                   "bound_moiety":boundM})
                    if additionalData[:14] == 'mobile_element':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        mobEmTyp = additionalData.split('||')[2]
                        mobEmTyp = mobEmTyp[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="mobile_element", strand=1, qualifiers={"note":note,
                                                                                                                                     "mobile element type":mobEmTyp})
                    if additionalData[:4] == 'mRNA':
                        gene = additionalData[8:]
                        gene = gene[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="mRNA", strand=1, qualifiers={"gene":gene})
                    if additionalData[:4] == 'tRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="tRNA", strand=1, qualifiers={"product":product})
                    if additionalData[:4] == 'rRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="rRNA", strand=1, qualifiers={"product":product})

                    plasmid.features.append(newPrimer)
                    #TODO  The feature.qualifier note shall also indicate in the whether there is only a partial (but perfect match to 15 bases of the 3’ end of the primer).
            if searchsequence[i:i+seqlen].seq == seqRev:
                if i - plasmidlen < 0:
                    print("gefunden an Stelle: " + str(plasmidlen - (i + seqlen)) + " im Gegenstrang")
                    if additionalData[:3] == 'LTR':
                        note = additionalData[7:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="LTR", strand=-1, qualifiers={"note":note})
                    if additionalData[:4] == 'oriT':
                        note = additionalData[8:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="oriT", strand=-1, qualifiers={"note":note})
                    if additionalData[:8] == 'enhancer':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="enhancer", strand=-1, qualifiers={"note":note})
                    if additionalData[:8] == 'promoter':
                        note = additionalData[12:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="promoter", strand=-1, qualifiers={"note":note})
                    if additionalData[:10] == 'rep_origin':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="rep_origin", strand=-1, qualifiers={"note":note})
                    if additionalData[:11] == 'primer_bind':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="primer_bind", strand=-1, qualifiers={"note":note})
                    if additionalData[:10] == 'terminator':
                        note = additionalData[14:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="terminator", strand=-1, qualifiers={"note":note})
                    if additionalData[:11] == 'misc_signal':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="misc_signal", strand=-1, qualifiers={"note":note})
                    if additionalData[:11] == 'misc_recomb':
                        note = additionalData[15:]
                        note = note[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="misc_recomb", strand=-1, qualifiers={"note":note})
                    if additionalData[:10] == '-35_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="-35_signal", strand=-1, qualifiers={"gene":gene,
                                                                                                                                "product":product})
                    if additionalData[:10] == '-10_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="-10_signal", strand=-1, qualifiers={"gene":gene,
                                                                                                                                "product":product})
                    if additionalData[:3] == 'RBS':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="RBS", strand=-1, qualifiers={"gene":gene,
                                                                                                                         "product":product})
                    if additionalData[:12] == 'polyA_signal':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="polyA_signal", strand=-1, qualifiers={"gene":gene,
                                                                                                                                  "product":product})
                    if additionalData[:11] == 'sig_peptide':
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="sig_peptide", strand=-1, qualifiers={"gene":gene,
                                                                                                                                 "product":product})
                    if additionalData[:3] == 'CDS': #gene and prod
                        gene = additionalData.split('||')[1]
                        gene = gene[2:-2]
                        product = additionalData.split('||')[2]
                        product = product[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="CDS", strand=-1, qualifiers={"gene":gene,
                                                                                                                         "product":product})
                    if additionalData[:12] == 'protein_bind':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="protein_bind", strand=-1, qualifiers={"note":note,
                                                                                                                                  "bound_moiety":boundM})
                    if additionalData[:12] == 'misc_binding':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        boundM = additionalData.split('||')[2]
                        boundM = boundM[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="misc_binding", strand=-1, qualifiers={"note":note,
                                                                                                                                  "bound_moiety":boundM})
                    if additionalData[:14] == 'mobile_element':
                        note = additionalData.split('||')[1]
                        note = note[2:-2]
                        mobEmTyp = additionalData.split('||')[2]
                        mobEmTyp = mobEmTyp[2:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="mobile_element", strand=-1, qualifiers={"note":note,
                                                                                                                                    "mobile element type":mobEmTyp})
                    if additionalData[:4] == 'mRNA':
                        gene = additionalData[8:]
                        gene = gene[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="mRNA", strand=-1, qualifiers={"gene":gene})
                    if additionalData[:4] == 'tRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="tRNA", strand=-1, qualifiers={"product":product})
                    if additionalData[:4] == 'rRNA':
                        product = additionalData[8:]
                        product = product[:-2]
                        newPrimer = SeqFeature(FeatureLocation(start=i, end=i+seqlen), type="rRNA", strand=-1, qualifiers={"product":product})

                    plasmid.features.append(newPrimer)
                    #TODO  The feature.qualifier note shall also indicate in the whether there is only a partial (but perfect match to 15 bases of the 3’ end of the primer).
            i += 1

print("ready")

newFile = open("CommonFeatures_TEST.gb", "w")
SeqIO.write(fiftyRandomPlasmids, newFile, "genbank")

