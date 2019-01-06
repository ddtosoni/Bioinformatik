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

#nur für Testzwecke
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



datapath = "C:/Users/daphne/Documents/fhnw/Bioinformatik/Projekt_Plasmidannotation/datenPlasmidannotation/"

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
                if i -  plasmidlen < 0:
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
