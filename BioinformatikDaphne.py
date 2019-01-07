#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmidannotation

-- part of member 3 --

@author: daphne
"""


#TODO Stellen sind jetzt immer ab 0! Evtl. Änderungsbedarf
#TODO Längen durch 3 teilbar machen


from Bio import SeqIO, SearchIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from random import *
import os
from Bio.SeqUtils import six_frame_translations

#nur für Testzwecke
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



datapath = "C:/Users/daphne/Documents/fhnw/Bioinformatik/Projekt_Plasmidannotation/datenPlasmidannotation/"

fiftyRandomPlasmids = []
i=0
for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
    x = randint(0, 1)
    if (x == 1 and i < 50):
        fiftyRandomPlasmids.append(record)
        i = i+1
        print("untersuche folgende Records in der angegebenen Reihenfolge:")
        print(record.id + "\t " + str(i))


if i < 50:
    for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
        fiftyRandomPlasmids.append(record)
        i += 1
        if i == 49:    break



#aus common_primer.mfasta Liste primerBindingSites erstellen
primerBindingSites = []
for record in SeqIO.parse(datapath + "common_primer.mfasta", "fasta"):
    recordEnding = record[-15:-1]+record[-1]
    primerBindingSites.append([record,recordEnding])


# aus tags_epitropes.mfasta Liste specialTranslatedFeatures erstellen
specialTranslatedFeatures = []
for record in SeqIO.parse(datapath + "tags_epitopes.mfasta", "fasta"):
    specialTranslatedFeatures.append(record)


#Teilaufgabe 1
for plasmid in fiftyRandomPlasmids:
    plasmidlen = len(plasmid)

    searchsequence = plasmid + plasmid
    searchsequencelen = len(searchsequence)

    longseqStart = 0
    longseqEnd = 0
    longseqStartRev = 0
    longseqEndRev = 0


    #Teilaufgabe1

    for prim in primerBindingSites:

        shortseq = prim[1].seq
        shortseqRev = prim[1].reverse_complement().seq
        # Länge des Motivs wird als Ausgangspunkt gewählt, d.h. bitte schon richtig übergeben
        shortseqlen = len(shortseq)

        longseq = prim[0].seq
        longseqRev = prim[0].reverse_complement().seq
        longseqlen = len(longseq)

        i = 0
        while i < searchsequencelen:

            if i - plasmidlen < 0:

                #Suche nach Übereinstimmungen auf dem Vorderstrang

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

                #Suche mit dem reversen Komplement, um Postionen auf dem Gegenstrang zu finden

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
    positions = [0,1,2]


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


            #Teilaufgabe 2
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
                        #print(specFeature.seq)
                        #print(pos.seq)
                    j += 1


            newFile2 = open("fileAfterPart2.gb", "w")
            SeqIO.write(fiftyRandomPlasmids, newFile2, "genbank")



            #Teilaufgabe 3

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

                    # nur wenn vorhergehend ein "M" gefunden wurde, das die Überprüfung der Sequenz zulässt,
                    # wird ein offenes Leseraster durch den ersten Stern geschlossen
                    # --> zwei "*" hintereinander haben keinen Einfluss, denn Leseraster wurde nicht geöffnet
                    if openedFrame == bool(1):
                        end = k

                        # das "*" wird nicht mehr dazugezählt, deswegen kann seine Stelle als Ende des Strangs angegeben werden
                        openReadingFrame = strang[begin*3+i:end*3+i]
                        search = bool(0)

                        if (takeOnlyCareOfLongSequences > 49):
                            #vorher überprüfen, ob bereits im gb-file annotiert
                            for anno in existingAnnotations:
                                if (begin*3+i) != anno.location.start and (end*3+i) != anno.location.end:
                                    search = bool(1)


                        # nehme openReadingFrame in Suchliste auf, aber nur
                        # wenn sie a) lang genug ist und b) noch nicht im gb-File annotiert ist
                        if search == bool(1):
                            if strang.name == "vorderstrang":   strand=1
                            if strang.name == "gegenstrang":    strand=-1

                            openReadingFrames.append([openReadingFrame, begin*3+1, end*3+i, strand])
                            print("langes, offenes Leseraster gefunden an Stelle " + str(begin*3+1))

                            # "schließe" die Sequenz wieder --> erst muss ein neuer Anfang gefunden werden,
                            #  der das Überprüfen einer neuen Sequenz zulässt
                            openedFrame = bool(0)
                k += 1
                takeOnlyCareOfLongSequences += 1

            print("Habe alle offenen Leseraster des Plasmids aufgenommen.")


    # alle offenen Leseraster, die bisher nicht annotiert waren, werden nun über die Blastsuche identifiziert
    #blast_result = open("my_blast.xml", "w")
    print("BLAST-Suche startet.")
    for opFr in openReadingFrames:

        # frame ist der SeqRecord der zu suchenden Sequenz
        # opFr[1] bzw. opFr[2] beinhalten die Position des Starts bzw. des Endes der Sequenz
        frame = opFr[0]

        # suche Sequenz über das Internet in der Protein-Datenbank
        result_handle = NCBIWWW.qblast("blastx", "refseq_protein", frame.seq)

        # schreibe die erhaltenen Ergebnisse in die Datei "my_blast.xml"
        #blast_result.write(result_handle.read())
        #blastResult = open("my_blast.xml")

        blastRecord = NCBIXML.read(result_handle)
        print("blastRecord")
        print(blastRecord)

        hit_title = blastRecord.alignments[0].title
        print(hit_title)
        hit_seq = blastRecord.alignments[0].hsps[0].query
        print(hit_seq)

        i=0
        sameletters=0
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