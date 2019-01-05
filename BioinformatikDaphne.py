#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmidannotation

-- part of member 3 --

@author: daphne
"""


#TODO Stellen sind jetzt immer ab 0! Evtl. Änderungsbedarf


from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from random import *
import os

datapath = "C:/Users/daphne/Documents/fhnw/Bioinformatik/Projekt_Plasmidannotation/datenPlasmidannotation/"

fiftyRandomPlasmids = []
i=0
#(TODO): aktuell keine Sicherheit, dass es wirklich 50 werden.. aber das ist gerade nicht so wichtig
for record in SeqIO.parse(datapath + "vectors.gb", "genbank"):
    x = 1
    #x = randint(0, 1)
    if (x == 1 and i < 50):
        fiftyRandomPlasmids.append(record)
        i = i+1
        print(record.id + "\t " + str(i))
    #todo später wegmachen
    break



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


""" Teilaufgabe 1 """
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
        while i < searchsequence.__len__():
            if searchsequence[i:i+seqlen].seq == seq:
                if i -  plasmidlen < 0:
                    print("gefunden an Stelle: " + str(i) + " im Vorderstrang")
            if searchsequence[i:i+seqlen].seq == seqRev:
                if i - plasmidlen < 0:
                    print("gefunden an Stelle: " + str(plasmidlen - (i + seqlen)) + " im Gegenstrang")
            i += 1

    #TODO: in gb-file schreiben, schauen, ob schon angegeben (in diesem Fall nicht dazuschreiben)
print("ready")


for plasmid in fiftyRandomPlasmids:
    vorderstrang = plasmid + plasmid
    gegenstrang = vorderstrang.reverse_complement()
    straenge = [vorderstrang, gegenstrang]
    positions = [0,1,2]
    plasmidlen = len(plasmid)

    openReadingFrames = []

    for strang in straenge:
        for i in positions:
            print(i)
            pos = strang[i:]
            transPlas = pos.translate()
            print(transPlas.seq)
            transPlaslen = len(transPlas)

            """ Teilaufgabe 2 """
            #todo, hab noch nicht ganz verstanden: muss es ein offenes Leseraster sein, in dem ich nach den specialTranslatedFeatures schaue? Eigentlich nicht...
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
                            #f = SeqFeature(FeatureLocation(featureLocation, featureLocation+specFeaturelen*3))
                        else:
                            print("Wiederholungstäter der Stelle " + str(featureLocation-plasmidlen))
                        print(specFeature.seq)
                        print(pos.seq)
                    j += 1


            """Teilaufgabe 3"""

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

                        #todo: später 199 auf 49 setzen
                        if (takeOnlyCareOfLongSequences > 49):
                            #TODO: vorher überprüfen, ob bereits im gb-file annotiert
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

                            os.remove("my_blast.xml")

                        else:
                            print("sequence is to short")
                        openedFrame = bool(0)
                    else:
                        print("counldn't find a new point to start for opening the reading frame")
                k += 1
                takeOnlyCareOfLongSequences += 1