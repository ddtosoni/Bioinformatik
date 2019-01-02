#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plasmidannotation

-- part of member 3 --

@author: daphne
"""


#TODO Stellen sind jetzt immer ab 0! Evtl. Änderungsbedarf


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
from random import *

# on Linux:
# datapath = "/home/daphne/Documents/fhnw/Bioinformatik/Projekt_Plasmidannotation/daten_plasmid_annotation/"

# on Windows:
datapath = "C:/Users/daphne/Documents/fhnw/Bioinformatik/Projekt_Plasmidannotation/datenPlasmidannotation/"

"""
TODO
- Encoding
- Complementbildung --> dazu muss es ein Se-Objekt sein, oder?
"""



"""
import sys
reload(sys)
sys.setdefaultencoding('utf-16')
"""
#encoding='cp1252'


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
import codecs
        
""" Versuch für Linux, die Datei korrekt auszulesen
# with codecs.open(datapath + "common_primer.mfasta", encoding='utf-8', errors='ignore') as test:
with codecs.open(datapath + "common_primer.mfasta", 'r', encoding='utf-16-le') as test:
    primerBindingSites = []
    print("hello")
#    for record in np.loadtxt(test, skiprows=8):
#        print("die erset Schleife geht! :)")
    for record in SeqIO.parse(test, "fasta"):
        print("do you reach this?")
        primerBindingSites.append(record)
        print(record.id)
"""

# primerBindingSites ist ein SeqRecord; Beweis: print(primerBindingSites)
primerBindingSites = []
# da die Sequenz nur auf den Inhalt der letzten 15 Basen überprüft werden muss,
# werden die Enden separat in der Liste primerBindingEndings gespeichert
primerBindingSiteEndings = []
#print("\nPrimer binding sites:")
for record in SeqIO.parse(datapath + "common_primer.mfasta", "fasta"):
    primerBindingSites.append(record)
    recordEnding = record[-15:-1]+record[-1]
    primerBindingSiteEndings.append(recordEnding)
    #print(record.seq)
    #print(recordEnding.seq)


# aus tags_epitropes.mfasta Liste specialTranslatedFeatures erstellen
#specialTranslatedFeatures ist ein SeqRecord; Beweis: print(specialTranslatedFeatures)
specialTranslatedFeatures = []
#print("\nSpecial translated features:")
for record in SeqIO.parse(datapath + "tags_epitopes.mfasta", "fasta"):
    specialTranslatedFeatures.append(record)
    #print(record.seq)
newfeature = SeqRecord(Seq("ILLVCF"), id="Testfeature")
print(newfeature.seq)
specialTranslatedFeatures.append(newfeature)


""" Teilaufgabe 1 """
#Erkennung funktioniert
#todo später Funktion daraus machen
for plasmid in fiftyRandomPlasmids:
    plasmidlen = plasmid.__len__()
    searchsequence = plasmid + plasmid
    for prim in primerBindingSiteEndings:

        seq = prim.seq
        seqRev = prim.reverse_complement().seq
        # Achtung, Länge des Motivs wird als Ausgangspunkt gewählt, d.h. bitte schon richtig übergeben
        seqlen = prim.__len__()

        #seq = "CTCTGGCC" #wäre am Übergang des zirkulären Plasmids
        #seqlen = 8
        #todo brauche für Verallgemeinerung noch motifRevComp, das überprüft wird
        #print("start")
        i = 0
        while i < searchsequence.__len__():
            if searchsequence[i:i+seqlen].seq == seq:
                #print("YES! Found it.")
                if i -  plasmidlen < 0:
                    print("gefunden an Stelle: " + str(i) + " im Vorderstrang")
            if searchsequence[i:i+seqlen].seq == seqRev:
                if i - plasmidlen < 0:
                    print("gefunden an Stelle: " + str(plasmidlen - (i + seqlen)) + " im Gegenstrang")
            i += 1

    #TODO: in gb-file schreiben, schauen, ob schon angegeben (in diesem Fall nicht dazuschreiben)
print("ready")



"""
specialTranslatedFeatures
die Sequenz in Aminosäuren umwandeln
dafür alle 6 Varianten durchgehen
merke, auf welchem Strang man ist
pro Variante Vergleich der Peptidsequenz mit der Variante der Auslesung
falls identisch auf doppelter Sequenz ergänze im gb-File mit Angabe des Strangs
-  können hier auch schon Eingaben existieren? --> dann Extramethode?
wichtig: genaue Stellenbeschreibung!
"""

"""
proteinCodingGenes
auslesen, welche Stellen annotiert sind
speichere Stellen, die noch nicht annotiert sind, auf doppelter Sequenz
wenn diese größer 50 (größer 200) sind, lese Sequenz aus
suche Sequenz via blastx oder blastp in der refseq_protein Datenbank
lese erstes Ergebnis aus
ergänze im gb-File die Beschreibung des besten Ergebnisses
??? Stellen abspecken oder
- stern ist stopcodon
"""




for plasmid in fiftyRandomPlasmids:
    vorderstrang = plasmid + plasmid
    gegenstrang = vorderstrang.reverse_complement()
    straenge = [vorderstrang, gegenstrang]
    positions = [0,1,2]
    plasmidlen = plasmid.__len__()

    openReadingFrames = []

    for strang in straenge:
        for i in positions:
            print(i)
            pos = strang[i:]
            transPlas = pos.translate()
            print(transPlas.seq)
            transPlaslen = transPlas.__len__()

            """ Teilaufgabe 2 """
            #todo, hab noch nicht ganz verstanden: muss es ein offenes Leseraster sein, in dem ich nach den specialTranslatedFeatures schaue? Eigentlich nicht...
            for specFeature in specialTranslatedFeatures:
                specFeaturelen = specFeature.__len__()
                j = 0
                while j < transPlaslen:
                    pos = transPlas[j:(j+specFeaturelen)]
                    if specFeature.seq == pos.seq:
                        print(plasmid.__len__())
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
                        #Test, ob die Stellenangaben mit der Übersetzung in die Aminosäuren stimmen
                        #print(transPlas[begin:end].seq)
                        #print(openReadingFrame.translate().seq)
                        #TODO überprüfe, ob dort schon etwas annotiert ist
                        #todo: später 199 auf 49 setzen
                        if (takeOnlyCareOfLongSequences > 49):
                            #TODO: vorher überprüfen, ob bereits im gb-file annotiert
                            openReadingFrames.append(openReadingFrame)
                            print("got it! :)")
                            #TODO: genügend große offene Leseraster werden gefunden, aber das Auslesen der Blast-Suche-Ergebnisse funktioniert noch nicht wie gewollt
                            result_handle = NCBIWWW.qblast("blastx", "refseq_protein", openReadingFrame.seq)
                            blast_records = NCBIXML.read(result_handle) #Ergebnis: <Bio.Blast.Record.Blast object at 0x0000028E25C88208>
                            #blast_record = next(blast_records)
                            #print(blast_record)

                            for alignment in blast_records.alignments:
                                print("alignment")
                                for hsp in alignment.hsps:
                                    print("hsp")
                                    if hsp.expect < 0.01:
                                        print("****Alignment****")
                                        print("sequence: " + alignment.title)
                                        print(hsp.query[0:75] + "...")

                            #mit newtry wird die ganze xml-Datei ausgelesen; möchte aber nur den ersten Hit
                            #newtry = result_handle.read()
                            #print(newtry.title())

                        else: print("sequence is to short")
                        openedFrame = bool(0)
                    else:
                        print("counldn't find a new point to start for opening the reading frame")
                k += 1
                takeOnlyCareOfLongSequences += 1
























"""
#Test, ob die Stelle richtig benannt wird im Gegenstrang
for plasmid in fiftyRandomPlasmids:
    plasmidlen = plasmid.__len__()
    searchsequence = plasmid + plasmid
    seq = "AGAAATGT" #"CTCTGGCC" wäre am Übergang des zirkulären Plasmids
    seqRev = "ACATTTCT" #"GGCCAGAG"
    seqlen = 8
    #todo brauche für Verallgemeinerung noch motifRevComp, das überprüft wird
    print("start")
    i = 0
    while i < searchsequence.__len__():
        if searchsequence[i:i+seqlen].seq == seq:
            if i -  plasmidlen < 0:
                print("gefunden an Stelle: " + str(i) + " im Vorderstrang")
        if searchsequence[i:i+seqlen].seq == seqRev:
            if i - plasmidlen < 0:
                print("gefunden an Stelle: " + str(plasmidlen - (i + seqlen)) + " im Gegenstrang")
        i += 1
print("ready")
"""


"""
        print(plasmid.id)
        print(plasmid.seq)
        beginnSequenzVorstrang = plasmid[0:15]
        print(beginnSequenzVorstrang.seq)
        endSequenzVorstrang = plasmid[-15:-1] + plasmid[-1]
        #TODO brauche das Komplement, nicht das reverse Komplement
        beginnSequenzGegenstrang = beginnSequenzGegenstrang.reverse_complement()
        print(beginnSequenzGegenstrang.seq)
        print(plasmid.reverse_complement().seq)
        # Zeile 76 darf ich noch nicht... herausfinden! :)
        plasmidstring = plasmid.seq
        print(plasmidstring)
        print(plasmid.features)


print(fiftyRandomPlasmids)
"""

"""Ergebnis der Blastsuche, mögliche Parameter

<Hit>
  <Hit_num>1</Hit_num>
  <Hit_id>gi|1325999206|ref|WP_101865506.1|</Hit_id>
  <Hit_def>hypothetical protein [Klebsiella quasipneumoniae]</Hit_def>
  <Hit_accession>WP_101865506</Hit_accession>
  <Hit_len>642</Hit_len>
  <Hit_hsps>
    <Hsp>
      <Hsp_num>1</Hsp_num>
      <Hsp_bit-score>35.8094</Hsp_bit-score>
      <Hsp_score>81</Hsp_score>
      <Hsp_evalue>2.12394</Hsp_evalue>
      <Hsp_query-from>10</Hsp_query-from>
      <Hsp_query-to>147</Hsp_query-to>
      <Hsp_hit-from>496</Hsp_hit-from>
      <Hsp_hit-to>542</Hsp_hit-to>
      <Hsp_query-frame>1</Hsp_query-frame>
      <Hsp_hit-frame>0</Hsp_hit-frame>
      <Hsp_identity>21</Hsp_identity>
      <Hsp_positive>29</Hsp_positive>
      <Hsp_gaps>5</Hsp_gaps>
      <Hsp_align-len>49</Hsp_align-len>
      <Hsp_qseq>YFGRPYEGIAVFDG---KKITVTGTLWNGNKIIDERLINPDGSLLFRVT</Hsp_qseq>
      <Hsp_hseq>YYGRAYKQ-AVIDGIKPKAITPKGVIWNGNKVTVQFEV-PNPPLVFDTT</Hsp_hseq>
      <Hsp_midline>Y+GR Y+  AV DG   K IT  G +WNGNK+  +  + P+  L+F  T</Hsp_midline>
    </Hsp>
  </Hit_hsps>
</Hit>
"""


""" StackOverflow Beispiel, wie man open reading frames findet

from Bio.Seq import Seq

seq = Seq("CCTCAGCGAGGACAGCAAGGGACTAGCCAGGAGGGAGAACAGAAACTCCAGAACATCTTGGAAATAGCTCCCAGAAAAGCAAGCAGCCAACCAGGCAGGTTCTGTCCCTTTCACTCACTGGCCCAAGGCGCCACATCTCCCTCCAGAAAAGACACCATGAGCACAGAAAGCATGATCCGCGACGTGGAACTGGCAGAAGAGGCACTCCCCCAAAAGATGGGGGGCTTCCAGAACTCCAGGCGGTGCCTATGTCTCAGCCTCTTCTCATTCCTGCTTGTGGCAGGGGCCACCACGCTCTTCTGTCTACTGAACTTCGGGGTGATCGGTCCCCAAAGGGATGAGAAGTTCCCAAATGGCCTCCCTCTCATCAGTTCTATGGCCCAGACCCTCACACTCAGATCATCTTCTCAAAATTCGAGTGACAAGCCTGTAGCCCACGTCGTAGCAAACCACCAAGTGGAGGAGCAGCTGGAGTGGCTGAGCCAGCGCGCCAACGCCCTCCTGGCCAACGGCATGGATCTCAAAGACAACCAACTAGTGGTGCCAGCCGATGGGTTGTACCTTGTCTACTCCCAGGTTCTCTTCAAGGGACAAGGCTGCCCCGACTACGTGCTCCTCACCCACACCGTCAGCCGATTTGCTATCTCATACCAGGAGAAAGTCAACCTCCTCTCTGCCGTCAAGAGCCCCTGCCCCAAGGACACCCCTGAGGGGGCTGAGCTCAAACCCTGGTATGAGCCCATATACCTGGGAGGAGTCTTCCAGCTGGAGAAGGGGGACCAACTCAGCGCTGAGGTCAATCTGCCCAAGTACTTAGACTTTGCGGAGTCCGGGCAGGTCTACTTTGGAGTCATTGCTCTGTGAAGGGAATGGGTGTTCATCCATTCTCTACCCAGCCCCCACTCTGACCCCTTTACTCTGACCCCTTTATTGTCTACTCCTCAGAGCCCCCAGTCTGTATCCTTCTAACTTAGAAAGGGGATTATGGCTCAGGGTCCAACTCTGTGCTCAGAGCTTTCAACAACTACTCAGAAACACAAGATGCTGGGACAGTGACCTGGACTGTGGGCCTCTCATGCACCACCATCAAGGACTCAAATGGGCTTTCCGAATTCACTGGAGCCTCGAATGTCCATTCCTGAGTTCTGCAAAGGGAGAGTGGTCAGGTTGCCTCTGTCTCAGAATGAGGCTGGATAAGATCTCAGGCCTTCCTACCTTCAGACCTTTCCAGATTCTTCCCTGAGGTGCAATGCACAGCCTTCCTCACAGAGCCAGCCCCCCTCTATTTATATTTGCACTTATTATTTATTATTTATTTATTATTTATTTATTTGCTTATGAATGTATTTATTTGGAAGGCCGGGGTGTCCTGGAGGACCCAGTGTGGGAAGCTGTCTTCAGACAGACATGTTTTCTGTGAAAACGGAGCTGAGCTGTCCCCACCTGGCCTCTCTACCTTGTTGCCTCCTCTTTTGCTTATGTTTAAAACAAAATATTTATCTAACCCAATTGTCTTAATAACGCTGATTTGGTGACCAGGCTGTCGCTACATCACTGAACCTCTGCTCCCCACGGGAGCCGTGACTGTAATCGCCCTACGGGTCATTGAGAGAAATAA")


table = 1
min_pro_len = 100

for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
    for frame in range(3):
        for pro in nuc[frame:].translate(table).split("*"):
            if len(pro) >= min_pro_len:
                print "%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame)

"""