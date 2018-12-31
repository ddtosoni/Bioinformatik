from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
import pandas as pd
from collections import Counter
import numpy as np
from collections import defaultdict
from itertools import chain
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqFeature

'''File to read'''
inputFile = "vectors-100.gb"


'''Initialise dictionary'''
featureDict = dict()


# Define interested tags!
interested_only_note = ['promoter', 'oriT','rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb',
                       'LTR', 'enhancer']  # note

interested_note_and_gen = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']  #note and gene

interested_gene_and_product = ['CDS']  # gene and product

interested_note_and_bound_moiety = ['protein_bind', 'misc_binding']  # note and bound moiety

interested_note_and_mobile_element_type = ['mobile_element']  # note and mobile element type

interested_gene = ['mRNA']  # gene

interested_product = ['tRNA', 'rRNA']  # product


#test_dict = SeqIO.to_dict(SeqIO.parse("vectors-100.gb", "genbank"))


record = SeqIO.parse(inputFile, "genbank")

#[featureDict.__setitem__(item, 1+featureDict.get(item, 1)) for item in interested_only_note]

x = []
testDict = dict()



for search in record:
    if len(search.seq) >= 1500:
        for feature in search.features:
            if feature.type in interested_only_note:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'])
                    a = feature.type + " " + feature.location.extract(search).seq
                    x.append(a)






                except:
                    break
                    #print("Couldn't finde any gene or note in ", interested_only_note)

            if feature.type in interested_note_and_gen:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                     feature.qualifiers['gene'])


                except:
                    break
                    #print("Couldn't find any gene or note in", interested_note_and_gen)

            if feature.type in interested_gene_and_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'],
                                                                 feature.qualifiers['gene'])
                    a = feature.type + " " + feature.location.extract(search).seq
                    x.append(a)

                except:
                    break
                    #print ("Couldn't find any product or gene in ", interested_gene_and_product)

            if feature.type in interested_note_and_bound_moiety:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                 feature.qualifiers['bound_moiety'])
                    a = feature.type + " " + feature.location.extract(search).seq
                    x.append(a)

                except:
                    break
                    #print ("Coulnd't find any note or bound_moiety in ", interested_note_and_bound_moiety)

            if feature.type in interested_gene:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['gene'])
                    a = feature.type + " " + feature.location.extract(search).seq
                    x.append(a)
                except:
                    break
                    #print("Couldn't find any gene in ", interested_gene)

            if feature.type in interested_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'])
                    a = feature.type + " " + feature.location.extract(search).seq
                    x.append(a)
                except:
                    break
                    #print("Couldn't find any product in ", interested_product)
print(featureDict)
print(len(featureDict))

print(len(testDict))

c = Counter(x)
print(c)
print(len(c))

