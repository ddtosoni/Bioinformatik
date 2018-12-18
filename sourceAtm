from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqFeature

'''File to read'''
inputFile = "vectors-100.gb"


# Define interested tags!
interested_tags = ['promoter', 'CDS', 'polyA_signal','rep_origin', 'primer_bind',
 'terminator', 'protein_bind', 'misc_binding', 'misc_recomb', 'oriT', 'LTR',
 'misc_signal', 'enhancer', 'mobile_element', 'RBS', 'sig_peptide', '-10_signal',
 '-35_signal', 'mRNA', 'tRNA', 'rRNA']


record = SeqIO.parse(inputFile, "genbank")


for search in record:
    if len(search.seq) >= 1500:
        for feature in search.features:
            if feature.type in interested_tags:
                start = feature.location.start.position+1
                end = feature.location.end.position
                d = dict([(feature.type, search.seq[start:end])])
                print(d)

