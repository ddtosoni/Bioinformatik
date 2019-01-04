from Bio import SeqIO
from collections import Counter


'''File to read'''
inputFile = "vectors-100.gb"

'''Initialise dictionary'''
featureDict = dict()

'''Initialise Containers'''
storage = []

# Define interested tags!
interested_only_note = ['promoter', 'oriT', 'rep_origin', 'primer_bind', 'terminator', 'misc_signal', 'misc_recomb',
                                    'LTR', 'enhancer']  # note

interested_note_and_gen = ['-35_signal', '-10_signal', 'RBS', 'polyA_signal', 'sig_peptide']  # note and gene

interested_gene_and_product = ['CDS']  # gene and product

interested_note_and_bound_moiety = ['protein_bind', 'misc_binding']  # note and bound moiety

interested_note_and_mobile_element_type = ['mobile_element']  # note and mobile element type

interested_gene = ['mRNA']  # gene

interested_product = ['tRNA', 'rRNA']  # product

# test_dict = SeqIO.to_dict(SeqIO.parse("vectors-100.gb", "genbank"))


record = SeqIO.parse(inputFile, "genbank")

for search in record:
    if len(search.seq) >= 1500:
        for feature in search.features:
            if feature.type in interested_only_note:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'])
                    value = feature.type + " " + feature.location.extract(search).seq
                    storage.append(value)

                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_note_and_gen:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['gene'])

                    value = feature.type + " " + feature.location.extract(search).seq
                    storage.append(value)

                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_gene_and_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'],
                                                                         feature.qualifiers['gene'])

                    value = feature.type + " " + feature.location.extract(search).seq
                    storage.append(value)

                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_note_and_bound_moiety:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['note'],
                                                                         feature.qualifiers['bound_moiety'])

                    value = feature.type + " " + feature.location.extract(search).seq
                    storage.append(value)

                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_gene:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['gene'])
                    value = feature.type + " " + feature.location.extract(search).seq
                    storage.append(value)

                except:  # (RuntimeError, TypeError, NameError):
                    break

            if feature.type in interested_product:
                try:
                    featureDict[feature.location.extract(search).seq] = (feature.type, feature.qualifiers['product'])
                    value = feature.type + " " + feature.location.extract(search).seq
                    storage.append(value)

                except:  # (RuntimeError, TypeError, NameError):
                    break

print(featureDict)
print(len(featureDict))


CounterInformation = Counter(storage)
print(CounterInformation)
print(len(CounterInformation))

