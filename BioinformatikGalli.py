from Bio import SeqIO
from collections import Counter
from itertools import chain
from collections import defaultdict


'''File to read'''
inputFile = "vectors.gb"

'''Initialise dictionary'''
featureDict = dict()
counterDict = dict()

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

x = []

for search in record:
    if len(search.seq) >= 1500:
        for feature in search.features:
            x.append(feature.type)
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

'''Manage the Dictionary'''
CounterInformation = Counter(storage)

for k in list(CounterInformation):
        if CounterInformation[k] <= 2:
            del CounterInformation[k]
            del featureDict[k]


counterDict.update(CounterInformation)


'''Create correct Dictionary'''
finalDictionary = defaultdict(list)
for a, b in chain(featureDict.items(), counterDict.items()):
    finalDictionary[a].append(b)


for k in list(finalDictionary.keys()):
    if len(k) <= 10:
        del finalDictionary[k]
        del featureDict[k]


print(len(finalDictionary))


a = Counter(x)
print(a)
