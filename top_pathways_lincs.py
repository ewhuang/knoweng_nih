### Author: Edward Huang

import numpy as np
from collections import OrderedDict
from scipy.stats import fisher_exact
import operator
import sys

### Gets the top pathways for each drug/cell-line using the LINCS data set.
### Can read in level 3 or 4 LINCS data.
### Usage: python top_pathways_lincs.py 3/4

Z_SCORE_MIN = 2.5
MAX_GENES_PER_DRUG = 500

### Create new list copy without duplicates.
def create_no_dup(lst):
    new_lst = []
    for e in lst:
        if e not in new_lst:
            new_lst += [e]
    return new_lst

if __name__ == '__main__':
    level = int(sys.argv[1])

    # Translate the drugs into English.
    f = open('./data/drug_translation.txt', 'r')
    trans_dct = {}
    for line in f:
        line = line.split()
        trans_dct[line[1]] = line[0]

    print 'Extracting NCI pathways...'
    path_file = open('./data/nci_pathway.txt', 'r')
    pathnames = []
    nci_path_dct = {}
    nci_genes = set([])
    for line in path_file:
        line = line.split('\t')
        path_name, path_gene = line[0], line[1][:-2]
        nci_genes.add(path_gene)
        if path_name not in nci_path_dct:
            nci_path_dct[path_name] = [path_gene]
        else:
            nci_path_dct[path_name] += [path_gene]
        if path_name not in pathnames:
            pathnames += [path_name]
    path_file.close()

    # Extract genes from old file...
    f = open('./data/lincs_zscore.txt', 'r')
    genes = []
    for i, line in enumerate(f):
        line = line.split()
        if i == 0:
            continue
        else:
            genes += [line[1]]
    f.close()

    print 'Extracting LINCS data...'
    if level == 4:
        f = open('./data/lincs_zscore_new.txt', 'r')
    elif level == 3:
        f = open('./data/lincs_zscore.txt', 'r')
    drugs = []
    gene_dct = OrderedDict({})
    inf = float('-inf')
    for i, line in enumerate(f):
        line = line.split()
        if i == 0:
            drugs = line
        elif genes[i-1] == '-666':
            continue
        else:
            if level == 4:
                gene, coeffs = genes[i-1], line[2:]
            elif level == 3:
                gene, coeffs = genes[i-1], line[3:]
            # #NAME? is negative infinity.
            coeffs = [inf if x == '#NAME?' else abs(float(x)) for x in coeffs]
            if gene not in gene_dct:
                gene_dct[gene] = coeffs
            else:
                # If a gene appears twice, then we want to get the max values.
                curr_vals = gene_dct[gene]
                for ci, value in enumerate(coeffs):
                    gene_dct[gene][ci] = max(curr_vals[ci], value)
    f.close()
    # Delete bad genes from the listself.
    while '-666' in genes:
        genes.remove('-666')

    # Remove the duplicate genes.
    genes = create_no_dup(genes)

    # Take out the lvl4_ prefix in all of the drug strings.
    if level == 4:
        drugs = [raw_string[5:] for raw_string in drugs]

    print 'Cleaning data and converting to dictionary...'
    drug_matrix = [drugs]
    for gene in gene_dct:
        drug_matrix += [gene_dct[gene]]
    drug_matrix = np.array(drug_matrix).transpose()

    # Change names of the drug ID's to English drug names.
    temp_drug_matrix = OrderedDict({})
    for i, row in enumerate(drug_matrix):
        drug, cell_line = row[0].split('_')
        drug, values = trans_dct[drug] + '_' + cell_line, row[1:]
        values = map(float, values)
        if drug in temp_drug_matrix:
            temp_drug_matrix[drug] += [values]
        else:
            temp_drug_matrix[drug] = [values]
    drug_matrix = temp_drug_matrix

    for drug in drug_matrix:
        values = drug_matrix[drug]
        values = [np.mean(x) for x in zip(*values)]
        top_gene_indices = []
        for i, z_score in enumerate(values):
            if z_score > Z_SCORE_MIN:
                top_gene_indices += [(i, z_score)]
        # Take only the top 500 genes at most.
        top_gene_indices = sorted(top_gene_indices, key=lambda k:k[1],
            reverse=True)[:MAX_GENES_PER_DRUG]
        top_genes = [genes[i] for i, z_score in top_gene_indices]
        drug_matrix[drug] = top_genes

    if level == 4:
        outfile = open('./results/top_pathways_lincs_lvl4.txt', 'w')
    elif level == 3:
        outfile = open('./results/top_pathways_lincs_lvl3.txt', 'w')
    # Fisher's test for every drug/cell-line and path pair.
    for drug in drug_matrix:
        fisher_dct = {}
        for path_index, path in enumerate(nci_path_dct):
            path_genes = set(nci_path_dct[path])
            n = len(path_genes)
            corr_genes = set(drug_matrix[drug])
            corr_and_path = len(corr_genes.intersection(path_genes))
            corr_not_path = len(corr_genes.difference(path_genes))
            path_not_corr = len(path_genes.difference(corr_genes))
            neither = len(nci_genes.union(genes)) - len(corr_genes.union(path_genes))
            odds_ratio, p_value = fisher_exact([[corr_and_path, corr_not_path],
                [path_not_corr, neither]])
            fisher_dct[(drug, path, len(corr_genes), len(path_genes), corr_and_path)] = p_value
        sorted_fisher = sorted(fisher_dct.items(), key=operator.itemgetter(1))
        # Write the drug's top pathways to file.
        for pair, score in sorted_fisher:
            drug, path, corr, path_len, inter = pair
            drug, cell_line = drug.split('_')
            combo = drug + '\t' + cell_line + '\t' + path
            outfile.write(combo + '\t' + str(score) + '\t' + str(corr) + '\t' + str(path_len) + '\t' + str(inter) + '\n')
    outfile.close()