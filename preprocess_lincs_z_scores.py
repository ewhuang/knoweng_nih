### Author: Edward Huang

from collections import OrderedDict
import json
import operator
import time

### This file dumps the processed dictionary from the LINCS data and stores it
### to file.

# Define -infinity
inf = float('-inf')

def mean(lst):
    return float(sum(lst)) / len(lst)

def get_lincs_genes():
    # Extract lincs_genes from all_map.txt, provided by Sheng.
    f = open('./data/all_map.txt', 'r')
    lincs_genes = []
    for i, line in enumerate(f):
        line = line.split()
        lincs_genes += [line[1]]
    f.close()
    return lincs_genes

def get_drugs_and_gene_matrix():
    f = open('./data/lvl4_Stuart_combinedPvalue_diff_normalize_DMSO.txt', 'r')
    drugs, gene_matrix = [], []
    for line in f:
        line = line.split()
        drug, raw_z_scores = line[0], line[1:]
        drugs += [drug]
        gene_matrix += [raw_z_scores]
    f.close()
    return drugs, gene_matrix

def main():
    lincs_genes = get_lincs_genes()

    drugs, gene_matrix = get_drugs_and_gene_matrix()

    # Transpose the gene matrix so we can filter out the lincs_genes.
    gene_matrix = zip(*gene_matrix) # Genes are now rows.
    assert len(gene_matrix) == len(lincs_genes)
    
    # Dictionary: keys = genes, values = lists of z-scores. Length of list is
    # number of drugs.
    gene_dct = OrderedDict({})
    while gene_matrix != []:
        gene_z_scores = gene_matrix.pop(0)
        assert len(gene_z_scores) == len(drugs)
        gene = lincs_genes.pop(0)
        # Skip invalid genes.
        if gene == '-666':
            continue
        # Convert the list of z-scores to floats.
        gene_z_scores = [inf if score == '#NAME?' else abs(float(
            score)) for score in gene_z_scores]
        if gene not in gene_dct:
            gene_dct[gene] = gene_z_scores
        else:
            # If a gene appears twice, then we want to get the max values.
            for ci, value in enumerate(gene_z_scores):
                gene_dct[gene][ci] = max(gene_dct[gene][ci], value)
    
    # Ensure that we popped off every gene.
    assert len(lincs_genes) == 0
    # Update lincs_genes to be just the valid ones in our dictionary.
    lincs_genes = gene_dct.keys()

    # Transposing back the z-scores.
    drug_matrix = gene_dct.values()[:]
    gene_dct.clear()
    drug_matrix = zip(*drug_matrix)

    # Write out to file.
    out = open('./data/processed_lincs_normalized_DMSO.txt', 'w')
    out.write('\t'.join(lincs_genes) + '\n')
    for i, z_scores in enumerate(drug_matrix):
        out.write('%s\t' % (drugs[i]))
        z_scores = map(str, z_scores)
        assert len(z_scores) == len(lincs_genes)
        out.write('\t'.join(z_scores) + '\n')
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))