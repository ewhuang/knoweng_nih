### Author: Edward Huang

from collections import OrderedDict
import file_operations
import numpy as np
import operator
from scipy.stats import fisher_exact
import sys
import time

### Gets the top pathways for each drug/cell-line using the LINCS data set.
### Usage: python top_pathways_lincs.py
### Run time: 50 minutes

# This variable doesn't affect the code. It is just helps count low p-values
# for eyeballing purposes.
LOW_P_THRESHOLD = 0.0001 # Count how many pathway-drug pairs are below this.

def get_gene_to_z_dct(gene_matrix, drugs, lincs_genes):
    '''
    Returns a dictionary.
    Key: gene -> str
    Value: list of z-scores -> list(float)
    '''
    gene_to_z_dct, gene_counter = OrderedDict({}), 0

    while gene_matrix != []:
        gene_z_scores = gene_matrix.pop(0)
        assert len(gene_z_scores) == len(drugs)
        gene = lincs_genes[gene_counter]

        gene_counter += 1
        # Skip undefined genes.
        if gene == '-666':
            continue

        # Take the absolute value of the z-scores.
        neg_inf = float('-inf')
        gene_z_scores = [neg_inf if score == '#NAME?' else abs(float(score)
            ) for score in gene_z_scores]
        if gene not in gene_to_z_dct:
            gene_to_z_dct[gene] = gene_z_scores
        else:
            # If a gene appears twice, then we want to get the max values.
            for ci, value in enumerate(gene_z_scores):
                gene_to_z_dct[gene][ci] = max(gene_to_z_dct[gene][ci], value)
    return gene_to_z_dct

def get_drug_to_z_dct(drugs, gene_to_z_dct):
    '''
    Returns a dictionary mapping drugs to their LINCS z-scores.
    Key: drug -> str
    Value: list of LINCS z-scores -> list(float)
    '''
    # Transposing back the z-scores.
    drug_to_z_dct = [drugs]
    for gene in gene_to_z_dct:
        drug_to_z_dct += [gene_to_z_dct[gene]]
    gene_to_z_dct.clear()
    drug_to_z_dct = zip(*drug_to_z_dct)

    # Make a new dictionary and reconstruct the drug to z-score mappings.
    temp_drug_to_z_dct = OrderedDict({})
    for row in drug_to_z_dct:
        drug, z_scores = row[0], row[1:]

        assert type(z_scores[0]) == float

        if drug in temp_drug_to_z_dct:
            temp_drug_to_z_dct[drug] += [z_scores]
        else:
            temp_drug_to_z_dct[drug] = [z_scores]
    return temp_drug_to_z_dct

def get_drug_to_top_genes_dct(drug_to_z_dct, lincs_genes):
    '''
    Takes the drug to z-score dictionary and returns a new dictionary. The
    higher the z-score, the more correlated the drug and the gene are.
    Key: drug -> str
    Value: list of genes with the highest z-scores, maximum max_genes_per_drug
                -> list(str)
    '''
    drug_top_genes_dct = {}
    for drug in drug_to_z_dct:
        z_scores = drug_to_z_dct[drug]
        # Average the z-scores for different experiments of the same drug-cell
        # line pair.
        z_scores = [np.mean(x) for x in zip(*z_scores)]

        top_gene_indices = []
        for i, z_score in enumerate(z_scores):
            if z_score >= z_score_min:
                top_gene_indices += [(i, z_score)]
        # Take only at most max_genes_per_drug.
        top_gene_indices = sorted(top_gene_indices, key=lambda k:k[1],
            reverse=True)[:max_genes_per_drug]
        # Update the drug with the top genes from LINCS.
        drug_top_genes_dct[drug] = [lincs_genes[i] for (i,
            z_score) in top_gene_indices]
    return drug_top_genes_dct

def compute_top_pathways_per_drug(lincs_genes, drug_top_genes_dct):
    '''
    Calculates the Fisher's exact test for every drug's top genes and every
    NCI pathway's genes. Returns the sorted results and the number of drug-path
    pairs that exceed LOW_P_THRESHOLD.
    '''
    # Read NCI pathway dictionary and the lincs_genes in the NCI pathways.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()

    # Fisher's test for every drug/cell-line and path pair.
    gene_universe = len(nci_genes.union(lincs_genes))
    fisher_p_val_dct, num_low_p = {}, 0
    for drug in drug_top_genes_dct:
        # Genes with the top z-scores for each drug, up to max_genes_per_drug.
        corr_genes = set(drug_top_genes_dct[drug])
        for path in nci_path_dct:
            path_genes = set(nci_path_dct[path])

            # Get the four relevant numbers for Fisher's test.
            corr_and_path = len(corr_genes.intersection(path_genes))
            corr_not_path = len(corr_genes.difference(path_genes))
            path_not_corr = len(path_genes.difference(corr_genes))
            neither = gene_universe - len(corr_genes.union(path_genes))
            # Compute Fisher's test.
            f_table = ([[corr_and_path, corr_not_path], [path_not_corr,
                neither]])
            o_r, p_value = fisher_exact(f_table, alternative='greater')
            if p_value < LOW_P_THRESHOLD:
                num_low_p += 1
            fkey = (drug, path, corr_and_path, corr_not_path, path_not_corr)
            fisher_p_val_dct[fkey] = p_value
    sorted_fisher = sorted(fisher_p_val_dct.items(), key=operator.itemgetter(1))

    return sorted_fisher, num_low_p

def write_top_pathways(sorted_fisher, num_low_p):
    # Write how many drug-pathway pairs have a very low p-value, determined by
    # LOW_P_THRESHOLD.
    subfolder = './results/lincs_top_pathway_files'
    out = open('%s/top_pathways_lincs_z%g_max%d.txt' % (subfolder, z_score_min,
        max_genes_per_drug), 'w')
    out.write('num_p_below_%s\t%d\n' % (str(LOW_P_THRESHOLD), num_low_p))
    # Write the drug's top pathways to file.
    out.write('drug\tcell_line\tpath\tp-value\tinter\tlincs_len\tpath_len\n')
    for fkey, p_val in sorted_fisher:
        drug, path, inter, corr_len, path_len = fkey
        drug, cell_line = drug.split('_')
        out.write('%s\t%s\t%s\t%g\t%d\t%d\t%d\n' % (drug, cell_line, path,
            p_val, inter, corr_len, path_len))
    out.close()

def main():
    if (len(sys.argv) < 3):
        print "Usage: " + sys.argv[0] + " z_score_min max_genes_per_drug"
        exit(1)

    global z_score_min, max_genes_per_drug
    z_score_min, max_genes_per_drug = float(sys.argv[1]), int(sys.argv[2])

    lincs_genes = file_operations.get_lincs_genes()
    drugs, gene_matrix = file_operations.get_drugs_and_gene_matrix()

    # Transpose the gene matrix so we can filter out the lincs_genes.
    gene_matrix = zip(*gene_matrix)
    assert len(gene_matrix) == len(lincs_genes)
    
    gene_to_z_dct = get_gene_to_z_dct(gene_matrix, drugs, lincs_genes)

    # Update lincs_genes to be just the valid ones in our dictionary. We lose
    # some genes because some are labeled as '-666'
    lincs_genes = gene_to_z_dct.keys()

    drug_to_z_dct = get_drug_to_z_dct(drugs, gene_to_z_dct)

    drug_top_genes_dct = get_drug_to_top_genes_dct(drug_to_z_dct, lincs_genes)

    sorted_fisher, num_low_p = compute_top_pathways_per_drug(lincs_genes,
        drug_top_genes_dct)

    write_top_pathways(sorted_fisher, num_low_p)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))