### Author: Edward Huang

import file_operations
import numpy as np
import operator
from scipy.stats.mstats import kruskalwallis
import time

### Find the top genes for each drug from gene expression.
### Uses Kruskal-Wallis instead of Fisher's test to then use the gene rankings
### to find most similar pathways.
### Run time: 3.4 hours.

# This variable is for counting purposes, only.
KRUSKAL_P_THRESH = 0.05

def get_gene_drug_correlations(drug, drug_response_vector, gene_expression_dct):
    '''
    Returns a dictionary.
    Each key is a (gene, drug, Pearson correlation) triple.
    Each value is the p-value corresponding to the correlation in the key.
    '''
    # Indices of None values in our drug response table.
    NA_i = [i for i, e in enumerate(drug_response_vector) if e == None]
    NA_i.reverse()

    gene_map_list = gene_expression_dct.keys()
    gene_expression_mat = gene_expression_dct.values()
    # Delete None elements in the drug response and gene expression vectors.
    for index in NA_i:
        del drug_response_vector[index]
        for i, row in enumerate(gene_expression_mat):
            del row[index]

    if len(drug_response_vector) == 0:
        return {}

    r = file_operations.generate_correlation_map(np.array(
        [drug_response_vector]), np.array(gene_expression_mat))
    num_valid_patients = len(drug_response_vector)
    p = [file_operations.compute_p_val(pcc, num_valid_patients) for pcc in r]

    gene_drug_correlation_dct = {}
    for i, gene in enumerate(gene_map_list):
        pcc, p_value = r[i], p[i]
        # We keep all p-values here because KW needs them. correlation_top_
        # pathways.py has a threshold because of Fisher's test.
        gene_drug_correlation_dct[(gene, drug, pcc)] = p_value
    return gene_drug_correlation_dct

def get_drug_path_kw_dct(drug, nci_path_dct, corr_gene_dct,
    gene_expression_dct):
    '''
    Returns a (dict, int) pair.
    Key: (drug, path) tuple -> (str, str)
    Value: p-value corresponding to Kruskal-Wallis test -> float

    The int counts the number of drug-pathway pairs that have p-value below 
    KRUSKAL_P_THRESH.
    '''
    drug_path_p_values, num_low_p = {}, 0

    # Compute the top pathways for each drug with the Kruskal-Wallis test.
    for path in nci_path_dct:
        # Make a copy of drug-gene correlations and then remove the genes
        # that are in each pathway.
        non_pathway_correlations = corr_gene_dct.copy()
        pathway_correlations = []
        path_genes = nci_path_dct[path]
        # Split the correlations into genes in the pathways and genes that
        # aren't.
        for gene in path_genes:
            # Skip genes in NCI pathways but not in gene expression data.
            if gene not in gene_expression_dct:
                continue
            pathway_correlations += [non_pathway_correlations[gene]]
            del non_pathway_correlations[gene]

        non_pathway_correlations = non_pathway_correlations.values()
        h_stat, p_value = kruskalwallis(pathway_correlations,
            non_pathway_correlations)
        drug_path_p_values[(drug, path, h_stat)] = p_value

        if p_value < KRUSKAL_P_THRESH:
            num_low_p += 1

    return drug_path_p_values, num_low_p

def write_gene_drug_correlations(gene_drug_correlations):
    '''
    Write to file all significant gene-drug correlations and their p-values.
    '''
    # Sort the top genes by value. Get the top genes.
    gene_out = open('./results/correlation_top_genes_kw.txt', 'w')
    gene_out.write('gene\tdrug\tcorrelation\tp_value\n')
    # gene_out.write('gene\tdrug\tcorrelation\n')
    gene_drug_correlations = sorted(gene_drug_correlations.items(),
        key=operator.itemgetter(1))
    for (gene, drug, pcc), p_value in gene_drug_correlations:
        gene_out.write('%s\t%s\t%f\t%g\n' % (gene, drug, pcc, p_value))
    gene_out.close()

def write_drug_path_correlations(drug_path_p_values, num_low_p):
    '''
    Writes out all drug-pathway Pearson p-values, sorted. Also writes the number
    of significant drug-pathway correlations in the first line.
    '''
    top_paths = sorted(drug_path_p_values.items(), key=operator.itemgetter(1))

    # Write out the results.
    path_out = open('./results/correlation_top_pathways_kw.txt', 'w')
    path_out.write('num_below_%f\t%d\n' % (KRUSKAL_P_THRESH, num_low_p))
    path_out.write('drug\tpath\tp_value\th_statistic\n')    
    for (drug, path, h_stat), p_val in top_paths:
        path_out.write('%s\t%s\t%g\t%g\n' % (drug, path, p_val, h_stat))
    path_out.close()

def write_top_pathways():
    '''
    The main function.
    '''
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
    drug_response_dct = file_operations.get_drug_response_dct()

    gene_drug_correlations = {}
    # Key: (drug, path, h-statistic). Value: p-value of KW test between drug's
    # correlated genes and path's correlated genes.
    drug_path_p_values = {}
    # Counts the number of drug-pathway pairs with KW p-value KRUSKAL_P_THRESH.
    num_low_p, num_drugs = 0, float(len(drug_response_dct))

    for drug in drug_response_dct:
        # Fetch this repeatedly because of deep copy issues with del.
        gene_expression_dct = file_operations.get_gene_expression_dct()
        
        drug_response_vector = drug_response_dct[drug]

        gene_drug_correlation_dct = get_gene_drug_correlations(drug,
            drug_response_vector, gene_expression_dct)

        if gene_drug_correlation_dct == {}:
            continue

        gene_drug_correlations.update(gene_drug_correlation_dct)

        corr_gene_dct = {}
        for (gene, drug, pcc) in gene_drug_correlation_dct:
            corr_gene_dct[gene] = gene_drug_correlations[(gene, drug, pcc)]

        curr_drug_path_p_values, curr_num_low_p = get_drug_path_kw_dct(drug,
            nci_path_dct, corr_gene_dct, gene_expression_dct)

        drug_path_p_values.update(curr_drug_path_p_values)
        num_low_p += curr_num_low_p

    write_drug_path_correlations(drug_path_p_values, num_low_p)
    write_gene_drug_correlations(gene_drug_correlations)

def main():
    write_top_pathways()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))