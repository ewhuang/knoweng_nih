### Author: Edward Huang

import cPickle as pickle
import file_operations
import json
import numpy as np
import operator
from scipy.stats import pearsonr, betai
from scipy.stats import fisher_exact
import time

### Find the top correlated genes from the gene expression data sets.
### Uses Fisher's test with the resulting gene rankings to find most correlated
### pathways.
### Run time: 75 minutes.

# The maximum Pearson's p-value between drug response and gene expression to be
# a significantly correlated gene-drug pair.
PEARSON_P_THRESH = 0.005
MAX_GENES_PER_DRUG = 250
FISHER_P_THRESH = 0.0001

def get_gene_drug_correlations(drug, drug_response_vector, gene_expression_dct):
    '''
    Returns a dictionary.
    Each key is a (gene, drug, Pearson correlation) triple.
    Each value is the p-value corresponding to the correlation in the key.
    Only records correlations with p-value below PEARSON_P_THRESH.
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
        if p_value < PEARSON_P_THRESH:
            gene_drug_correlation_dct[(gene, drug, pcc)] = p_value
    return gene_drug_correlation_dct

def get_drug_path_fisher_dct(drug, gene_universe, corr_genes, nci_path_dct):
    '''
    Returns a (dictionary, int) pair.
    Each key is a drug-path tuple along with their corresponding Fisher's test
    table information.
    Each value is the p-value corresponding to the Fisher's test.
    The int counts the number of drug-pathway pairs that have p-value below 
    FISHER_P_THRESH.
    '''
    drug_path_p_values, num_low_p = {}, 0

    # Compute the top pathways for the drug with Fisher's test.
    for path in nci_path_dct:
        path_genes = set(nci_path_dct[path])
        corr_and_path = len(corr_genes.intersection(path_genes))
        corr_not_path = len(corr_genes.difference(path_genes))
        path_not_corr = len(path_genes.difference(corr_genes))
        neither = len(gene_universe) - len(corr_genes.union(path_genes))

        # o_r = odds ratio.
        o_r, p_value = fisher_exact([[corr_and_path, corr_not_path], 
            [path_not_corr, neither]])
        
        # Count the number of significant p-values.
        if p_value < FISHER_P_THRESH:
            num_low_p += 1
        drug_path_p_values[(drug, path, corr_and_path, corr_not_path,
            path_not_corr, neither)] = p_value

    return drug_path_p_values, num_low_p

def write_gene_drug_correlations(gene_drug_correlations):
    '''
    Write to file all significant gene-drug correlations and their p-values.
    '''
    # Sort the top genes by value. Get the top genes.
    gene_out = open('./results/top_genes_correlation_hgnc.txt', 'w')
    gene_out.write('gene\tdrug\tcorrelation\tp_value\n')
    # gene_out.write('gene\tdrug\tcorrelation\n')
    gene_drug_correlations = sorted(gene_drug_correlations.items(),
        key=operator.itemgetter(1))
    # TODO
    for (gene, drug, pcc), p_value in gene_drug_correlations:
        gene_out.write('%s\t%s\t%f\t%g\n' % (gene, drug, pcc, p_value))
    # for (gene, drug), pcc in gene_drug_correlations:
    #     gene_out.write('%s\t%s\t%f\n' % (gene, drug, pcc))
    gene_out.close()

def write_drug_path_correlations(drug_path_p_values, num_low_p):
    '''
    Writes out all drug-pathway Pearson p-values, sorted. Also writes the number
    of significant drug-pathway correlations in the first line.
    '''
    sorted_drug_path_p_values = sorted(drug_path_p_values.items(),
        key=operator.itemgetter(1))

    path_out = open('./results/top_pathways_correlation_hgnc.txt', 'w')
    path_out.write('num_below_%f\t%d\n' % (FISHER_P_THRESH, num_low_p))
    path_out.write('drug\tpath\tp_value\tinter\tcorr\tpath\tneither\n')    
    for key, p_val in sorted_drug_path_p_values:
        drug, path, inter, corr_len, path_len, neither = key
        string = '%s\t%s\t%g\t%d\t%d\t' % (drug, path, p_val, inter, corr_len)
        string += '%d\t%d\n' % (path_len, neither)
        path_out.write(string)
    path_out.close()

def main():
    # Extract the NCI pathway data.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
    # Get the drug response dictionary.
    drug_response_dct = file_operations.get_drug_response_dct()

    # Each key is a (drug, gene, pearson_cc) triple. Each value is the p-value.
    gene_drug_correlations = {}
    # Dictionary of the top drug-pathway pairs.
    drug_path_p_values = {}
    # Number of drug-pathway pairs with Fisher p-value below FISHER_P_THRESH.
    num_low_p = 0

    progress_counter = 0
    num_drugs = float(len(drug_response_dct))
    gene_universe = []

    for drug in drug_response_dct:
        # Fetch this repeatedly because of deep copy issues with del[].
        gene_expression_dct = file_operations.get_gene_expression_dct()
        if gene_universe == []:
            gene_universe = nci_genes.union(gene_expression_dct.keys())
        drug_response_vector = drug_response_dct[drug]

        curr_gene_drug_correlation_dct = get_gene_drug_correlations(drug,
            drug_response_vector, gene_expression_dct)

        if curr_gene_drug_correlation_dct == {}:
            continue

        gene_drug_correlations.update(curr_gene_drug_correlation_dct)

        # Get up to MAX_GENES_PER_DRUG correlated genes for the drug.
        # TODO
        corr_genes = sorted(curr_gene_drug_correlation_dct.items(),
            key=operator.itemgetter(1))[:MAX_GENES_PER_DRUG]
        # corr_genes = sorted(curr_gene_drug_correlation_dct.items(),
        #     key=operator.itemgetter(1), reverse=True)[:MAX_GENES_PER_DRUG]

        # # Get just the genes now, without the corresponding scores.
        # # TODO
        corr_genes = set([gene for (gene, drug, pcc), p_value in corr_genes])
        # corr_genes = set([gene for (gene, drug), pcc in corr_genes])

        curr_drug_path_p_values, curr_num_low_p = get_drug_path_fisher_dct(drug,
            gene_universe, corr_genes, nci_path_dct)
        drug_path_p_values.update(curr_drug_path_p_values)
        num_low_p += curr_num_low_p

        progress_counter += 1.0
        print 'Progress: %f%%' % (progress_counter / num_drugs * 100)

    write_drug_path_correlations(drug_path_p_values, num_low_p)
    write_gene_drug_correlations(gene_drug_correlations)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))