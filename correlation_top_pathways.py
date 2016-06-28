### Author: Edward Huang

import file_operations
import operator
from scipy.stats.stats import pearsonr
from scipy.stats import fisher_exact
import time

### Find the top correlated genes from the gene expression data sets.
### Uses Fisher's test with the resulting gene rankings to find most correlated
### pathways.
### Run time: 5.3 hours.

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
    gene_drug_correlation_dct = {}

    # Indices of None values in our drug response table.
    NA_i = [i for i, e in enumerate(drug_response_vector) if e == None]
    # Delete None elements in the drug response vectors.
    drug_response_vector = [e for i, e in enumerate(drug_response_vector
        ) if i not in NA_i]

    if len(drug_response_vector) == 0:
        return {}

    for gene in gene_expression_dct:
        # Remove indices for elements that were None in drug response.
        gene_expression = [e for i, e in enumerate(gene_expression_dct[gene]
            ) if i not in NA_i]

        # Compute the Pearson correlation between these drug response and gene
        # expression.
        pcc, p_value = pearsonr(drug_response_vector, gene_expression)
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
    # Getting gene expression dictionary.
    gene_expression_dct = file_operations.get_gene_expression_dct()
    # Extract the NCI pathway data.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
    # The set of all genes in the gene expression data in addition to NCI data.
    gene_universe = nci_genes.union(gene_expression_dct.keys())
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
    for drug in drug_response_dct:
        print 'Progress: %f%%' % (progress_counter / num_drugs * 100)

        drug_response_vector = drug_response_dct[drug]

        curr_gene_drug_correlation_dct = get_gene_drug_correlations(drug,
            drug_response_vector, gene_expression_dct)

        if curr_gene_drug_correlation_dct == {}:
            continue

        gene_drug_correlations.update(curr_gene_drug_correlation_dct)

        # Get up to MAX_GENES_PER_DRUG correlated genes for the drug.
        corr_genes = sorted(curr_gene_drug_correlation_dct.items(),
            key=operator.itemgetter(1))[:MAX_GENES_PER_DRUG]

        # Get just the genes now, without the corresponding scores.
        corr_genes = set([gene for (gene, drug, pcc), p_value in corr_genes])

        curr_drug_path_p_values, curr_num_low_p = get_drug_path_fisher_dct(drug,
            gene_universe, corr_genes, nci_path_dct)
        drug_path_p_values.update(curr_drug_path_p_values)
        num_low_p += curr_num_low_p

        progress_counter += 1.0

    write_drug_path_correlations(drug_path_p_values, num_low_p)
    write_gene_drug_correlations(gene_drug_correlations)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))