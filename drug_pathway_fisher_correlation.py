### Author: Edward Huang

import file_operations
import numpy as np
import operator
from scipy.stats import fisher_exact
import sys
import time

### For each drug, find the k most correlated genes, where k is the command line
### argument. Write this out to file. Then, for each drug, find the pathways
### that have the highest overlap with the correlated genes. Write this out to
### file.
### Run time: 6 minutes.

# FISHER_P_THRESH = 0.0001

def get_gene_drug_pearson_dct(drug, dr_vector, gene_expr_mat, gene_list):
    '''
    Returns a dictionary. Each value is the p-value corresponding to the
    correlation in the key. Only records correlations with p-value below
    pearson_thresh.
    Key: (gene, drug, Pearson correlation) -> (str, str, float)
    Value: p-value -> float
    '''
    # Indices of -666 patients in the drug response.
    good_indices = np.where(dr_vector != -666)[0]

    # Delete -666 drug response patients from both matrices.
    sliced_dr = np.array([dr_vector[good_indices]]) # Convert to matrix.
    sliced_mat = gene_expr_mat[:,good_indices]

    # Compute the correlation between the DR vector and every expression vector.
    r = file_operations.generate_correlation_map(sliced_dr, sliced_mat)
    # Get the corresponding p-value for each Pearson correlation coefficient.
    p = [file_operations.compute_p_val(pcc, len(good_indices)) for pcc in r]

    # Build the dictionary mapping gene-drug pairs to their correlation p-value.
    gene_drug_pearson_dct = {}
    for i, gene in enumerate(gene_list):
        pcc, p_value = r[i], p[i]
        if p_value < pearson_thresh:
            if sort_value == 'sortCorr':
                gene_drug_pearson_dct[(gene, drug)] = pcc
            else:
                gene_drug_pearson_dct[(gene, drug, pcc)] = p_value
    return gene_drug_pearson_dct

# def get_drug_path_fisher_dct(drug, gene_universe, corr_genes, path_to_gene_dct):
#     '''
#     Returns a dictionary mapping each pathway's gene's overlap with the input
#     drug's most correlated genes.
#     Key: (drug, path, Fisher's test table) -> (str, str, float, float, float,
#         float)
#     Value: Fisher's test p-value -> float
#     Also returns the number of drug-pathway pairs that have p-value below 
#     FISHER_P_THRESH.
#     '''
#     drug_path_fisher_dct, num_low_p = {}, 0

#     # Compute the top pathways for the drug with Fisher's test.
#     for path in path_to_gene_dct:
#         path_genes = path_to_gene_dct[path]
#         corr_and_path = len(corr_genes.intersection(path_genes))
#         corr_not_path = len(corr_genes.difference(path_genes))
#         path_not_corr = len(path_genes.difference(corr_genes))
#         neither = gene_universe - len(corr_genes.union(path_genes))

#         # o_r = odds ratio.
#         o_r, p_value = fisher_exact([[corr_and_path, corr_not_path], 
#             [path_not_corr, neither]], alternative='greater')
        
#         # Count the number of significant p-values.
#         if p_value < FISHER_P_THRESH:
#             num_low_p += 1
#         drug_path_fisher_dct[(drug, path, corr_and_path, corr_not_path,
#             path_not_corr, neither)] = p_value

#     return drug_path_fisher_dct, num_low_p

# def write_drug_path_correlations(drug_path_fisher_dct, num_low_p):
#     '''
#     Writes out all drug-pathway Pearson p-values, sorted. Also writes the number
#     of significant drug-pathway correlations in the first line.
#     '''
#     # Sort the drug-pathway overlaps by p-value.
#     fish_dct = sorted(drug_path_fisher_dct.items(), key=operator.itemgetter(1))

#     path_out = open('./results/drug_path_fisher_%d.txt' % top_k, 'w')
#     # Write out header lines.
#     path_out.write('Number of drug-pathway pairs with Fisher\'s p-value < %f'
#         '\t%d\n' % (FISHER_P_THRESH, num_low_p))
#     path_out.write('drug\tpath\tp_value\tinter\tcorr\tpath\tneither\n')
#     for key, p_val in fish_dct:
#         drug, path, inter, corr_len, path_len, neither = key
#         path_out.write('%s\t%s\t%g\t%d\t%d\t%d\t%d\n' % (drug, path, p_val,
#             inter, corr_len, path_len, neither))
#     path_out.close()

def write_gene_drug_pearson_dct(gene_drug_pearson_dct):
    '''
    Write to file all significant gene-drug correlations and their p-values.
    '''
    # Sort the top genes by value. Get the top genes.
    reverse = sort_value == 'sortCorr'
    gene_drug_pearson_dct = sorted(gene_drug_pearson_dct.items(),
        key=operator.itemgetter(1), reverse=reverse)

    gene_out = open('./results/gene_drug_pearson_%s_%g.txt' % (sort_value,
        pearson_thresh), 'w')
    # gene_out.write('gene\tdrug\tpearson_correlation\tp_value\n')
    gene_out.write('gene\tdrug\tpearson_correlation\n')
    # for (gene, drug, pcc), p_value in gene_drug_pearson_dct:
    #     gene_out.write('%s\t%s\t%g\t%g\n' % (gene, drug, pcc, p_value))
    if sort_value == 'sortCorr':
        for (gene, drug), pcc in gene_drug_pearson_dct:
            gene_out.write('%s\t%s\t%g\n' % (gene, drug, pcc))
    else:
        for (gene, drug, pcc), p_val in gene_drug_pearson_dct:
            gene_out.write('%s\t%s\t%g\n' % (gene, drug, pcc))
    gene_out.close()

def main():
    if len(sys.argv) != 4:
        print 'Usage: %s top_k sortCorr/sortP pearson_thresh' % sys.argv[0]
        exit(1)
    global top_k, sort_value, pearson_thresh
    top_k, sort_value = int(sys.argv[1]), sys.argv[2]
    pearson_thresh = float(sys.argv[3])
    assert sort_value in ['sortCorr', 'sortP']

    # Drug response and gene expression information.
    drug_to_dr_dct = file_operations.get_drug_to_dr_dct()
    gene_expr_mat, gene_list = file_operations.get_gene_expr_mat()
    # NCI pathway information.
    path_to_gene_dct, nci_genes = file_operations.get_path_to_gene_dct()
    gene_universe = len(nci_genes.union(gene_list))

    # Number of drug-pathway pairs with Fisher p-value below FISHER_P_THRESH.
    # num_low_p = 0

    # Dictionaries to write out to file.
    gene_drug_pearson_dct, drug_path_fisher_dct = {}, {}
    for drug in drug_to_dr_dct:
        dr_vector = drug_to_dr_dct[drug]

        # Get the all gene correlations to the current drug.
        curr_corr_dct = get_gene_drug_pearson_dct(drug, dr_vector,
            gene_expr_mat, gene_list)

        # Sort correlated genes by increasing p-value, then take top K.
        # TODO: sorting by PCC, rather than p-value.
        # top_genes = sorted(curr_corr_dct.items(), key=operator.itemgetter(1))
        # top_genes = sorted(curr_corr_dct.items(), key=operator.itemgetter(1),
        #     reverse=True)
        # top_genes = top_genes[:top_k]
        # # Keep only the gene names.
        # # top_genes = set([gene for (gene, drug, pcc), p_value in top_genes])
        # top_genes = set([gene for (gene, drug), pcc in top_genes])

        # Compute Fisher's test between pathway genes and current drug's most
        # correlated genes.
        # curr_fish_dct, curr_num_low_p = get_drug_path_fisher_dct(drug,
        #     gene_universe, top_genes, path_to_gene_dct)
        # Update dictionary for all drugs.
        # drug_path_fisher_dct.update(curr_fish_dct)
        # num_low_p += curr_num_low_p

        # Lastly, update the master dictionary for all drugs.
        gene_drug_pearson_dct.update(curr_corr_dct)

    # Write dictionaries out to file.
    # write_drug_path_correlations(drug_path_fisher_dct, num_low_p)
    write_gene_drug_pearson_dct(gene_drug_pearson_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))