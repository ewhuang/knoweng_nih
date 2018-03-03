### Author: Edward Huang

import argparse
import file_operations
from multiprocessing import Pool
import numpy as np
import operator
import pandas as pd
from scipy.stats import fisher_exact, pearsonr
# import sys
import time

### For each drug, find the k most correlated genes, where k is the command line
### argument. Write this out to file. Then, for each drug, find the pathways
### that have the highest overlap with the correlated genes. Write this out to
### file.
### Run time: 6 minutes.

# FISHER_P_THRESH = 0.0001

def rename_cols(dr_df, gene_df):
    '''
    Gene essentiality samples need to be parsed.
    '''
    # Map old column names to new column names.
    rename_mappings = {}
    for col_name in list(gene_df):
        # Take only the first string before an underscore as sample name.
        rename_mappings[col_name] = col_name.split('_')[0].upper()
    # Rename dataframe with the mappings.
    gene_df.rename(rename_mappings, axis='columns', inplace=True)
    
    # Only keep the samples in common for genes and drugs.
    samples_in_common = list(set(gene_df).intersection(list(dr_df)))

    # Reform the columns for both matrices.
    dr_df = dr_df[samples_in_common]
    gene_df = gene_df[samples_in_common]

    return dr_df, gene_df

# def get_gene_drug_pearson_dct(drug, dr_vector, gene_expr_mat, gene_list):
def get_gene_drug_pearson_dct(drug):
    '''
    Returns a dictionary. Each value is the p-value corresponding to the
    correlation in the key. Only records correlations with p-value below
    pearson_thresh.
    Key: (gene, drug, Pearson correlation) -> (str, str, float)
    Value: p-value -> float
    '''
    # drug = dr_vector.index.values[0]
    dr_vector = dr_df.loc[[drug]]
    num_genes = len(gene_df.index.values)
    # Get list of bad sample names.
    bad_dr_sample_lst = dr_vector.columns[dr_vector.isna().any()].tolist()

    gene_drug_pearson_dct = {}
    for gene, row in gene_df.iterrows():
        bad_ge_sample_lst = [sample for sample in row.axes[0] if pd.isnull(row[sample])]
        bad_samples_both = list(set(bad_dr_sample_lst).union(bad_ge_sample_lst))

        filtered_dr = np.array(dr_vector.drop(bad_samples_both, axis=1))[0]
        filtered_row = np.array(row.drop(bad_samples_both))

        pcc, p_value = pearsonr(filtered_dr, filtered_row)
        p_value *= num_genes
        if p_value < args.pearson_thresh:
            if args.sort_value == 'sortCorr':
                gene_drug_pearson_dct[(gene, drug)] = pcc
            else:
                gene_drug_pearson_dct[(gene, drug, pcc)] = p_value

    #     # bad_ge_sample_lst = row.columns[row.isna().any()].tolist()
    #     # print bad_ge_sample_lst
    # # Drop bad samples from drug response vector.
    # dr_matrix = np.array(dr_vector.drop(bad_dr_sample_lst, axis=1))
    # # Drop bad samples from gene expression matrix.
    # ge_matrix = np.array(gene_expr_df.drop(bad_dr_sample_lst, axis=1))
    # # # Indices of -666 patients in the drug response.
    # # good_indices = np.where(dr_vector != -666)[0]

    # # Delete -666 drug response patients from both matrices.
    # # sliced_dr = np.array([dr_vector[good_indices]]) # Convert to matrix.
    # # sliced_mat = gene_expr_mat[:,good_indices]

    # # Compute the correlation between the DR vector and every expression vector.
    # r = file_operations.generate_correlation_map(dr_matrix, ge_matrix)
    # # r = file_operations.generate_correlation_map(sliced_dr, sliced_mat)
    # # Get the corresponding p-value for each Pearson correlation coefficient.
    # # Number of good indices?
    # p = [file_operations.compute_p_val(pcc, len(dr_matrix[0])) for pcc in r]

    # gene_list = gene_expr_df.index.values
    # # Build the dictionary mapping gene-drug pairs to their correlation p-value.
    # gene_drug_pearson_dct = {}
    # for i, gene in enumerate(gene_list):
    #     # Multiplying by length of gene list to do Bonferroni correction.
    #     pcc, p_value = r[i], p[i] * len(gene_list)
    #     if p_value < args.pearson_thresh:
    #         if args.sort_value == 'sortCorr':
    #             gene_drug_pearson_dct[(gene, drug)] = pcc
    #         else:
    #             gene_drug_pearson_dct[(gene, drug, pcc)] = p_value
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
    reverse = args.sort_value == 'sortCorr'
    sorted_pearson_dct = sorted(gene_drug_pearson_dct.items(),
        key=operator.itemgetter(1), reverse=reverse)

    gene_out = open('./results/%s_gene_drug_pearson_%s_%g.txt' % (args.input_file,
        args.sort_value, args.pearson_thresh), 'w')
    # gene_out.write('gene\tdrug\tpearson_correlation\tp_value\n')
    gene_out.write('gene\tdrug\tpearson_correlation\tp_value\n')
    # for (gene, drug, pcc), p_value in gene_drug_pearson_dct:
    #     gene_out.write('%s\t%s\t%g\t%g\n' % (gene, drug, pcc, p_value))
    if args.sort_value == 'sortCorr':
        for (gene, drug), pcc in sorted_pearson_dct:
            gene_out.write('%s\t%s\t%g\n' % (gene, drug, pcc))
    else:
        for (gene, drug, pcc), p_val in sorted_pearson_dct:
            gene_out.write('%s\t%s\t%g\t%g\n' % (gene, drug, pcc, p_val))
    gene_out.close()

def write_num_dcg(gene_drug_pearson_dct):
    '''
    For each drug, write the number of highly correlated genes.
    '''
    dcg_count_dct = {}
    for key in gene_drug_pearson_dct:
        gene, drug, pcc = key
        # Initialize a drug not in the dictionary.
        if drug not in dcg_count_dct:
            dcg_count_dct[drug] = 0
        dcg_count_dct[drug] += 1
    sorted_dcg_count_dct = sorted(dcg_count_dct.items(), key=operator.itemgetter(1),
        reverse=True)

    # Write out the file. # TODO: Change output name for each metohd. Same with regular thinig.
    out = open('./results/%s_drug_dcg_counts_%g.txt' % (args.input_file,
        args.pearson_thresh), 'w')
    for drug, count in sorted_dcg_count_dct:
        out.write('%s\t%d\n' % (drug, count))
    out.close()

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sort_value', help='Value by which to sort.',
        required=True, choices=['sortCorr', 'sortP'])
    parser.add_argument('-p', '--pearson_thresh', required=True, type=float,
        help='Pearson threshold that determines a significant correlation.')
    parser.add_argument('-i', '--input_file', required=True, type=str,
        choices=['demeter', 'rsa', 'ataris', 'ge'], help='Type of file to use for correlations.')
    args = parser.parse_args()

def main():
    parse_args()
    # if len(sys.argv) != 4:
    #     print 'Usage: %s top_k sortCorr/sortP pearson_thresh' % sys.argv[0]
    #     exit(1)
    # global top_k, sort_value, pearson_thresh
    # top_k, sort_value = int(sys.argv[1]), sys.argv[2]
    # pearson_thresh = float(sys.argv[3])
    # assert sort_value in ['sortCorr', 'sortP']

    if args.input_file == 'ge':
        fname = 'gene_expression_hgnc.tsv'
    elif args.input_file == 'demeter':
        fname = 'DEMETER_data.txt'
    elif args.input_file == 'ataris':
        fname = 'DRIVE_ATARiS_data.txt'
    elif args.input_file == 'rsa':
        fname = 'DRIVE_RSA_data.txt'

    # Drug response and gene expression information.
    global dr_df, gene_df
    dr_df = file_operations.get_dr_df()
    gene_df = file_operations.get_gene_df(fname)

    # Gene expression doesn't have any null values.
    if args.input_file == 'ge':
        assert not gene_df.isnull().values.any()
    else:
        # Rename columns for the gene essentiality samples.
        dr_df, gene_df = rename_cols(dr_df, gene_df)

    # Ensure the samples are aligned.
    assert list(dr_df) == list(gene_df)

    # NCI pathway information.
    # path_to_gene_dct, nci_genes = file_operations.get_path_to_gene_dct()
    # gene_universe = len(nci_genes.union(gene_list))

    # Number of drug-pathway pairs with Fisher p-value below FISHER_P_THRESH.
    # num_low_p = 0

    # Dictionaries to write out to file.
    gene_drug_pearson_dct, drug_path_fisher_dct = {}, {}

    pool = Pool(processes=20)
    dcts = pool.map(get_gene_drug_pearson_dct, dr_df.index.values)
    pool.close()
    pool.join()
    for dct in dcts:
        gene_drug_pearson_dct.update(dct)
    # # for drug in drug_to_dr_dct:
    # for drug in dr_df.index.values:
    #     # dr_vector = drug_to_dr_dct[drug]
    #     dr_vector = dr_df.loc[[drug]]

    #     # Get the all gene correlations to the current drug.
    #     # curr_corr_dct = get_gene_drug_pearson_dct(drug, dr_vector,
    #     #     gene_expr_mat, gene_list)
    #     curr_corr_dct = get_gene_drug_pearson_dct(dr_vector, gene_df)

    #     # Sort correlated genes by increasing p-value, then take top K.
    #     # TODO: sorting by PCC, rather than p-value.
    #     # top_genes = sorted(curr_corr_dct.items(), key=operator.itemgetter(1))
    #     # top_genes = sorted(curr_corr_dct.items(), key=operator.itemgetter(1),
    #     #     reverse=True)
    #     # top_genes = top_genes[:top_k]
    #     # # Keep only the gene names.
    #     # # top_genes = set([gene for (gene, drug, pcc), p_value in top_genes])
    #     # top_genes = set([gene for (gene, drug), pcc in top_genes])

    #     # Compute Fisher's test between pathway genes and current drug's most
    #     # correlated genes.
    #     # curr_fish_dct, curr_num_low_p = get_drug_path_fisher_dct(drug,
    #     #     gene_universe, top_genes, path_to_gene_dct)
    #     # Update dictionary for all drugs.
    #     # drug_path_fisher_dct.update(curr_fish_dct)
    #     # num_low_p += curr_num_low_p

    #     # Lastly, update the master dictionary for all drugs.
    #     gene_drug_pearson_dct.update(curr_corr_dct)

    # Write dictionaries out to file.
    # write_drug_path_correlations(drug_path_fisher_dct, num_low_p)
    write_gene_drug_pearson_dct(gene_drug_pearson_dct)

    # For each drug, get the number of correlated genes.
    write_num_dcg(gene_drug_pearson_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))