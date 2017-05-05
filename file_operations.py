### Author: Edward Huang

from collections import OrderedDict
import math
import numpy as np
import random
from scipy.stats import betai

### This file contains functions that parse the data files and return the 
### data objects that we work with in our scripts.

# lincs_top_pathways.py
def get_lincs_genes():
    '''
    Each line of all_map.txt is %s\t%s\n. The second item is the LINCS gene.
    Returns a list of genes -> list(str)
    '''
    lincs_genes = []
    f = open('./data/all_map.txt', 'r')
    for line in f:
        lincs_genes += [line.split()[1]]
    f.close()
    return lincs_genes

# lincs_top_pathways.py
def get_drugs_and_gene_matrix():
    '''
    Reads the LINCS data, and returns a (list, list) tuple.
    list_0: list of drug names -> list(str)
    list_1: 2D list of pre-normalization z-scores -> list(list(str))
    Each row of list_1 corresponds to a gene in list_0.
    '''
    f = open('./data/new_lvl4_Stuart_combinedPvalue_diff_normalize_DMSO.txt',
        'r')
    drugs, gene_matrix = [], []
    for line in f:
        line = line.split()
        drug, raw_z_scores = line[0], line[1:]
        drugs += [drug]
        gene_matrix += [raw_z_scores]
    f.close()
    return drugs, gene_matrix

# lincs_top_pathways.py
# drug_pathway_fisher_correlation.py
# correlation_top_pathways_kw.py
# embedding_top_pathways.py
def get_path_to_gene_dct():
    '''
    Returns a (dictionary, set) pair.
    Dictionary: the NCI pathways in the form of a dictionary. Each key is the
    name of a pathway. Each value is a list of genes associated with the pathway
    key.
    Set: All genes that appear in all pathways.
    '''
    # Pair to return.
    path_to_gene_dct, nci_gene_set = {}, set([])

    path_file = open('./data/nci_pathway_hgnc.txt', 'r')
    for line in path_file:
        pathway_name, gene = line.strip().split('\t')        
        # Update pathway in the dictionary.
        if pathway_name not in path_to_gene_dct:
            path_to_gene_dct[pathway_name] = set([])
        path_to_gene_dct[pathway_name].add(gene)
        nci_gene_set.add(gene)
    path_file.close()

    return path_to_gene_dct, nci_gene_set

# drug_pathway_fisher_correlation.py
# correlation_top_pathways_kw.py
def get_drug_to_dr_dct():
    '''
    Returns a dictionary mapping drugs to drug response values.
    Key: drug -> str
    Value: list of drug responses for the drug key -> list(float)
    '''
    drug_to_dr_dct = OrderedDict({})
    resp_file = open('./data/auc_hgnc.tsv', 'r')
    for i, line in enumerate(resp_file):
        line = line.split()
        # Header line contains patient ID's.
        if i == 0:
            patient_id_list = line[1:]
            continue
        # Each row is one drug's performance on each patient.
        drug, resp_line = line[0], line[1:]
        # Skip drugs that are unavailable for all patients.
        if resp_line == ['NA'] * len(patient_id_list):
            continue
        assert len(resp_line) == len(patient_id_list)
        assert 'BRD-' in drug
        # Convert 'NA' strings to -666.
        resp_line = [-666 if val == 'NA' else float(val) for val in resp_line]
        assert drug not in drug_to_dr_dct
        drug_to_dr_dct[drug] = np.array(resp_line)
    resp_file.close()
    return drug_to_dr_dct

# drug_pathway_fisher_correlation.py
# correlation_top_pathways_kw.py
def get_gene_expr_mat():
    '''
    Returns a dictionary mapping genes to gene expression values.
    Key: gene -> str
    Value: list of gene expression values -> list(float)
    '''
    gene_expr_mat, gene_list = [], []
    exp_file = open('./data/gene_expression_hgnc.tsv', 'r')
    for i, line in enumerate(exp_file):
        line = line.split()
        # Header row contains patient ID's.
        if i == 0:
            # The first item is 'gid'. We skip it.
            patient_id_list = line[1:]
            continue
        gene, expression_value_list = line[0], map(float, line[1:])
        assert len(expression_value_list) == len(patient_id_list)

        # Take care of duplicate genes. Average existing expression rows.
        if gene in gene_list:
            # Only 'TTL' is the duplicate.
            assert gene == 'TTL'
            dup_idx = gene_list.index(gene)
            # Average old row with current one.
            # old_row = gene_expr_mat[dup_idx]
            # TODO:
            # mean_row = list(np.mean([expression_value_list, old_row], axis=0))
            # Update the existing row.
            # gene_expr_mat[dup_idx] = mean_row
            gene_expr_mat[dup_idx] = expression_value_list
        else:
            # Only add a row and gene if it's a new gene.
            gene_expr_mat += [expression_value_list]
            gene_list += [gene]
    exp_file.close()
    assert len(gene_expr_mat) == len(gene_list)
    return np.array(gene_expr_mat), gene_list

# def min_p_exp(p_val_lst):
#     return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))

# compare_methods_with_lincs.py
def get_lincs_drug_path_dct(lincs_z, lincs_max_num):
    '''
    Returns a dictionary.
    Key: (drug, pathway) pairs -> (str, str)
    Value: p-values of the Fisher's test for LINCS -> float
    '''
    subfolder = './results/lincs_top_pathway_files'    
    lincs_drug_path_dct = {}

    f = open('%s/top_pathways_lincs_z%s_max%s.txt' % (subfolder, lincs_z,
        lincs_max_num), 'r')
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        drug, cell_line, path, score = line.strip().split('\t')[:4]
        # Add the drug-pathway p-value to the dictionary.
        score = float(score)
        if (drug, path) in lincs_drug_path_dct:
            lincs_drug_path_dct[(drug, path)] += [score]
        else:
            lincs_drug_path_dct[(drug, path)] = [score]
    f.close()
    # Aggregate p-values by cell lines.
    for (drug, path) in lincs_drug_path_dct:
        p_val_lst = lincs_drug_path_dct[(drug, path)]
        # Change the function for other aggregation functions.
        # lincs_drug_path_dct[(drug, path)] = min_p_exp(p_val_lst)
        lincs_drug_path_dct[(drug, path)] = min(p_val_lst)
    return lincs_drug_path_dct

# def get_lincs_drug_path_dct_kw():
#     # print 'Extracting level 4 LINCS top pathways...'
#     f = open('./results/top_pathways_lincs_diff_normalize_DMSO_kw.txt', 'r')
#     lincs_drug_path_dct = {}
#     for i, line in enumerate(f):
#         # Skip header lines.
#         if i < 2:
#             continue
#         drug, cell_line, path, score = line.strip().split('\t')[:4]
#         # Add the drug-pathway p-value to the dictionary.
#         score = float(score)
#         if (drug, path) in lincs_drug_path_dct:
#             lincs_drug_path_dct[(drug, path)] += [score]
#         else:
#             lincs_drug_path_dct[(drug, path)] = [score]
#     f.close()
#     # Aggregate p-values by cell lines.
#     for (drug, path) in lincs_drug_path_dct:
#         p_val_lst = lincs_drug_path_dct[(drug, path)]
#         # Change the function for other aggregation functions.
#         lincs_drug_path_dct[(drug, path)] = min_p_exp(p_val_lst)
#     return lincs_drug_path_dct

# Get the superdrug's p-values for all the pathways.
def get_superdrug_pathway_p_values():
    superdrug_pathway_p_values = {}
    f = open('./results/superdrug_pathway_p_values.txt', 'r')
    for line in f:
        pathway, p_value = line.strip().split('\t')
        superdrug_pathway_p_values[pathway] = float(p_value)
    f.close()
    return superdrug_pathway_p_values

# Gets the top genetic global pathways from random walk rankings.
def get_top_global_pathways():
    # Get the index to dictionary mapping.
    index_to_pathway_dct = {}
    f = open('./data/nci_pathway_name.txt', 'r')
    for i, line in enumerate(f):
        index_to_pathway_dct[i] = line.strip()
    f.close()

    # Get the rankings of indices.
    top_global_pathways = []
    f = open('./data/genetic.networktop.pathway', 'r')
    for line in f:
        top_global_pathways += [index_to_pathway_dct[int(line.strip()) - 1]]
    f.close()

    return top_global_pathways

# embedding_top_pathways.py
def get_brd_drug_to_name_dct():
    '''
    Returns a dictionary mapping BRD drug ID's to their common English names.
    Key: BRD ID -> str
    Value: common drug name -> str
    '''
    brd_drug_to_name_dct = {}
    f = open('./data/brd_to_name.txt', 'r')
    for i, line in enumerate(f):
        # Skip header line.
        if i == 0 or line == '\n':
            continue
        line = line.strip().split('\t')
        assert len(line) == 3
        drug_name = line[0]
        drug_id = line[2]
        assert drug_id not in brd_drug_to_name_dct
        brd_drug_to_name_dct[drug_id] = drug_name
    f.close()
    return brd_drug_to_name_dct

# embedding_top_pathways.py
def get_emb_node_lst():
    '''
    Returns a list of nodes that appear in the embedding data. Entities include
    both genes and pathways.
    '''
    emb_node_lst = []
    f = open('./data/embedding/gene_pathway_id.txt', 'r')
    for line in f:
        # Node is either a gene or a pathway.
        emb_node_lst += [line.strip()]
    f.close()
    return emb_node_lst

# embedding_top_pathways.py
def get_drug_corr_genes_dct(top_k, emb_node_lst):
    '''
    Returns a dictionary finding top correlated genes for each drug.
    Key: drug -> str
    Value: a nested dictionary -> {}
        Key: gene -> str
        Value: Pearson correlation coefficient -> float
    Also returns set of genes that appear in both embedding and expression data.
    '''
    all_genes = set([])
    drug_corr_genes_dct = {}

    f = open('./results/gene_drug_pearson.txt', 'r')
    for i, line in enumerate(f):
        if i == 0: # Skip header.
            continue
        # TODO: no p-value info right now.
        # gene, drug, correlation, p_val = line.split()
        gene, drug, correlation = line.split()
        if gene not in emb_node_lst:
            continue
        # The gene appears in both expression and embedding.
        all_genes.add(gene)
        # Initialize the dictionary.
        if drug not in drug_corr_genes_dct:
            drug_corr_genes_dct[drug] = {}
        if len(drug_corr_genes_dct[drug]) == top_k:
            continue
        drug_corr_genes_dct[drug][gene] = float(correlation)
    f.close()

    return all_genes, drug_corr_genes_dct

def get_corr_drug_random_genes(top_k, emb_node_lst):
    '''
    Returns a dictionary where each key is a drug, and each value is another
    dictionary d_2, where each key is a gene, and each value is a p-value
    corresponding to the correlation between the drug and the gene. We pick
    top_k random genes for each drug.
    '''
    gene_expr_mat = get_gene_expr_mat()
    embedding_and_expression_genes = list(set(emb_node_lst
        ).intersection(gene_expr_mat.keys()))

    random.shuffle(embedding_and_expression_genes)

    return embedding_and_expression_genes[:top_k]

# drug_pathway_fisher_correlation.py
# correlation_top_pathways_kw.py
def generate_correlation_map(x, y):
    '''
    Correlate each n with each m, where x is an N x T matrix, and y is an M x T
    matrix.
    From: http://stackoverflow.com/questions/30143417/
    '''
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x, y.T) - n * np.dot(mu_x[:, np.newaxis], mu_y[np.newaxis, :])
    return (cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :]))[0]

# drug_pathway_fisher_correlation.py
# correlation_top_pathways_kw.py
def compute_p_val(r, n):
    '''
    Given a Pearson correlation coefficient and a sample size, compute the p-
    value corresponding to that coefficient.
    '''
    df = n - 2
    if abs(r) == 1.0:
        prob = 0.0
    else:
        t_squared = r**2 * (df / ((1.0 - r) * (1.0 + r)))
        prob = betai(0.5*df, 0.5, df/(df+t_squared))
    return prob