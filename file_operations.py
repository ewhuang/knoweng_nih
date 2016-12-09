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
# correlation_top_pathways_fisher.py
# correlation_top_pathways_kw.py
# embedding_top_pathways.py
def get_nci_path_dct():
    '''
    Returns a (dictionary, set) pair.
    Dictionary: the NCI pathways in the form of a dictionary. Each key is the
    name of a pathway. Each value is a list of genes associated with the pathway
    key.
    Set: All genes that appear in all pathways.
    '''
    # Pair to return.
    nci_path_dct = {}

    path_file = open('./data/nci_pathway_hgnc.txt', 'r')
    for line in path_file:
        pathway_name, gene = line.strip().split('\t')
        
        # Update pathway in the dictionary.
        if pathway_name in nci_path_dct:
            nci_path_dct[pathway_name] += [gene]
        else:
            nci_path_dct[pathway_name] = [gene]
    path_file.close()

    # Flatten the values in the dictionary to get the set of all NCI genes.
    nci_gene_set = set([gene for pathway in nci_path_dct.values(
        ) for gene in pathway])

    return nci_path_dct, nci_gene_set

# correlation_top_pathways_fisher.py
# correlation_top_pathways_kw.py
def get_drug_response_dct():
    '''
    Returns a dictionary mapping drugs to drug response values.
    Key: drug -> str
    Value: list of drug responses for the drug key -> list(float)
    '''
    drug_response_dct = OrderedDict({})
    resp_file = open('./data/auc_hgnc.tsv', 'r')
    for i, line in enumerate(resp_file):
        # Header line contains patient ID's.
        if i == 0:
            patient_id_list = line.split()[1:]
            continue
        # Each row is one drug's performance on each patient.
        line = line.split()
        drug, resp_line = line[0], line[1:]
        assert len(resp_line) == len(patient_id_list)
        assert 'BRD-' in drug
        # Convert 'NA' strings to None values.
        resp_line = [None if resp == 'NA' else float(resp
            ) for resp in resp_line]
        assert drug not in drug_response_dct
        drug_response_dct[drug] = resp_line
    resp_file.close()
    return drug_response_dct

# correlation_top_pathways_fisher.py
# correlation_top_pathways_kw.py
def get_gene_expression_dct():
    '''
    Returns a dictionary mapping genes to gene expression values.
    Key: gene -> str
    Value: list of gene expression values -> list(float)
    '''
    gene_expression_dct = OrderedDict({})
    exp_file = open('./data/gene_expression_hgnc.tsv', 'r')
    for i, line in enumerate(exp_file):
        # Header row contains patient ID's.
        if i == 0:
            # The first item is 'gid'. We skip it.
            patient_id_list = line.split()[1:]
            continue
        line = line.split()
        gene, expression_value_list = line[0], line[1:]
        assert len(expression_value_list) == len(patient_id_list)

        gene_expression_dct[gene] = map(float, expression_value_list)
    exp_file.close()
    return gene_expression_dct

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
def get_embedding_gene_pathway_lst():
    '''
    Returns a list of 'entities' that appear in the embedding data. Entities
    include both genes and pathways.
    '''
    embedding_gene_pathway_lst = []
    f = open('./data/embedding/gene_pathway_id.txt', 'r')
    for line in f:
        # entity is either a gene or a pathway.
        embedding_gene_pathway_lst += [line.strip()]
    f.close()
    return embedding_gene_pathway_lst

# embedding_top_pathways.py
def get_corr_drug_top_genes(top_k, embedding_gene_pathway_lst):
    '''
    Returns a dictionary finding top correlated genes for each drug.
    Each key is a drug name.
    Each value is another dictionary. Each key of sub-dictionary is a gene name.
    Each value is correlation coefficient.
    '''
    embedding_and_expression_genes = set([])
    drug_top_genes_dct = {}
    f = open('./results/correlation_top_genes_kw.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        gene, drug, correlation, p_val = line.split()
        if gene not in embedding_gene_pathway_lst:
            continue
        # The gene appears in both expression and embedding.
        embedding_and_expression_genes.add(gene)
        correlation = float(correlation)
        if drug not in drug_top_genes_dct:
            drug_top_genes_dct[drug] = {gene: correlation}
            assert (len(drug_top_genes_dct[drug]) < top_k)
        elif len(drug_top_genes_dct[drug]) == top_k:
            # We have added enough genes for the drug.
            continue
        else:
            drug_top_genes_dct[drug][gene] = correlation
    f.close()

    # Check that each drug has at most top_k genes.
    for drug in drug_top_genes_dct:
        assert len(drug_top_genes_dct[drug]) <= top_k

    return embedding_and_expression_genes, drug_top_genes_dct

def get_corr_drug_random_genes(top_k, embedding_gene_pathway_lst):
    '''
    Returns a dictionary where each key is a drug, and each value is another
    dictionary d_2, where each key is a gene, and each value is a p-value
    corresponding to the correlation between the drug and the gene. We pick
    top_k random genes for each drug.
    '''
    gene_expression_dct = get_gene_expression_dct()
    embedding_and_expression_genes = list(set(embedding_gene_pathway_lst
        ).intersection(gene_expression_dct.keys()))

    random.shuffle(embedding_and_expression_genes)

    return embedding_and_expression_genes[:top_k]

# correlation_top_pathways_fisher.py
# correlation_top_pathways_kw.py
def generate_correlation_map(x, y):
    '''
    Correlate each n with each m, where x is an N x T matrix, and y is an M x T
    matrix.
    '''
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                  mu_y[np.newaxis, :])
    return (cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :]))[0]

# correlation_top_pathways_fisher.py
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