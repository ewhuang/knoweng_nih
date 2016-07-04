### Author: Edward Huang

from collections import OrderedDict
import math
import random

### This file contains functions that parse the data files and return the 
### data objects that we work with in our scripts.

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

    nci_gene_set = set([gene for pathway in nci_path_dct.values(
        ) for gene in pathway])

    return nci_path_dct, nci_gene_set

def get_drug_response_dct():
    '''
    Returns a dictionary.
    Each key is a drug, and each value is a list of drug responses for that drug
    across all patients.
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

def get_gene_expression_dct():
    '''
    Returns a dictionary.
    Each key is a gene, each value is a list of gene expression values across
    all patients.
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

# Get the top pathways for LINCS. Keys are (drug, pathway) tuples, and values
# are floats of the p-values.
def min_p_exp(p_val_lst):
    return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))
def get_lincs_drug_path_dct():
    # print 'Extracting level 4 LINCS top pathways...'
    f = open('./results/top_pathways_lincs_diff_normalize_DMSO.txt', 'r')
    lincs_drug_path_dct = {}
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
        lincs_drug_path_dct[(drug, path)] = min_p_exp(p_val_lst)
    return lincs_drug_path_dct

def get_lincs_drug_path_dct_kw():
    # print 'Extracting level 4 LINCS top pathways...'
    f = open('./results/top_pathways_lincs_diff_normalize_DMSO_kw.txt', 'r')
    lincs_drug_path_dct = {}
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
        lincs_drug_path_dct[(drug, path)] = min_p_exp(p_val_lst)
    return lincs_drug_path_dct

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

# Dictionary where keys are BRD drug id's, and values are the common names.
def get_brd_drug_to_name_dct():
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

def get_corr_drug_top_genes(top_k, embedding_gene_pathway_lst):
    '''
    Returns a dictionary finding top correlated genes for each drug.
    Each key is a drug name.
    Each value is another dictionary. Each key of sub-dictionary is a gene name.
    Each value is correlation coefficient.
    '''
    embedding_and_expression_genes = set([])
    drug_top_genes_dct = {}
    f = open('./results/top_genes_correlation_hgnc.txt', 'r')
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

    # Check that each drug has exactly top_k number of genes.
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