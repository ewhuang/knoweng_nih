### Author: Edward Huang

from collections import OrderedDict
import math

### This file contains functions that parse the data files and return the 
### data objects that we work with in our scripts.

# Returns the NCI pathways in the form of a dictionary. The keys are the names
# of the pathways, and the values are lists of the genes corresponding to the
# pathway keys. Also returns a set of all genes that appear in the pathways
# as the second element of the returned tuple.
def get_nci_path_dct():
    path_file = open('./data/nci_pathway_hgnc.txt', 'r')
    nci_path_dct = {}
    # Set of all genes that appear in all NCI pathways.
    nci_genes = set([])
    for line in path_file:
        line = line.strip().split('\t')
        assert len(line) == 2
        path_name, path_gene = line
        nci_genes.add(path_gene)
        if path_name in nci_path_dct:
            nci_path_dct[path_name] += [path_gene]
        else:
            nci_path_dct[path_name] = [path_gene]
    path_file.close()
    return nci_path_dct, nci_genes

# Get the drug responses from the spreadsheet file. Keys are drugs, values are
# lists of drug responses across all patients.
def get_drug_resp_dct():
    drug_resp_dct = OrderedDict({})
    resp_file = open('./data/auc_hgnc.tsv', 'r')
    for i, line in enumerate(resp_file):
        # Header line contains patient ID's.
        if i == 0:
            drug_response_patients = line.split()[1:]
            continue
        # Each row is one drug's performance on each patient.
        line = line.split()
        drug, resp_line = line[0], line[1:]
        assert len(resp_line) == len(drug_response_patients)
        assert 'BRD-' in drug
        # Convert 'NA' strings to None values.
        resp_line = [None if resp == 'NA' else float(resp) for resp in resp_line]
        drug_resp_dct[drug] = resp_line
    resp_file.close()
    return drug_resp_dct

# Get expression information. Keys are genes, and values are lists of floats,
# which are the expression values across patients.
def get_exp_dct():
    exp_dct = OrderedDict({})
    exp_file = open('./data/gene_expression_hgnc.tsv', 'r')
    for i, line in enumerate(exp_file):
        # Header row contains patient ID's.
        if i == 0:
            expression_patients = line.split()[1:]
            continue
        line = line.split()
        gene, exp_line = line[0], line[1:]
        assert len(exp_line) == len(expression_patients)
        exp_dct[gene] = map(float, exp_line)
    exp_file.close()
    return exp_dct

# Get the top pathways for LINCS. Keys are (drug, pathway) tuples, and values
# are floats of the p-values.
def min_p_exp(p_val_lst):
    return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))
def get_lincs_drug_path_dct():
    # print 'Extracting level 4 LINCS top pathways...'
    f = open('./results/top_pathways_lincs_Aft_3_hgnc.txt', 'r')
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
    index_rankings = []
    f = open('./data/genetic.networktop.pathway', 'r')
    for line in f:
        index_rankings += [index_to_pathway_dct[int(line.strip()) - 1]]
    f.close()

    return index_rankings