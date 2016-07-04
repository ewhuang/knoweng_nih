### Author: Edward Huang

from collections import OrderedDict
import file_operations
import math
import numpy as np
import operator
from scipy import spatial, linalg
import sys
import time

### This script opens the embedding data and performs drug/pathway analysis.
### We first translate the rows to the right genes and pathways, given by the
### file gene_pathway_id.txt. For each drug-pathway pair, we find a score.
### This score is given by the summation of the most correlated genes
### for a drug, taking the cosine similarity between the pathway and the
### gene vectors in the embedding file, multiplied by the correlation between
### the gene and drug response for the drug.

def create_dimension_list(network):
    '''
    Creates the list of dimensions for which embedding was run. PPI was run open
    twice as many dimensions.
    '''
    dimension_list = [50, 100, 500]
    if network == 'ppi':
        # We ran more dimensions for the ppi network.
        dimension_list += [1000, 1500, 2000]
    return map(str, dimension_list)

def get_entity_vector(network, dimension, suffix, embedding_gene_pathway_lst,
    random_genes, nci_path_dct):
    '''
    Returns a dictionary.
    Each key is an entity (either a gene or a pathway).
    Each value is the low-dimensional vector representation of the entity as
    computed by embedding.
    '''
    entity_vector_dct = OrderedDict({})

    # First, process the input network filename.
    subfolder = './data/embedding/'
    extension = '%s_0.8.%s' % (dimension, suffix)
    if network == 'ppi':
        filename = '%sppi_6_net_%s' % (subfolder, extension)
    else:
        filename = '%s%s.network_net_%s' % (subfolder, network, extension)

    f = open(filename, 'r')
    for i, line in enumerate(f):
        # Each row in the file maps to the row in gene_pathway_id.txt.
        entity = embedding_gene_pathway_lst[i]
        # Skip genes and pathways that don't appear in expression genes.
        if (entity not in random_genes) and (entity not in nci_path_dct):
            continue
        # Convert lines to floats.
        entity_vector_dct[entity] = map(float, line.split())
    f.close()
    return entity_vector_dct, extension

def create_pathway_vector_matrix(nci_pathways, entity_vector_dct):
    pathway_vector_matrix = []
    for pathway in nci_pathways:
        pathway_vec = entity_vector_dct[pathway]
        pathway_vector_matrix += [pathway_vec]
    return np.matrix(pathway_vector_matrix)

def create_gene_vector_matrix(expression_genes, entity_vector_dct):
    gene_vector_matrix = []
    for gene in expression_genes:
        gene_vec = entity_vector_dct[gene]
        gene_vector_matrix += [gene_vec]
    return np.matrix(gene_vector_matrix)

def compute_drug_path_score(random_genes, nci_pathways, cosine_matrix,
    path_score_dct):
    '''
    Given a dictionary where keys are genes and values are p-values of the
    correlation between the key and the input drug, we return a dictionary.
    Each key is a (drug, pathway) pair.
    Each value is the NetPath score, given by the sum of the cosine between
    each path and each gene of the highly correlated gene set, multiplied by
    that gene's p-value for the input drug.
    '''
    for path_i, pathway in enumerate(nci_pathways):
        # Initialize the score.
        drug_path_score = 0.0
        for gene_i, gene in enumerate(random_genes):
            # Correlation between gene and drug response.
            # cos = cosine_dct[gene]
            # drug_path_score += abs(cos * gene_drug_corr)
            drug_path_score += cosine_matrix[path_i, gene_i]
        path_score_dct[pathway] += [drug_path_score]

def write_top_pathway_file(network, extension, top_k, path_score_dct, subfolder,
    out_fname):
    out = open(subfolder + out_fname, 'w')
    out.write('filler\npath\tscores\n')
    for pathway in path_score_dct:
        out.write('%s\t%s\n' % (pathway, '\t'.join(map(str,
            path_score_dct[pathway]))))
    out.close()

def find_top_pathways(network, top_k):
    # Extract the NCI pathway data.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
    # Find genes and pathways that appear in embedding.
    embedding_gene_pathway_lst = (
        file_operations.get_embedding_gene_pathway_lst())

    # Initialize the path score dictionary with all pathways as keys and empty
    # lists as values.
    path_score_dct = {}
    for pathway in nci_path_dct:
        path_score_dct[pathway] = []

    for i in range(10000):
        # Find the most significantly correlated genes for each drug.
        random_genes = file_operations.get_corr_drug_random_genes(top_k,
            embedding_gene_pathway_lst)

        # dimension_list = create_dimension_list(network)
        entity_vector_dct, extension = get_entity_vector(network, '50', 'U',
            embedding_gene_pathway_lst, random_genes, nci_path_dct)

        # Compute cosine matrix.
        nci_pathways = nci_path_dct.keys()
        pathway_vector_matrix = create_pathway_vector_matrix(nci_pathways,
            entity_vector_dct)
        gene_vector_matrix = create_gene_vector_matrix(random_genes,
            entity_vector_dct)
        cosine_matrix = (pathway_vector_matrix * gene_vector_matrix.T
            ) / linalg.norm(pathway_vector_matrix, axis=1
            )[:, np.newaxis] / linalg.norm(gene_vector_matrix, axis=1)

        # Calculate the score for each drug-pathway pair.
        compute_drug_path_score(random_genes, nci_pathways, cosine_matrix,
            path_score_dct)
    
    subfolder = './results/random_embedding/'

    out_fname = 'random_%s_top_pathways_%s_top_%d.txt' % (network, extension, top_k)
    write_top_pathway_file(network, extension, top_k, path_score_dct,
        subfolder, out_fname)

def main():
    find_top_pathways('ppi', 250)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))