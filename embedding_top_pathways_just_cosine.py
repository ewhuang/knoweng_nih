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
    embedding_and_expression_genes, nci_path_dct):
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
        if (entity not in embedding_and_expression_genes) and (entity
            not in nci_path_dct):
            continue
        # Convert lines to floats.
        assert entity not in entity_vector_dct
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

def compute_drug_path_score(drug, gene_p_val_dct, entity_vector_dct,
    nci_path_dct):
    '''
    Given a dictionary where keys are genes and values are p-values of the
    correlation between the key and the input drug, we return a dictionary.
    Each key is a (drug, pathway) pair.
    Each value is the NetPath score, given by the sum of the cosine between
    each path and each gene of the highly correlated gene set, multiplied by
    that gene's p-value for the input drug.
    '''
    drug_path_score_dct = {}

    nci_pathways = nci_path_dct.keys()
    expression_genes = gene_p_val_dct.keys()

    pathway_vector_matrix = create_pathway_vector_matrix(nci_pathways,
        entity_vector_dct)
    gene_vector_matrix = create_gene_vector_matrix(expression_genes,
        entity_vector_dct)
    cosine_matrix = (pathway_vector_matrix * gene_vector_matrix.T
        ) / linalg.norm(pathway_vector_matrix, axis=1
        )[:, np.newaxis] / linalg.norm(gene_vector_matrix, axis=1)

    for path_i, pathway in enumerate(nci_pathways):
        # Initialize the score.
        drug_path_score = 0.0
        for gene_i, gene in enumerate(expression_genes):
            # Correlation between gene and drug response.
            gene_drug_corr = gene_p_val_dct[gene]
            # cos = cosine_dct[gene]
            cos = cosine_matrix[path_i, gene_i]
            # drug_path_score += abs(cos * gene_drug_corr)
            drug_path_score += cos
        drug_path_score_dct[(drug, pathway)] = drug_path_score
    return drug_path_score_dct

def write_top_pathway_file(network, extension, top_k, drug_path_score_dct,
    subfolder, out_fname):
    out = open(subfolder + out_fname, 'w')
    out.write('filler\ndrug\tpath\tscore\n')
    for (drug, pathway), score in drug_path_score_dct:
        out.write('%s\t%s\t%f\n' % (drug, pathway, score))
    out.close()

# Writes out to file the top drug-pathway pairs with score equal to the inverse
# of their overall ranking. 
def write_inverse_rankings(subfolder, filename):   
    brd_drug_to_name_dct = file_operations.get_brd_drug_to_name_dct()

    # print 'Taking out top 20 global pathways'
    # top_global_pathways = file_operations.get_top_global_pathways()[:20]

    drug_path_dct = {}
    f = open(subfolder + filename, 'r')
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        line = line.strip().split('\t')
        assert len(line) == 3
        drug, path, score = line
        # if path in top_global_pathways:
        #     continue
        drug_path_dct[(drug, path)] = float(score)
    f.close()

    # Sort ppi dictionary by value.
    ranked_drug_path_dct = {}
    drug_path_dct = sorted(drug_path_dct.items(), key=operator.itemgetter(1),
        reverse=True)
    for i, ((drug, path), score) in enumerate(drug_path_dct):
        inverse_ranking = (i + 1) / float(len(drug_path_dct))
        ranked_drug_path_dct[(drug, path)] = inverse_ranking

    ranked_drug_path_dct = sorted(ranked_drug_path_dct.items(),
        key=operator.itemgetter(1), reverse=True)

    out = open(subfolder + 'inverse_rankings/' + 'inverse_' + filename, 'w')
    out.write('drug\tpath\tinverse_rank\n')
    for (drug, path), score in ranked_drug_path_dct:
        out.write('%s\t%s\t%g\n' % (brd_drug_to_name_dct[drug], path, score))
    out.close()

def find_top_pathways(network, top_k):
    # Extract the NCI pathway data.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
    # Find genes and pathways that appear in embedding.
    embedding_gene_pathway_lst = (
        file_operations.get_embedding_gene_pathway_lst())
    # Find the most significantly correlated genes for each drug.
    (embedding_and_expression_genes, drug_top_genes_dct
        ) = file_operations.get_corr_drug_top_genes(top_k,
        embedding_gene_pathway_lst)

    for dimension in ['50']:
        for suffix in ['U']:
            entity_vector_dct, extension = get_entity_vector(network, dimension,
                suffix, embedding_gene_pathway_lst,
                embedding_and_expression_genes, nci_path_dct)

            # Calculate the score for each drug-pathway pair.
            drug_path_score_dct = {}
            for drug in drug_top_genes_dct:
                gene_p_val_dct = drug_top_genes_dct[drug]
                drug_path_score_dct.update(compute_drug_path_score(drug,
                    gene_p_val_dct, entity_vector_dct, nci_path_dct))

            # Sort the score dictionary, and write to file.
            drug_path_score_dct = sorted(drug_path_score_dct.items(),
                key=operator.itemgetter(1), reverse=True)
            
            subfolder = './results/random_embedding/'

            out_fname = '%s_top_pathways_%s_top_%d_just_cosine.txt' % (network, extension,
                top_k)
            write_top_pathway_file(network, extension, top_k,
                drug_path_score_dct, subfolder, out_fname)

            write_inverse_rankings(subfolder, out_fname)

def main():    
    find_top_pathways('ppi', 250)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))