### Author: Edward Huang

# from collections import OrderedDict
import file_operations
import numpy as np
import operator
import os
from scipy import linalg
import sys
import time

### For each drug, we take the top K most correlated genes, as computed by the
### drug_pathway_fisher_correlation.py script. Then, for each drug-pathway pair,
### we compute a NetPath score by multiplying each gene's correlation to that 
### drug with the gene-pathway cosine similarity scores.
### Run time: 6 minutes per embedding file.

def get_embedding_dct(dimension, suffix, fraction, emb_node_lst, all_genes,
    path_to_gene_dct):
    '''
    Returns a dictionary.
    Key: node (either a gene or a pathway) -> str
    Value: key's corresponding embedding vector -> list(float)
    '''
    node_embedding_dct = {}

    if isPpi:
        data_folder = './data/embedding' # With PPI.
        filename = '%s/ppi_6_net_%s_%s.%s' % (data_folder, dimension, fraction,
            suffix)
    else:
        data_folder = './embedding' # string_experimental
        filename = '%s/string_experimental_%s_%s.%s' % (data_folder, dimension,
        fraction, suffix)

    f = open(filename, 'r')
    for i, line in enumerate(f):
        # Each row in the file maps to the row in gene_pathway_id.txt.
        node = emb_node_lst[i]
        # Skip genes and pathways that aren't in expression or NCI pathways.
        if (node not in all_genes) and (node not in path_to_gene_dct):
            continue
        # Convert lines to floats.
        assert node not in node_embedding_dct
        node_embedding_dct[node] = map(float, line.split())
    f.close()
    return node_embedding_dct

def create_embedding_matrix(node_list, node_embedding_dct):
    '''
    Given a list of nodes (either expression genes or pathways), output a
    matrix where each row is a node's embedding vector.
    '''
    embedding_matrix = []
    for node in node_list:
        embedding_matrix += [node_embedding_dct[node]]
    return np.matrix(embedding_matrix)

def compute_drug_path_score(drug, gene_p_val_dct, node_embedding_dct,
    path_to_gene_dct):
    '''
    Given a dictionary where keys are genes and values are p-values of the
    correlation between the key and the input drug, we return a dictionary.
    Key: (drug, pathway) pair -> (str, str)
    Value: NetPath score, the sum of the cosine between each path and each gene
        of the highly correlated gene set, multiplied by that gene's p-value for
        the input drug -> float
    '''
    drug_path_score_dct = {}

    # Get embedding matrices for pathways and genes.
    nci_pathways = path_to_gene_dct.keys()
    expr_genes = gene_p_val_dct.keys()
    path_emb_matrix = create_embedding_matrix(nci_pathways, node_embedding_dct)
    gene_emb_matrix = create_embedding_matrix(expr_genes, node_embedding_dct)

    # Compute the cosine matrix between pathways and genes.
    cosine_matrix = (path_emb_matrix * gene_emb_matrix.T) / linalg.norm(
        path_emb_matrix, axis=1)[:, np.newaxis] / linalg.norm(gene_emb_matrix,
        axis=1)

    for path_i, pathway in enumerate(nci_pathways):
        # Initialize the score.
        drug_path_score = 0.0
        for gene_i, gene in enumerate(expr_genes):
            # Correlation between gene and drug response.
            gene_drug_corr = gene_p_val_dct[gene]
            cos = cosine_matrix[path_i, gene_i]
            if cos_abs:
                drug_path_score += abs(cos * gene_drug_corr)
            else:
                drug_path_score += cos * gene_drug_corr
        drug_path_score_dct[(drug, pathway)] = drug_path_score
    return drug_path_score_dct

def write_top_pathway_file(drug_path_score_dct, results_folder, out_fname):
    out = open(results_folder + out_fname, 'w')
    # Description of method.
    out.write('cosine * gene_drug_corr, unnormalized')
    out.write('\ndrug\tpath\tscore\n')
    for (drug, pathway), score in drug_path_score_dct:
        out.write('%s\t%s\t%f\n' % (drug, pathway, score))
    out.close()

# Writes out to file the top drug-pathway pairs with score equal to the inverse
# of their overall ranking. 
def write_inverse_rankings(results_folder, filename):   
    brd_drug_to_name_dct = file_operations.get_brd_drug_to_name_dct()

    drug_path_dct = {}
    f = open(results_folder + filename, 'r')
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        line = line.strip().split('\t')
        assert len(line) == 3
        drug, path, score = line
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

    inverse_folder = '%sinverse_rankings/' % results_folder
    if not os.path.exists(inverse_folder):
        os.makedirs(inverse_folder)

    out = open(inverse_folder + 'inverse_' + filename, 'w')
    out.write('drug\tpath\tinverse_rank\n')
    for (drug, path), score in ranked_drug_path_dct:
        out.write('%s\t%s\t%g\n' % (brd_drug_to_name_dct[drug], path, score))
    out.close()

def compute_drug_pathway_scores():
    # Extract the NCI pathway data.
    path_to_gene_dct, nci_genes = file_operations.get_path_to_gene_dct()
    # Find genes and pathways that appear in embedding.
    emb_node_lst = file_operations.get_emb_node_lst()
    # Find the most significantly correlated genes for each drug.
    all_genes, drug_corr_genes_dct = file_operations.get_drug_corr_genes_dct(
        top_k, emb_node_lst, sort_value, pearson_thresh)

    dimension_list = map(str, [50, 100, 500, 1000])
    fraction_list, suffix_list = map(str, [0.3, 0.5, 0.8]), ['U', 'US']
    if isPpi:
        dimension_list, fraction_list, suffix_list = ['50'], ['0.8'], ['U']

    for dimension in dimension_list:
        for fraction in fraction_list:
            for suffix in suffix_list:
                # Get the embedding vectors for each node (gene or pathway).
                node_embedding_dct = get_embedding_dct(dimension, suffix,
                    fraction, emb_node_lst, all_genes, path_to_gene_dct)

                # Calculate the score for each drug-pathway pair.
                drug_path_score_dct = {}
                for drug in drug_corr_genes_dct:
                    gene_p_val_dct = drug_corr_genes_dct[drug]
                    drug_path_score_dct.update(compute_drug_path_score(drug,
                        gene_p_val_dct, node_embedding_dct, path_to_gene_dct))

                # Sort the score dictionary, and write to file.
                drug_path_score_dct = sorted(drug_path_score_dct.items(),
                    key=operator.itemgetter(1), reverse=True)
                
                results_folder = './results/embedding/'
                if not os.path.exists(results_folder):
                    os.makedirs(results_folder)

                input_suffix = '%s_%g' % (sort_value, pearson_thresh)
                if cos_abs:
                    input_suffix += '_cosAbs'
                else:
                    input_suffix += '_cosNoAbs'
                if isPpi:
                    out_fname = 'top_pathways_ppi_%s_%s_%s_embed_%d_%s.txt' % (
                        dimension, fraction, suffix, top_k, input_suffix)
                else:
                    out_fname = 'top_pathways_%s_%s_%s_embed_%d_%s.txt' % (
                        dimension, fraction, suffix, top_k, input_suffix)
                write_top_pathway_file(drug_path_score_dct, results_folder,
                    out_fname)

                write_inverse_rankings(results_folder, out_fname)

def main():
    if (len(sys.argv) not in [5, 6]):
        print ('Usage: %s top_k sortCorr/sortP pearson_thresh cos_abs '
            'ppi<optional>' % (sys.argv[0]))
        exit(1)
    global isPpi, top_k, sort_value, pearson_thresh, cos_abs
    isPpi = len(sys.argv) == 6
    top_k = int(sys.argv[1])
    sort_value = sys.argv[2]
    assert sort_value in ['sortCorr', 'sortP']
    pearson_thresh = float(sys.argv[3])
    cos_abs = sys.argv[4] == 'cos_abs'
    
    compute_drug_pathway_scores()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))