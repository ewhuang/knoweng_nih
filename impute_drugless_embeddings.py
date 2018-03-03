### Author: Edward Huang

import file_operations
import numpy as np
import os

### This script takes the embedding files that do not contain pathway embeddings,
### and adds pathway embeddings by averaging, for each pathway, its member genes'
### embedding vectors.

def read_embedding_file(fname, emb_node_lst):
    '''
    Returns a dictionary.
    Key: gene -> str
    Value: embedding vector -> list(float)
    '''
    node_embedding_dct = {}
    f = open(fname, 'r')
    for i, line in enumerate(f):
        # Each row in the file maps to the row in gene_pathway_id.txt.
        node = emb_node_lst[i]
        assert ' ' not in node # Pathways have spaces.
        node_embedding_dct[node] = map(float, line.split())
    f.close()
    assert len(node_embedding_dct) == 18362
    return node_embedding_dct

def get_pathway_embeddings(node_embedding_dct, path_to_gene_dct):
    '''
    Get the pathways' embeddings by averaging their member genes' embedding
    vectors.
    '''
    pathway_embedding_dct = {}
    for pathway in path_to_gene_dct:
        gene_set = path_to_gene_dct[pathway]
        # Get the embedding vectors for the gene set.
        embedding_vector_lst = []
        for gene in gene_set:
            if gene in node_embedding_dct:
                embedding_vector_lst += [node_embedding_dct[gene]]
        # Average the vectors for the gene set.
        pathway_embedding_dct[pathway] = np.mean(embedding_vector_lst, axis=0)
    return pathway_embedding_dct

def write_embedding_file(fname, pathway_embedding_dct, node_embedding_dct, emb_node_lst):
    '''
    Given the embedding dictionary for the pathways, write it out to file by
    appending to the original gene vectors.
    '''
    # Rename the new output file.
    new_fname = fname.replace('no_drug', '')
    out = open(new_fname, 'w')
    for node in emb_node_lst:
        if node in node_embedding_dct:
            embedding_vector = node_embedding_dct[node]
        else:
            embedding_vector = pathway_embedding_dct[node]
        out.write('\t'.join(map(str, embedding_vector)) + '\n')
    out.close()

def main():
    emb_node_lst = file_operations.get_emb_node_lst()
    path_to_gene_dct, nci_gene_set = file_operations.get_path_to_gene_dct()

    for subdir, dirs, files in os.walk('./data/embedding_new'):
        for fname in files:
            if 'no_drug' not in fname:
                continue
            fname = '%s/%s' % (subdir, fname)

            node_embedding_dct = read_embedding_file(fname, emb_node_lst)
            
            pathway_embedding_dct = get_pathway_embeddings(node_embedding_dct,
                path_to_gene_dct)

            write_embedding_file(fname, pathway_embedding_dct,
                node_embedding_dct, emb_node_lst)

if __name__ == '__main__':
    main()