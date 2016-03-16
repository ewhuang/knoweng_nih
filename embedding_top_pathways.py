### Author: Edward Huang

from collections import OrderedDict
import operator
import math
import sys

### This script opens the embedding data and performs drug/pathway analysis.
### We first translate the rows to the right genes and pathways, given by the
### file gene_pathway_id.txt. For each drug-pathway pair, we find a score.
### This score is given by the summation of the top 5 most important genes
### for a drug, taking the cosine similarity between the pathway and the
### gene vectors in the embedding file, multiplied by the correlation between
### the gene and drug response for the drug.

def cosine_similarity(v1, v2):
    "compute cosine similarity of v1 to v2: (v1 dot v2)/{||v1||*||v2||)"
    sumxx, sumxy, sumyy = 0, 0, 0
    for i in range(len(v1)):
        x = v1[i]; y = v2[i]
        sumxx += x * x
        sumyy += y * y
        sumxy += x * y
    return sumxy / math.sqrt(sumxx * sumyy)

def find_top_pathways(network, top_k):
    pathways = set([])
    f = open('./data/nci_pathway_hgnc.txt', 'r')
    for line in f:
        pathway, gene = line.strip().split('\t')
        pathways.add(pathway)
    f.close()

    # Find genes and pathways that appear in embedding.
    gene_pathway_lst = []
    f = open('./data/embedding/gene_pathway_id.txt', 'r')
    for line in f:
        # entity is either a gene or a pathway.
        gene_pathway_lst += [line.strip()]
    f.close()

    # Find the top genes for each drug.
    shared_genes = set([])
    drug_top_genes_dct = {}
    f = open('./results/top_genes_exp_hgnc.txt', 'r')
    for line in f:
        gene, drug, p_val = line.split()
        if gene not in gene_pathway_lst:
            continue
        # The gene appears in both expression and embedding.
        shared_genes.add(gene)
        p_val = float(p_val)
        if drug not in drug_top_genes_dct:
            drug_top_genes_dct[drug] = {gene: p_val}
            assert (len(drug_top_genes_dct[drug]) < top_k)
        elif len(drug_top_genes_dct[drug]) == top_k:
            # We have added enough genes for the drug.
            continue
        else:
            drug_top_genes_dct[drug][gene] = p_val
    f.close()

    # Check that each drug has exactly top_k number of genes.
    for drug in drug_top_genes_dct:
        assert len(drug_top_genes_dct[drug]) <= top_k

    # Loop through the files.
    if network == 'ppi':
        # We ran more dimensions for the ppi network.
        dimension_list = [50, 100, 500, 1000, 1500, 2000]
    else:
        dimension_list = [50, 100, 500]
    for dim in map(str, dimension_list):
        for suffix in ['U', 'US']:
            entity_vector_dct = OrderedDict({})

            extension = '%s_0.8.%s' % (dim, suffix)
            if network == 'ppi':
                filename = './data/embedding/ppi_6_net_%s' % extension
            else:
                filename = './data/embedding/%s.network_net_%s' % (network, 
                    extension)
            f = open(filename, 'r')
            for i, line in enumerate(f):
                # There should be the same number of rows in the embedding files
                # as the mapping file.
                entity = gene_pathway_lst[i]
                # Skip genes and pathways we have not previously used.
                if entity not in shared_genes and entity not in pathways:
                    continue
                # Convert lines to floats.
                entity_vector_dct[entity] = map(float, line.split())
            f.close()

            # Calculate the score for each drug-pathway pair.
            out = open('./results/embedding/%s_top_pathways_%s_top_%d.txt' 
                % (network, extension, top_k), 'w')
            score_dct = {}
            out.write('filler\n')
            out.write('drug\tpath\tscore\n')
            for drug in drug_top_genes_dct:
                gene_p_val_dct = drug_top_genes_dct[drug]

                # Get the cosine similarity between each pathway vector and gene
                # vector.
                for pathway in pathways:
                    # Initialize the score.
                    drug_path_score = 0.0
                    pathway_vec = entity_vector_dct[pathway]

                    for gene in gene_p_val_dct:
                        gene_vec = entity_vector_dct[gene]
                        cos = cosine_similarity(pathway_vec, gene_vec)
                        # Also get the correlation between gene and drug
                        # response.
                        gene_drug_corr = gene_p_val_dct[gene]
                        drug_path_score += cos * gene_drug_corr
                    score_dct[(drug, pathway)] = drug_path_score
            # Sort the score dictionary, and write out to file.
            score_dct = sorted(score_dct.items(), key=operator.itemgetter(1), 
                reverse=True)
            for (drug, pathway), score in score_dct:
                out.write('%s\t%s\t%f\n' % (drug, pathway, score))
            out.close()

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print "Usage: " + sys.argv[0] + " ppi/genetic/literome/sequence top_k"
        exit(1)
    network = sys.argv[1]
    top_k = int(sys.argv[2])
    assert network in ['ppi', 'genetic', 'literome', 'sequence']
    find_top_pathways(network, top_k)