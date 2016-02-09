### Author: Edward Huang

from scipy import spatial
from collections import OrderedDict
import operator

### This script opens the embedding data and performs drug/pathway analysis.
### We first translate the rows to the right genes and pathways, given by the
### file gene_pathway_id.txt. For each drug-pathway pair, we find a score.
### This score is given by the summation of the top 5 most important genes
### for a drug, taking the cosine similarity between the pathway and the
### gene vectors in the embedding file, multiplied by the correlation between
### the gene and drug response for the drug.

TOP_GENES_PER_DRUG = 5

if __name__ == '__main__':
    pathways = set([])
    f = open('./data/nci_pathway_hgnc.txt', 'r')
    for line in f:
        pathway, gene = line.strip().split('\t')
        pathways.add(pathway)
        # assert (len(line) == 2)
    f.close()

    # Find the shared genes.
    exp_genes = set([])
    f = open('./results/top_genes_exp_hgnc.txt', 'r')
    for line in f:
        gene, drug, p_val = line.split()
        exp_genes.add(gene)
    f.close()

    embed_genes = set([])
    f = open('./data/embedding/gene_pathway_id.txt', 'r')
    for line in f:
        embed_genes.add(line.strip())
    f.close()

    shared_genes = exp_genes.intersection(embed_genes)

    # Find the top genes for each drug.
    drug_top_genes_dct = {}
    f = open('./results/top_genes_exp_hgnc.txt', 'r')
    for line in f:
        gene, drug, p_val = line.split()
        if gene not in shared_genes:
            continue
        p_val = float(p_val)
        if drug not in drug_top_genes_dct:
            drug_top_genes_dct[drug] = {}
            drug_top_genes_dct[drug][gene] = p_val
        elif len(drug_top_genes_dct[drug]) == TOP_GENES_PER_DRUG:
            continue
        else:
            drug_top_genes_dct[drug][gene] = p_val
    f.close()

    gene_pathway_lst = []
    f = open('./data/embedding/gene_pathway_id.txt', 'r')
    for line in f:
        # entity is either a gene or a pathway.
        gene_pathway_lst += [line.strip()]
    f.close()

    # Loop through the files.
    for num in [50, 100, 500, 1000, 1500, 2000]:
        num = str(num)
        for suffix in ['U', 'US']:
            entity_vector_dct = OrderedDict({})

            extension = '%s_0.8.%s' % (num, suffix)
            filename = './data/embedding/ppi_6_net_%s' % extension
            f = open(filename, 'r')
            for i, line in enumerate(f):
                # There should be the same number of rows in the embedding files
                # as the mapping file.
                entity = gene_pathway_lst[i]
                # Skip genes and pathways we have not previously used.
                if entity not in shared_genes and entity not in pathways:
                    continue
                # Convert lines to floats.
                vector = map(float, line.split())
                entity_vector_dct[entity] = vector
            f.close()

            # Calculate the score for each drug-pathway pair.
            out = open('./results/embedding/top_pathways_%s' % extension, 'w')
            score_dct = {}
            out.write('filler\n')
            out.write('drug\tpath\tscore\n')
            for drug in drug_top_genes_dct:
                gene_p_val_dct = drug_top_genes_dct[drug]

                # Get the cosine similarity between each pathway vector and gene
                # vector.
                for pathway in pathways:
                    # Initialize the score.
                    drug_path_score = 0
                    pathway_vec = entity_vector_dct[pathway]

                    for gene in gene_p_val_dct:
                        gene_vec = entity_vector_dct[gene]
                        cos = 1 - spatial.distance.cosine(pathway_vec, gene_vec)
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