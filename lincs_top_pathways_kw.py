### Author: Edward Huang

from collections import OrderedDict
import file_operations
import operator
from scipy.stats.mstats import kruskalwallis
import time

### Gets the top pathways for each drug/cell-line using the LINCS data set.
### Usage: python top_pathways_lincs.py
### Run time: 42 minutes

KRUSKAL_P_THRESH = 0.05
# Define -infinity
inf = float('-inf')

def mean(lst):
    return float(sum(lst)) / len(lst)

def get_lincs_genes():
    # Extract lincs_genes from all_map.txt, provided by Sheng.
    f = open('./data/all_map.txt', 'r')
    lincs_genes = []
    for i, line in enumerate(f):
        line = line.split()
        lincs_genes += [line[1]]
    f.close()
    return lincs_genes

def get_drugs_and_gene_matrix():
    f = open('./data/lvl4_Stuart_combinedPvalue_diff_normalize_DMSO.txt', 'r')
    drugs, gene_matrix = [], []
    for line in f:
        line = line.split()
        drug, raw_z_scores = line[0], line[1:]
        drugs += [drug]
        gene_matrix += [raw_z_scores]
    f.close()
    return drugs, gene_matrix

def main():
    # Read NCI pathway dictionary and the lincs_genes in the NCI pathways.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()

    lincs_genes = get_lincs_genes()

    drugs, gene_matrix = get_drugs_and_gene_matrix()

    # Transpose the gene matrix so we can filter out the lincs_genes.
    gene_matrix = zip(*gene_matrix)
    assert len(gene_matrix) == len(lincs_genes)
    
    gene_dct = OrderedDict({})
    gene_counter, num_genes = 0, len(lincs_genes)
    while gene_matrix != []:
        gene_z_scores = gene_matrix.pop(0)
        assert len(gene_z_scores) == len(drugs)
        gene = lincs_genes[gene_counter]
        if gene == '-666':
            gene_counter += 1
            continue
        gene_z_scores = [inf if score == '#NAME?' else abs(float(
            score)) for score in gene_z_scores]
        if gene not in gene_dct:
            gene_dct[gene] = gene_z_scores
        else:
            # If a gene appears twice, then we want to get the max values.
            for ci, value in enumerate(gene_z_scores):
                gene_dct[gene][ci] = max(gene_dct[gene][ci], value)
        gene_counter += 1

    # Update lincs_genes to be just the valid ones in our dictionary.
    lincs_genes = gene_dct.keys()

    # Transposing back the z-scores.
    drug_matrix = [drugs]
    for gene in gene_dct:
        drug_matrix += [gene_dct[gene]]
    gene_dct.clear()
    drug_matrix = zip(*drug_matrix)

    # Make a new dictionary, with keys as drugs, and values as lists of LINCS
    # z-scores.
    temp_drug_matrix = OrderedDict({})
    for row in drug_matrix:
        drug, z_scores = row[0], row[1:]
        assert len(z_scores) == len(lincs_genes)
        assert type(z_scores[0]) == float

        if drug in temp_drug_matrix:
            temp_drug_matrix[drug] += [z_scores]
        else:
            temp_drug_matrix[drug] = [z_scores]
    drug_matrix = temp_drug_matrix

    # Dictionary of the top drug-pathway pairs.
    top_drug_path_pairs = {}
    # Counts the number of drug-pathway pairs with KW p-value KRUSKAL_P_THRESH.
    num_low_p = 0
    # Program status variables.
    progress_counter = 0
    num_drugs = float(len(drug_matrix))

    for drug in [drug_matrix.keys()[0]]:
        z_scores = drug_matrix[drug]
        # Average across samples.
        z_scores = [mean(x) for x in zip(*z_scores)]

        # Finding the correlations between each gene and the drug.
        drug_gene_z_scores = {}
        for gene_index, gene in enumerate(lincs_genes):
            drug_gene_z_scores[gene] = z_scores[gene_index]

        # Compute the top pathways for each drug with Kruskal-Wallis test.
        for path in nci_path_dct:
            # Make a copy of drug-gene correlations and then remove the genes
            # that are in each pathway.
            non_pathway_z_scores = drug_gene_z_scores.copy()
            pathway_z_scores = []
            path_genes = nci_path_dct[path]
            # Split the correlations into genes in the pathways and genes that
            # aren't.
            for gene in path_genes:
                # Skip genes in NCI pathways but not in gene expression data.
                if gene not in lincs_genes:
                    continue
                pathway_z_scores += [non_pathway_z_scores[gene]]
                del non_pathway_z_scores[gene]

            non_pathway_z_scores = non_pathway_z_scores.values()
            h_stat, p_value = kruskalwallis(pathway_z_scores,
                non_pathway_z_scores)
            top_drug_path_pairs[(drug, path)] = p_value

            if p_value < KRUSKAL_P_THRESH:
                num_low_p += 1
        progress_counter += 1.0
        print 'Progress: %f%%' % (progress_counter / num_drugs * 100)

    top_paths = sorted(top_drug_path_pairs.items(), key=operator.itemgetter(1))

    # Write out the results.
    out = open('./results/top_pathways_lincs_diff_normalize_DMSO_kw.txt', 'w')
    path_out.write('num_below_%f\t%d\n' % (KRUSKAL_P_THRESH, num_low_p))
    path_out.write('drug\tcell_line\tpath\tscore\n')    
    for (drug, path), p_val in top_paths:
        drug, cell_line = drug.split('_')
        path_out.write('%s\t%s\t%s\t%g' % (drug, cell_line, path, p_val))
    path_out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))