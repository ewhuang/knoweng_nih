### Author: Edward Huang

from collections import OrderedDict
import file_operations
import fisher_test
import operator
import time

### Gets the top pathways for each drug/cell-line using the LINCS data set.
### Usage: python top_pathways_lincs.py
### Run time: 42 minutes

Z_SCORE_MIN = 2
MAX_GENES_PER_DRUG = 250
LOW_P_THRESHOLD = 0.0001 # Count how many pathway-drug pairs are below this.
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

        # # Convert to floats again because np.transpose() changes to strings.
        # z_scores = map(float, z_scores)
        assert type(z_scores[0]) == float

        if drug in temp_drug_matrix:
            temp_drug_matrix[drug] += [z_scores]
        else:
            temp_drug_matrix[drug] = [z_scores]
    drug_matrix = temp_drug_matrix

    for drug in drug_matrix:
        z_scores = drug_matrix[drug]
        # Average the z-scores for different experiments of the same drug-cell
        # line pair.
        # z_scores = [np.mean(x) for x in zip(*z_scores)]
        z_scores = [mean(x) for x in zip(*z_scores)]
        top_gene_indices = []
        for i, z_score in enumerate(z_scores):
            if z_score >= Z_SCORE_MIN:
                top_gene_indices += [(i, z_score)]
        # Take only at most MAX_GENES_PER_DRUG.
        top_gene_indices = sorted(top_gene_indices, key=lambda k:k[1],
            reverse=True)[:MAX_GENES_PER_DRUG]
        top_genes = [lincs_genes[i] for i, z_score in top_gene_indices]
        drug_matrix[drug] = top_genes

    out = open('./results/top_pathways_lincs_diff_normalize_DMSO.txt', 'w')
    # Fisher's test for every drug/cell-line and path pair.
    total_num_genes = len(nci_genes.union(lincs_genes))
    fish_dct = {}
    num_low_p = 0
    for drug in drug_matrix:
        # Genes with the top z-scores for each drug, up to MAX_GENES_PER_DRUG.
        corr_genes = set(drug_matrix[drug])
        for path in nci_path_dct:
            path_genes = set(nci_path_dct[path])

            # Get the four relevant numbers for Fisher's test.
            corr_and_path = len(corr_genes.intersection(path_genes))
            corr_not_path = len(corr_genes.difference(path_genes))
            path_not_corr = len(path_genes.difference(corr_genes))
            neither = total_num_genes - len(corr_genes.union(path_genes))
            # Compute Fisher's test.
            f_table = ([[corr_and_path, corr_not_path],
                [path_not_corr, neither]])
            ft = fisher_test.FishersExactTest(f_table)
            p_value = ft.two_tail_p()
            if p_value < LOW_P_THRESHOLD:
                num_low_p += 1
            fkey = (drug, path, corr_and_path, corr_not_path, path_not_corr)
            fish_dct[fkey] = p_value
    sorted_fisher = sorted(fish_dct.items(), key=operator.itemgetter(1))
    # Write how many drug-pathway pairs have a very low p-value, determined by
    # LOW_P_THRESHOLD.
    out.write('num_p_below_%s\t%d\n' % (str(LOW_P_THRESHOLD), num_low_p))
    # Write the drug's top pathways to file.
    out.write('drug\tcell_line\tpath\tp-value\tinter\tlincs\tpath\n')
    for fkey, score in sorted_fisher:
        drug, path, inter, corr_len, path_len = fkey
        drug, cell_line = drug.split('_')
        out.write('%s\t%s\t%s\t' % (drug, cell_line, path))
        out.write('%g\t%d\t%d\t%d\n' % (score, inter, corr_len, path_len))
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))