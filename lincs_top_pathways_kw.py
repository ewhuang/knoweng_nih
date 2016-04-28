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

def main():
    # Read NCI pathway dictionary and the lincs_genes in the NCI pathways.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()

    # Counts the number of drug-pathway pairs with KW p-value KRUSKAL_P_THRESH.
    num_low_p = 0
    # Dictionary of the top drug-pathway pairs.
    top_drug_path_pairs = {}
    # Reading in preprocessed data.
    lincs_genes = []
    num_drugs = 9298.0
    drug_counter = 0.0

    f = open('./data/processed_lincs_normalized_DMSO.txt', 'r')
    for i, line in enumerate(f):
        line = line.split()
        # First line is the LINCS genes.
        if i == 0:
            lincs_genes = line
            continue
        # Otherwise, read in the drug and all of its expressed genes.
        drug, z_scores = line[0], line[1:]
        z_scores = map(float, z_scores)
        assert len(z_scores) == len(lincs_genes)

        # Compute the top pathways for each drug with Kruskal-Wallis test.
        for path in nci_path_dct:
            path_genes = nci_path_dct[path]
            # Split the into genes in the pathways and genes that aren't.
            pathway_z_scores, non_pathway_z_scores = [], []
            for gene_index, gene in enumerate(lincs_genes):
                if gene in path_genes:
                    pathway_z_scores += [z_scores[gene_index]]
                else:
                    non_pathway_z_scores += [z_scores[gene_index]]

            h_stat, p_value = kruskalwallis(pathway_z_scores,
                non_pathway_z_scores)
            top_drug_path_pairs[(drug, path, h_stat)] = p_value

            if p_value < KRUSKAL_P_THRESH:
                num_low_p += 1
        drug_counter += 1
        print 'Progres:%f%%' % (drug_counter / num_drugs * 100)
    f.close()

    top_paths = sorted(top_drug_path_pairs.items(), key=operator.itemgetter(1))

    # Write out the results.
    path_out = open('./results/top_pathways_lincs_diff_normalize_DMSO_kw.txt',
        'w')
    path_out.write('num_below_%f\t%d\n' % (KRUSKAL_P_THRESH, num_low_p))
    path_out.write('drug\tcell_line\tpath\th_statistic\tp_value\n')    
    for (drug, path, h_stat), p_val in top_paths:
        drug, cell_line = drug.split('_')
        path_out.write('%s\t%s\t%s\t%g\t%g\n' % (drug, cell_line, path,
            h_stat, p_val))
    path_out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))