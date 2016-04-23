### Author: Edward Huang

import file_operations
from scipy.stats.mstats import kruskalwallis
from scipy.stats.stats import pearsonr
import operator
import time

### Find top genes from gene expression and mutation data sets.
### Uses Krusakl-Wallis instead of Fisher's test to then use the gene rankings
### to find most similar pathways.
### Run time: 3.4 hours.

# The maximum p-value to for Pearson's between drug response and gene expression
# to allow to be a significantly correlated gene for a drug.
KRUSKAL_P_THRESH = 0.05

# Extract the NCI pathway data.
nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
# Get the drug response dictionary.
drug_resp_dct = file_operations.get_drug_resp_dct()

def write_genes_pathways(exp_dct, method):
    # Dictionary of the top drug-pathway pairs.
    top_drug_path_pairs = {}
    # Counts the number of drug-pathway pairs with KW p-value KRUSKAL_P_THRESH.
    num_low_p = 0
    # Program status variables.
    progress_counter = 0
    num_drugs = float(len(drug_resp_dct))

    for drug in drug_resp_dct:
        print 'Progress: %f%%' % (progress_counter / num_drugs * 100)
        
        # Clean drug responses for each drug.
        drug_resp = drug_resp_dct[drug]
        # Indices of None values in our drug response table.
        NA_i = [i for i, e in enumerate(drug_resp) if e == None]
        # Get rid of None elements.
        drug_resp = [e for i, e in enumerate(drug_resp) if i not in NA_i]
        if len(drug_resp) == 0:
            assert [ele == None for ele in drug_resp_dct[drug]]
            continue

        # Finding the correlations between each gene and the drug.
        drug_gene_correlations = {}
        for gene in exp_dct:
            # Gene data is either gene expression values, or mutation counts.
            gene_data = exp_dct[gene]
            # Remove indices for elements that were None in drug response.
            gene_data = [e for i, e in enumerate(gene_data) if i not in NA_i]
            # Find the pearson coefficient between these two lists.
            pcc, p_value = pearsonr(drug_resp, gene_data)
            drug_gene_correlations[gene] = p_value

        # Compute the top pathways for each drug with Kruskal-Wallis test.
        for path in nci_path_dct:
            # Make a copy of drug-gene correlations and then remove the genes
            # that are in each pathway.
            non_pathway_correlations = drug_gene_correlations.copy()
            pathway_correlations = []
            path_genes = nci_path_dct[path]
            # Split the correlations into genes in the pathways and genes that
            # aren't.
            for gene in path_genes:
                # Skip genes in NCI pathways but not in gene expression data.
                if gene not in exp_dct:
                    continue
                pathway_correlations += [non_pathway_correlations[gene]]
                del non_pathway_correlations[gene]

            non_pathway_correlations = non_pathway_correlations.values()
            h_stat, p_value = kruskalwallis(pathway_correlations,
                non_pathway_correlations)
            top_drug_path_pairs[(drug, path)] = p_value

            if p_value < KRUSKAL_P_THRESH:
                num_low_p += 1

        progress_counter += 1.0

    top_paths = sorted(top_drug_path_pairs.items(), key=operator.itemgetter(1))

    # Write out the results.
    path_out = open('./results/top_pathways_%s_kw.txt' % method, 'w')
    path_out.write('num_below_%f\t%d\n' % (KRUSKAL_P_THRESH, num_low_p))
    path_out.write('drug\tpath\tscore\n')    
    for (drug, path), p_val in top_paths:
        path_out.write('%s\t%s\t%g\n' % (drug, path, p_val))
    path_out.close()

def main():
    # Getting gene expression dictionary.
    exp_dct = file_operations.get_exp_dct()

    # Write the top pathways for gene expression and mutation.
    write_genes_pathways(exp_dct, 'exp')

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))