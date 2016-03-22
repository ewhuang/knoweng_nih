### Author: Edward Huang

import file_operations
from scipy.stats.stats import pearsonr
from sklearn.decomposition import PCA

### This script finds the superdrug by PCA's first principal component, and then
### writes the component vector out to file.

if __name__ == '__main__':
    # Get the drug response dictionary.
    drug_resp_dct = file_operations.get_drug_resp_dct()

    # Getting gene expression dictionary.
    exp_dct = file_operations.get_exp_dct()

    num_drugs = float(len(drug_resp_dct))
    pca_matrix = []
    print "currently writing out pearson coefficients"
    for i, drug in enumerate(drug_resp_dct.keys()):
        print '%f%% done...' % (i / num_drugs * 100)
        # These are lists of drug responses for the drug.
        drug_resp = drug_resp_dct[drug]
        # Indices of None values in our drug response table.
        NA_i = [i for i, e in enumerate(drug_resp) if e == None]
        # Get rid of None elements.
        drug_resp = [e for i, e in enumerate(drug_resp) if i not in NA_i]
        if len(drug_resp) == 0:
            assert [ele == None for ele in drug_resp_dct[drug]]
            continue
        # Finding the correlations between each gene and the drug.
        drug_gene_p_values = []
        for gene in exp_dct:
            gene_exp = exp_dct[gene]
            # Remove indices for elements that were None in drug response.
            gene_exp = [e for i, e in enumerate(gene_exp) if i not in NA_i]
            # Find the pearson coefficient between these two lists.
            pcc, p_value = pearsonr(drug_resp, gene_exp)
            drug_gene_p_values += [pcc]
        pca_matrix += [drug_gene_p_values]

    pca = PCA()
    pca.fit(pca_matrix)
    # Superdrug is most principal component.
    superdrug_genes = pca.components_[0]

    # Write out the superdrug to file.
    out = open('./results/superdrug_gene_correlation_values.txt', 'w')
    for gene in superdrug_genes:
        out.write(str(gene) + '\n')
    out.close()