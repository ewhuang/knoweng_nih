### Author: Edward Huang

# from sklearn.decomposition import PCA
from collections import OrderedDict

### This script finds the superdrug by PCA's first principal component, and 
### runs Fisher's test between its top 250 genes and each pathway. Makes a file
### with the p-values.

if __name__ == '__main__':
    # Extract the NCI pathway data.
    path_file = open('./data/nci_pathway_hgnc.txt', 'r')
    nci_path_dct = {}
    # Set of all genes that appear in all NCI pathways.
    nci_genes = set([])
    for line in path_file:
        line = line.strip().split('\t')
        assert len(line) == 2
        path_name, path_gene = line
        nci_genes.add(path_gene)
        if path_name in nci_path_dct:
            nci_path_dct[path_name] += [path_gene]
        else:
            nci_path_dct[path_name] = [path_gene]
    path_file.close()

    # Get the drug responses from the spreadsheet file.
    # Keys are drugs, values are lists of drug responses across all patients.
    drug_resp_dct = OrderedDict({})
    resp_file = open('./data/auc_hgnc.tsv', 'r')
    for i, line in enumerate(resp_file):
        # Header line contains patient ID's.
        if i == 0:
            drug_response_patients = line.split()[1:]
            continue
        # Each row is one drug's performance on each patient.
        line = line.split()
        drug, resp_line = line[0], line[1:]
        assert len(resp_line) == len(drug_response_patients)
        assert 'BRD-' in drug
        # Convert 'NA' strings to None values.
        resp_line = [None if resp == 'NA' else float(resp) for resp in resp_line]
        drug_resp_dct[drug] = resp_line
    resp_file.close()

    # Getting gene expression vectors.
    exp_dct = OrderedDict({})
    exp_file = open('./data/gene_expression_hgnc.tsv', 'r')
    for i, line in enumerate(exp_file):
        # Header row contains patient ID's.
        if i == 0:
            expression_patients = line.split()[1:]
            assert drug_response_patients == expression_patients
            continue
        line = line.split()
        gene, exp_line = line[0], line[1:]
        assert len(exp_line) == len(expression_patients)
        exp_dct[gene] = map(float, exp_line)
    exp_file.close()

    num_drugs = float(len(drug_resp_dct))
    for drug in drug_resp_dct:
        drug_top_correlated_genes = {}
        # These are lists of drug responses for the drug.
        drug_resp = drug_resp_dct[drug]
        # Indices of None values in our drug response table.
        NA_i = [i for i, e in enumerate(drug_resp) if e == None]
        # Get rid of None elements.
        drug_resp = [e for i, e in enumerate(drug_resp) if i not in NA_i]
        if len(drug_resp) == 0:
            assert [ele == None for ele in drug_resp_dct[drug]]
            continue
        # Finding the most correlated genes for the drug.
        for gene in data_dct:
            # Gene data is either gene expression values, or mutation counts.
            gene_data = data_dct[gene]
            # Remove indices for elements that were None in drug response.
            gene_data = [e for i, e in enumerate(gene_data) if i not in NA_i]
            # Find the pearson coefficient between these two lists.
            pcc, p_value = pearsonr(drug_resp, gene_data)