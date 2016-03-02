### Author: Edward Huang

from scipy.stats.stats import pearsonr
from scipy.stats import fisher_exact
import operator
from collections import OrderedDict

### Find top genes from gene expression and mutation data sets.
### Uses Fisher's test to then use the gene rankings to find most similar
### pathways.

# The maximum p-value to for Pearson's between drug response and gene expression
# to allow to be a significantly correlated gene for a drug.
PEARSON_P_THRESH = 0.05
MAX_GENES_PER_DRUG = 250
FISHER_P_THRESH = 0.0001

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

def write_genes_pathways(data_dct, method):
    gene_drug_correlations = {}
    # Dictionary of the top drug-pathway pairs.
    top_drug_path_pairs = {}
    # Counts the number of drug-pathway pairs with Fisher p-value below
    # FISHER_P_THRESH.
    num_low_p = 0
    # The set of all genes in either expression or mutation, as well as NCI.
    gene_universe = nci_genes.union(data_dct.keys())
    progress_counter = 0
    num_drugs = float(len(drug_resp_dct))
    for drug in drug_resp_dct:
        print 'Progress: %f%%' % (progress_counter / num_drugs)
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
            if p_value < PEARSON_P_THRESH:
                gene_drug_correlations[(gene, drug)] = pcc
                drug_top_correlated_genes[gene] = p_value
        # These are the top correlated genes for the drug.
        top_genes = sorted(drug_top_correlated_genes.items(),
            key=operator.itemgetter(1))[:MAX_GENES_PER_DRUG]
        # Get just the genes now, without the corresponding scores.
        corr_genes = set([gene for gene, pcc in top_genes])
        # Compute the top pathways for the drug with Fisher's test.
        for path in nci_path_dct:
            path_genes = set(nci_path_dct[path])
            corr_and_path = len(corr_genes.intersection(path_genes))
            corr_not_path = len(corr_genes.difference(path_genes))
            path_not_corr = len(path_genes.difference(corr_genes))
            neither = len(gene_universe) - len(corr_genes.union(path_genes))
            assert neither == len((gene_universe.difference(
                corr_genes)).difference(path_genes))
            # o_r = odds ratio.
            o_r, p_value = fisher_exact([[corr_and_path, corr_not_path], 
                [path_not_corr, neither]])
            # Count the number of significant p-values.
            if p_value < FISHER_P_THRESH:
                num_low_p += 1
            top_drug_path_pairs[(drug, path, corr_and_path, 
                corr_not_path, path_not_corr, neither)] = p_value
        progress_counter += 1.0
    top_paths = sorted(top_drug_path_pairs.items(), key=operator.itemgetter(1))

    # Write out the results.
    path_out = open('./results/top_pathways_%s_hgnc.txt' % method, 'w')
    path_out.write('num_below_%f\t%d\n' % (FISHER_P_THRESH, num_low_p))
    path_out.write('drug\tpath\tscore\tinter\tcorr\tpath\tneither\n')    
    for key, p_val in top_paths:
        drug, path, inter, corr_len, path_len, neither = key
        string = '%s\t%s\t%g\t%d\t%d\t' % (drug, path, p_val, inter, corr_len)
        string += '%d\t%d\n' % (path_len, neither)
        path_out.write(string)
    path_out.close()

    # Sort the top genes by value. Get the top genes.
    gene_out = open('./results/top_genes_%s_hgnc.txt' % method, 'w')
    print 'Writing top genes for ' + method + '...'
    gene_drug_correlations = sorted(gene_drug_correlations.items(),
        key=operator.itemgetter(1), reverse=True)
    for (gene, drug), score in gene_drug_correlations:
        gene_out.write(gene + '\t' + drug + '\t' + str(score) + '\n') 
    gene_out.close()

if __name__ == '__main__':
    # Keys are genes, values are lists of gene expression across all patients.
    exp_dct = OrderedDict({})
    # mut_dct = OrderedDict({})

    print 'Extracting the gene expression vectors...'
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

    # print 'Extracting the mutation vectors...'
    # mut_file = open('./data/gene2SNPu50SumParse.txt', 'r')
    # for i, line in enumerate(mut_file):
    #     if i == 0:
    #         continue
    #     line = line.split()
    #     gene, mut_line = line[0], line[1:]
    #     mut_dct[gene] = map(int, mut_line)
    # mut_file.close()

    # Write the top pathways for gene expression and mutation.
    write_genes_pathways(exp_dct, 'exp')
    # write_genes_pathways(mut_dct, 'mut')