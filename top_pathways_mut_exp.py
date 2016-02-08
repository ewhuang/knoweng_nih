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
P_THRESHOLD = 0.05
MAX_GENES_PER_DRUG = 250
LOW_P_THRESHOLD = 0.0001

# Extract the NCI pathway data.
path_file = open('./data/nci_pathway_hgnc.txt', 'r')
nci_path_dct = {}
nci_genes = set([])
for line in path_file:
    line = line.split('\t')
    path_name, path_gene = line[0], line[1][:-2]
    nci_genes.add(path_gene)
    if path_name not in nci_path_dct:
        nci_path_dct[path_name] = [path_gene]
    else:
        nci_path_dct[path_name] += [path_gene]
path_file.close()

# Get the drug responses from the spreadsheet file.
# Keys are drugs, values are lists of drug responses across all patients.
drug_resp_dct = OrderedDict({})
resp_file = open('./data/auc_hgnc.tsv', 'r')
for i, line in enumerate(resp_file):
    if i == 0:
        continue
    # Each row is one drug's performance on each patient.
    line = line.split()
    drug, resp_line = line[0], line[1:]
    resp_line = [None if resp == 'NA' else float(resp) for resp in resp_line]
    drug_resp_dct[drug] = resp_line
resp_file.close()

def write_genes_pathways(data_dct, run):
    print 'Computing Pearson coefficients for ' + run + '...'
    all_top_genes = {}
    # Dictionary of the top drug-pathway pairs.
    top_pairs = {}
    num_low_p = 0
    # The set of all genes in either expression or mutation, as well as NCI.
    all_genes = nci_genes.union(data_dct.keys())
    for drug in drug_resp_dct:
        drug_top_genes = {}
        drug_resp = drug_resp_dct[drug]
        # List of indices of N/A values in our drug response table.
        NA_i = [i for i, e in enumerate(drug_resp) if e == None]
        drug_resp = [e for i, e in enumerate(drug_resp) if i not in NA_i]
        if len(drug_resp) == 0:
            continue
        # Finding the top genes for each drug.
        for gene in data_dct:
            gene_data = data_dct[gene]
            gene_data = [e for i, e in enumerate(gene_data) if i not in NA_i]
            # Find the pearson coefficient between these two lists.
            pcc, p_value = pearsonr(drug_resp, gene_data)
            if p_value < P_THRESHOLD:
                all_top_genes[(gene, drug)] = pcc
                drug_top_genes[gene] = p_value
        # These are the top correlated genes for each drug.
        top_genes = sorted(drug_top_genes.items(),
            key=operator.itemgetter(1))[:MAX_GENES_PER_DRUG]
        # Get just the genes now, without the corresponding scores.
        corr_genes = set([gene for gene, pcc in top_genes])
        # Compute the top pathways for each drug with Fisher's test.
        for path in nci_path_dct:
            path_genes = set(nci_path_dct[path])
            corr_and_path = len(corr_genes.intersection(path_genes))
            corr_not_path = len(corr_genes.difference(path_genes))
            path_not_corr = len(path_genes.difference(corr_genes))
            neither = len(all_genes) - len(corr_genes.union(path_genes))
            assert len((corr_genes.union(path_genes)).difference(all_genes)) == 0
            # o_r = odds ratio.
            o_r, p_value = fisher_exact([[corr_and_path, corr_not_path], 
                [path_not_corr, neither]])
            # Count the number of significant p-values.
            if p_value < LOW_P_THRESHOLD:
                num_low_p += 1
            top_pairs[(drug, path, corr_and_path, 
                corr_not_path, path_not_corr, neither)] = p_value
    top_paths = sorted(top_pairs.items(), key=operator.itemgetter(1))

    # Write out the results.
    path_out = open('./results/top_pathways_%s_hgnc.txt' % run, 'w')
    path_out.write('num_below_%f\t%d\n' % (LOW_P_THRESHOLD, num_low_p))
    path_out.write('drug\tpath\tscore\tinter\tcorr\tpath\tneither\n')    
    for info, p_val in top_paths:
        drug, path, inter, corr, path_len, neither = info
        combo = '%s\t%s\t%g\t%d\t%d\t' % (drug, path, p_val, inter, corr)
        combo += '%d\t%d\n' % (path_len, neither)
        path_out.write(combo)
    path_out.close()

    # Sort the top genes by value. Get the top genes.
    gene_out = open('./results/top_genes_%s_hgnc.txt' % run, 'w')
    print 'Writing top genes for ' + run + '...'
    all_top_genes = sorted(all_top_genes.items(), key=operator.itemgetter(1), reverse=True)
    for (gene, drug), score in all_top_genes:
        gene_out.write(gene + '\t' + drug + '\t' + str(score) + '\n') 
    gene_out.close()

if __name__ == '__main__':
    # Keys are genes, values are lists of gene expression across all patients.
    exp_dct = OrderedDict({})
    # mut_dct = OrderedDict({})

    print 'Extracting the gene expression vectors...'
    exp_file = open('./data/gene_expression_hgnc.tsv', 'r')
    for i, line in enumerate(exp_file):
        # Skip the header row.
        if i == 0:
            continue
        line = line.split()
        gene, exp_line = line[0], line[1:]
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