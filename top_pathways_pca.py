### Author: Edward Huang

from collections import OrderedDict
from sklearn.decomposition import PCA
from scipy.stats.stats import spearmanr
import numpy as np

### This script uses PCA to find the most important genes. For each pathway p,
### we will have a G*N gene expression matrix A. N is the number of patients, G
### is the number of genes in p. We appl PCA on A, giving us a G*G matrix B.
### We take the most principal eigenvector of B, X, which has length G.
### For each patient, we multiply X with his or her gene expression, getting a
### score for each patient-pathway pair. We compute the Spearman correlation
### between the scores of all patients and drug response of all patients of a
### drug D. This correlation is the final score of pathway p to drug D.

# Threshold for significant p-value to count.
LOW_P_VALUE = 1e-4
genes = []
# Keys are drugs, values are lists of drug responses across all patients.
drug_resp_dct = OrderedDict({})
# Keys are genes, values are lists of gene expression across all patients.
gene_exp_dct = OrderedDict({})

if __name__ == '__main__':
    # Convert the pathway genes to their indices in the coexpression network.
    gene_to_index_dct = {}
    medgene_f = open('./data/medGene.txt', 'r')
    for i, line in enumerate(medgene_f):
        gene = line.split()[0]
        gene_to_index_dct[gene] = i + 1
    medgene_f.close()

    # Gets the features vectors, which are the gene expressions across patients
    print 'Extracting the feature vector...'
    exp_file = open('./data/gene2medProbeExpr.txt', 'r')
    for i, line in enumerate(exp_file):
        if i == 0:
            continue
        line = line.split()
        gene, gene_exp_line = gene_to_index_dct[line[0].lower()], line[1:]
        num_patients = len(gene_exp_line)
        genes += [gene]
        gene_exp_dct[gene] = map(float, gene_exp_line)
    exp_file.close()

    # Get the drug responses from the spreadsheet file.
    print 'Extracting drug response labels...'
    resp_file = open('./data/auc.tsv', 'r')
    for i, line in enumerate(resp_file):
        if i == 0:
            continue
        # Each row is one drug's performance on each patient.
        line = line.split()
        drug, resp_line = line[0], line[1:]
        resp_line = [None if resp=='NA' else float(resp) for resp in resp_line]
        drug_resp_dct[drug] = resp_line
    resp_file.close()

    print 'Extracting NCI pathways...'
    path_file = open('./data/nci_pathway.txt', 'r')
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

    # Translate the pathway genes to their corresponding indices.
    for path in nci_path_dct:
        for i, gene in enumerate(nci_path_dct[path]):
            gene = gene.lower()
            if gene in gene_to_index_dct:
                nci_path_dct[path][i] = gene_to_index_dct[gene]
        nci_path_dct[path] = [x for x in nci_path_dct[path] if type(x) == int]

    # Build the matrix for each pathway, and find the pathway's most principal
    # components
    most_principal_component_dct = {}
    drug_pathway_score_dct = OrderedDict({})
    num_low_p_first = 0 # Number of p-values below threshold for first component
    num_low_p_sec = 0 # Number of significant p-values for second component.
    print 'Computing PCA...'
    for path in nci_path_dct:
        patient_scores_first = []
        patient_scores_second = []
        pca_matrix = []
        for gene in nci_path_dct[path]:
            # Add in the gene expression values for the patients across
            # the pathway genes.
            pca_matrix += [gene_exp_dct[gene]]
        pca_matrix = np.array(pca_matrix).transpose()
        pca = PCA()
        pca.fit(pca_matrix)

        for i in range(len(pca_matrix)):
            # pca.components_ is an array. 0th element is the most principal
            # component.
            patient_scores_first += [np.dot(pca_matrix[i], pca.components_[0])]
            if len(pca.components_) > 1:
                patient_scores_second += [np.dot(pca_matrix[i], pca.components_[1])]
        for drug in drug_resp_dct:
            NA_indices = [i for i, e in enumerate(drug_resp_dct[drug]) if e == None]
            drug_resps = [e for i, e in enumerate(drug_resp_dct[drug]) if i not in NA_indices]
            
            temp_scores_first = [e for i, e in enumerate(patient_scores_first) if i not in NA_indices]
            rcc_f, p_value_f = spearmanr(temp_scores_first, drug_resps)
            p_value_f = abs(p_value_f)

            # Repeat if we wish to also use the second most principal component.
            if len(pca.components_) > 1:
                temp_scores_second = [e for i, e in enumerate(patient_scores_second) if i not in NA_indices]
                rcc_s, p_value_s = spearmanr(temp_scores_second, drug_resps)
                p_value_s = abs(p_value_s)
                
            if drug not in drug_pathway_score_dct:
                drug_pathway_score_dct[drug] = OrderedDict({})
            # If we have a p-value below threshold, count it.
            if p_value_f <= LOW_P_VALUE:
                num_low_p_first += 1
            if len(pca.components_) > 1:
                drug_pathway_score_dct[drug][path] = (p_value_f, p_value_s)
                # Count number of low p-values for second principal component.
                if p_value_s <= LOW_P_VALUE:
                    num_low_p_sec += 1
            else:
                drug_pathway_score_dct[drug][path] = [p_value_f]
    # Output the results to file.
    out = open('./results/top_pathways_pca.txt', 'w')
    out.write('num_low_p_first\t%d\tnum_low_p_second\t%d\n' % (num_low_p_first, num_low_p_sec))
    out.write('drug\tpath\tp-value\n')
    for drug in drug_pathway_score_dct:
        for path in drug_pathway_score_dct[drug]:
            for i, score in enumerate(drug_pathway_score_dct[drug][path]):
                # The i+1 is to differentiate between 1st and 2nd principal
                # components.
                out.write('%s\t%s_%d\t%s\n' % (drug, path, i+1, str(score)))
    out.close()