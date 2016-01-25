# knoweng_nih
KNOWENG Project.
Author: Edward Huang

______________________________________
Finds the top pathways for each drug, cell-line pair from LINCS dataset.
>>> python top_pathways_lincs.py AFT_NUM
AFT_NUM is the suffix for various p-values of LINCS z-score files.
Reads in the LINCS data for level 4, and then finds the top pathway, drug
and cell-line pairs by using Fisher's hypergeometric test. The top genes for
each drug are found by taking z-scores the absolute values of which are higher
than 2.0. Output file is ./results/top_pathways_lincs_Aft_AFT_NUM.txt.

>>> python top_pathways_lincs_positive_control.py AFT_NUM
Generates drug-cell_line-pathway combinations to check correctness of LINCS 
data.

______________________________________
Finds top pathways for the gene expression data.
>>> python top_pathways_mut_exp.py
Uses Pearson to find the correlation between gene expression and drug response.
Outputs two files, top_pathways_exp.txt and top_pathways_mut.txt.
They show the top drug-pathway pairs through Fisher's test.

______________________________________
Uses PCA to find the most highly correlated drug-pathway pairs.
>>> python top_pathways_pca.py
Outputs results into pca_pathway_scores.txt.
Uses Spearman. First line of output file contains the number of pathway-drug
pairs that have p-value below 1e-4 (0.0001) for both first and second most
principal components.

______________________________________
Compare the previous drug-pathway ranking methods (PCA, L1 linear regression)
with LINCS as a baseline.
>>> python compare_methods_with_lincs.py pca/l1.exp
Outputs file, compare_lincs_and_METHOD.txt, depending on the method of choice. 
Uses Fisher's test to find how similar each method's top pathways are with
LINC's top pathways.