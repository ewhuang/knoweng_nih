# knoweng_nih
KNOWENG Project.
Author: Edward Huang

______________________________________
Finds the top pathways for each drug, cell-line pair from LINCS dataset.
>>> python top_pathways_lincs.py 3
or
>>> python top_pathways_lincs.py 4
Reads in the LINCS data for level 3 or 4, and then finds the top pathway, drug
and cell-line pairs by using Fisher's hypergeometric test. The top genes for
each drug are found by taking z-scores the absolute values of which are higher
than 2.5. Output file is ./results/top_pathways_lincs_lvl4.txt.

______________________________________
Uses PCA to find the most highly correlated drug-pathway pairs.
>>> python pca_path_drug_scores.py
Outputs results into pca_pathway_scores.txt.
Uses Spearman. First line of output file contains the number of pathway-drug
pairs that have p-value below 1e-4 (0.0001) for both first and second most
principal components.

______________________________________
Compare the previous drug-pathway ranking methods (PCA, L1 linear regression)
with LINCS as a baseline.
>>> python compare_pca_l1_with_lincs.py pca/l1
Outputs file, compare_lincs_and_pca.txt or compare_lincs_and_l1.txt, depending
on the method of choice.