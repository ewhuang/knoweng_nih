# knoweng_nih
KNOWENG Project.
Author: Edward Huang

______________________________________
>>> python ensg_to_hgnc_conversion.py
Converts the Stuart data that Sheng sent to the old data, compatible with the 
scripts meant to run on Mayo clinic data.

______________________________________
Finds the top pathways for each drug, cell-line pair from LINCS dataset.
>>> python top_pathways_lincs.py AFT_NUM
AFT_NUM is the suffix for various p-values of LINCS z-score files.
Reads in the LINCS data for level 4, and then finds the top pathway, drug
and cell-line pairs by using Fisher's hypergeometric test. The top genes for
each drug are found by taking z-scores the absolute values of which are higher
than 2.0. Output file is ./results/top_pathways_lincs_Aft_AFT_NUM.txt.
Currently uses the Stuart data.

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
>>> python compare_methods_with_lincs.py AFT_NUM pca/l1.exp
Outputs file, compare_lincs_and_METHOD.txt, depending on the method of choice. 
Uses Fisher's method to find how similar each method's top pathways are with
LINC's top pathways.
In line 123, we choose which method we use, defined by the methods created at
the top of the script.

______________________________________
Run embedding_top_pathways.py to find the top pathways for the embedding method.

>>> python embedding_top_pathways.py

Same format as the other top pathways. However, we can tune the top k pathways
to keep.