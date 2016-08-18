# knoweng_nih
KNOWENG Project.
Author: Edward Huang

______________________________________
>>> python ensg_to_hgnc_conversion.py
Converts the Stuart data that Sheng sent to the old data, compatible with the 
scripts meant to run on Mayo clinic data.

______________________________________
Finds the top pathways for each drug, cell-line pair from LINCS dataset.
$ python lincs_top_pathways.py z_score_min max_genes_per_drug

Original tuning parameters:
Z_SCORE_MIN = 2
MAX_GENES_PER_DRUG = 250

$ python top_pathways_lincs_positive_control.py AFT_NUM
Generates drug-cell_line-pathway combinations to check correctness of LINCS 
data.

__________________NETPATH____________________
1.
Finds top pathways for the gene expression data.
$ python correlation_top_pathways.py
Uses Pearson to find the correlation between gene expression and drug response.
Also finds the most correlated drug-pathway pairs through Fisher's test.

2.
Here, we can decide whether to use abs(cos) * corr or just (cos * corr) when
computing a drug-pathway score (line 110).

$ python embedding_top_pathways.py METHOD top_k

Same format as the other top pathways. However, we can tune the top k pathways
to keep.

3. Testing NetPath
$ python embedding_top_pathways_just_cosine.py
Creates an embedding file, ppi_top_pathways_50_0.8.U_top_250_just_cosine.txt,
that uses scores of only cosine similarity, not cosine * correlation of top
k most correlated genes per drug.

$ python random_embedding_top_pathways.py
Instead of taking the top 250 most correlated genes for each drug, we get the
a random 250 genes. Only runs for ppi_top_pathways_50_0.8.U_top_250. Produces
random_ppi_top_pathways_50_0.8.U_top_250.txt.

$ python compute_p_values_netpath.py
Compares the two files created by the two scripts prior to this one.

______________________________________
4.
Compare the previous drug-pathway ranking methods (PCA, L1 linear regression,
Pearson, multiple embedding networks)
with LINCS as a baseline.
$ python compare_methods_with_lincs.py correlation/ppi/genetic/literome/
        sequence top_k
Compares how each method finds drug-pathway pairs to LINCS, the ground truth, by
using Fisher's test.

5.
$ python summary_comparison_with_lincs.py top_k

Prints out summary tables for different p-values how each method compares to the
baseline, Pearson correlation.

______________SUPERDRUG______________
Superdrug - finding the drug that is most representative of all drugs to find
genes that have good drug response for all drugs.

$ python superdrug_principal_component.py

Outputs a file containing the most principal component vector, with length
equal to the number of genes in our data.

$ python superdrug_top_pathways.py TOP_K

Outputs a file containing the pathway p-values when compared to the TOP_K genes
of the superdrug.

Then 
$ python compare_methods_with_lincs.py ppi/genetic/literome/sequence TOP_K
This removes pathways below a threshold for the superdrug. We want to remove
the pathways that are too similar to too many drugs.

Then
>>> summary_comparison_with_lincs.py TOP_K
This script computes the 4x4 summary table for all of our methods.

$ python subtract_superdrug_from_pathways.py exp p_val
Removes the superdrug pathways from the expression pathways that are below the
p_val arugment. Outputs top_pathways_exp_hgnc_subtract_superdrug.txt.

Running random pathways.
$ python random_control_genes.py
Randomly samples pathways for each drug, where the number of samples equals
the number of pathway-drug pairs for expression below a certain threshold, with
this threshold in the range of [0.001, 0.005, 0.01, 0.05, 0.1].

$ python reformat_superdrug_genes.py
Gets the top 250 superdrug genes.

_____KRUSKAL-WALLIS TESTS_____
Recompute the top pathway rankings using the Kruskal-Wallis H-test instead of
Fisher's exact test.

$ python kw_top_pathways.py

First, preprocess LINCS file.
$ python preprocess_lincs_z_scores.py

$ python lincs_top_pathways_kw.py

Embedding doesn't use Fisher's test.

_____VARIOUS SCRIPTS_____
$ python kegg_lincs_intersection.py
Finds how many pathways are in common between KEGG pathway genes and LINCS
genes.

For the paper.
This script finds the drugs the NetPath has below the threshold that Fisher's
test does not retrieve for gene expression values. PPI, for 0.01 and 0.005
find 13 drug-path pairs. Expression only finds 10.

$ python various_paper_scripts.py