# knoweng_nih
###Author: Edward Huang
KnowEng project for NetPath.

## Preprocessing
Converts the Stuart data that Sheng sent to the old data, compatible with the 
scripts meant to run on Mayo clinic data.

```bash
$ python ensg_to_hgnc_conversion.py
```

## LINCS top pathways (ground truth)
1.  Finds the top pathways for each drug, cell-line pair from LINCS dataset.
    
    ```bash
    $ python lincs_top_pathways.py z_score_min max_genes_per_drug
    ```

    Original tuning parameters:
    Z_SCORE_MIN = 2
    MAX_GENES_PER_DRUG = 250

2.  Positive control for LINCS.

    ```bash
    $ python top_pathways_lincs_positive_control.py AFT_NUM
    ```
    
    Generates drug-cell_line-pathway combinations to check correctness of LINCS
    data.

## NetPath
1.  Finds top pathways for the gene expression data.

    ```bash
    <!-- $ python correlation_top_pathways_fisher.py top_k -->
    $ python2.7 drug_pathway_fisher_correlation.py -s sortP -p 0.05 -i ge
    $ python correlation_top_pathways_kw.py
    ```

    Makes two files. One shows p-values between drugs and pathways (Fisher's
    test between most highly correlated genes for the drug and genes in the
    pathway). Another shows p-values between drugs and genes (Pearson between
    gene expression and drug response).

2.  Make embedding files that only have gene embedding vectors. Kind of a 
    misnomer. Doesn't add drug embeddings, but rather the missing pathway
    embeddings.

    ```bash
    python impute_drugless_embeddings.py
    ```

2.  Here, we can decide whether to use abs(cos) * corr or just (cos * corr) when
    computing a drug-pathway score (line 110).

    ```bash
    $ python embedding_top_pathways.py top_k
    ```

    Same format as the other top pathways. However, we can tune the top k
    pathways to keep.

## Testing NetPath

1.  Creates an embedding file, ppi_top_pathways_50_0.8.U_top_250_just_cosine.txt,
    that uses scores of only cosine similarity, not cosine * correlation of top
    k most correlated genes per drug.

    ```bash
    $ python embedding_top_pathways_just_cosine.py
    ```

4.  Instead of taking the top 250 most correlated genes for each drug, we get
    the a random 250 genes. Only runs for ppi_top_pathways_50_0.8.U_top_250.
    Produces random_ppi_top_pathways_50_0.8.U_top_250.txt.

    ```bash
    $ python random_embedding_top_pathways.py
    ```
5.  Compares the two files created by the two scripts prior to this one.

    ```bash
    $ python compute_p_values_netpath.py
    ```

## NetPath Evaluation
1.  Compare the previous drug-pathway ranking methods with LINCS as a baseline.

    ```bash
    $ python compare_methods_with_lincs.py corr_fisher/corr_kw/ppi/genetic/literome/sequence lincs_z lincs_max_num top_k
    ```
        
    Compares how each method finds drug-pathway pairs to LINCS, the ground truth, 
    by using Fisher's test. lincs_z and lincs_max_num determine which LINCS top
    pathways file to compare to.

5.  Prints out summary tables for different p-values how each method compares to
    the baseline, Pearson correlation.

    ```bash
    $ python summary_comparison_with_lincs.py lincs_z lincs_max_num top_k
    ```

## Superdrug experiments
Superdrug - finding the drug that is most representative of all drugs to find
genes that have good drug response for all drugs.

1.  Outputs a file containing the most principal component vector, with length
    equal to the number of genes in our data.

    ```bash
    $ python superdrug_principal_component.py
    ```

2.  Outputs a file containing the pathway p-values when compared to the TOP_K
    genes of the superdrug.

    ```bash
    $ python superdrug_top_pathways.py TOP_K
    ```

3.  This removes pathways below a threshold for the superdrug. We want to remove
    the pathways that are too similar to too many drugs.

    ```bash
    $ python compare_methods_with_lincs.py ppi/genetic/literome/sequence TOP_K
    ```

4.  This script computes the 4x4 summary table for all of our methods.
    
    ```bash
    $ summary_comparison_with_lincs.py TOP_K
    ```

5.  Removes the superdrug pathways from the expression pathways that are below
    the p_val arugment. Outputs top_pathways_exp_hgnc_subtract_superdrug.txt.

    ```bash
    $ python subtract_superdrug_from_pathways.py exp p_val
    ```

6.  Running random pathways. Randomly samples pathways for each drug, where the
    number of samples equals the number of pathway-drug pairs for expression
    below a certain threshold, with this threshold in the range of [0.001,
    0.005, 0.01, 0.05, 0.1].

    ```bash
    $ python random_control_genes.py
    ```

7.  Gets the top 250 superdrug genes.

    ```bash
    $ python reformat_superdrug_genes.py
    ```

## Kruskal-Wallis baseline tests
1.  Recompute the top pathway rankings using the Kruskal-Wallis H-test instead
    of Fisher's exact test.

    ```bash
    $ python kw_top_pathways.py
    ```

2.  Preprocess LINCS file.

    ```bash
    $ python preprocess_lincs_z_scores.py
    ```

3.  Embedding doesn't use Fisher's test.

    ```bash
    $ python lincs_top_pathways_kw.py
    ```

## Assorted scripts
1.  Finds how many pathways are in common between KEGG pathway genes and LINCS
    genes.

    ```bash
    $ python kegg_lincs_intersection.py
    ```
2. For the paper. This script finds the drugs the NetPath has below the
    threshold that Fisher's test does not retrieve for gene expression values.
    PPI, for 0.01 and 0.005. Finds 13 drug-path pairs. Expression only finds 10.

    ```bash
    $ python various_paper_scripts.py
    ```