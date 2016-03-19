### Author: Edward Huang

import fisher_test
import file_operations
import operator
import sys

### This script takes the most principal component computed for the Pearson
### p-values (the superdrug), takes the top_k genes with lowest values from
### the superdrug, and then computes a Fisher's test between the genes and each
### pathway, and writes out the pathway p-values to file.

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print "Usage: " + sys.argv[0] + " top_k"
        exit(1)
    top_k = int(sys.argv[1])

    superdrug_gene_vector = []
    f = open('./results/superdrug_gene_correlation_values.txt', 'r')
    for line in f:
        superdrug_gene_vector += [float(line.strip())]
    f.close()
    # Get the largest values in the most principal component.
    sorted_gene_vector = sorted(superdrug_gene_vector,
        reverse=True)[:top_k]
    # Get the indices of the component corresponding to these indices.
    super_gene_indices = [i for i, e in enumerate(superdrug_gene_vector) if e
        in sorted_gene_vector]
    
    expression_genes = file_operations.get_exp_dct().keys()
    # Get the gene names corresponding to the smallest values.
    super_genes = set([expression_genes[i] for i in super_gene_indices])

    # Get the NCI pathway dictionary.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()

    gene_universe = nci_genes.union(expression_genes)

    pathway_p_values = {}
    # Compute the top pathways for the superdrug with Fisher's test.
    for path in nci_path_dct:
        path_genes = set(nci_path_dct[path])
        sup_and_path = len(super_genes.intersection(path_genes))
        sup_not_path = len(super_genes.difference(path_genes))
        path_not_sup = len(path_genes.difference(super_genes))
        neither = len(gene_universe) - len(super_genes.union(path_genes))
        assert neither == len((gene_universe.difference(
            super_genes)).difference(path_genes))
        
        f_table = [[sup_and_path, sup_not_path], [path_not_sup, neither]]
        ft = fisher_test.FishersExactTest(f_table)
        pathway_p_values[path] = ft.two_tail_p()

    sorted_paths = sorted(pathway_p_values.items(), key=operator.itemgetter(1))

    # Write out to file.
    out = open('./results/superdrug_pathway_p_values.txt', 'w')
    for path, p_val in sorted_paths:
        out.write('%s\t%g\n' % (path, p_val))
    out.close()