### Author: Edward Huang

import file_operations
import operator

### This script takes the vector of correlation coefficients after PCA for the
### superdrug genes, and prints them out with their corresponding gene names.


def get_superdrug_coefficients():
    superdrug_coefficients = []
    f = open('./results/superdrug_gene_correlation_values.txt', 'r')
    for line in f:
        superdrug_coefficients += [line.strip()]
    f.close()
    return superdrug_coefficients

def main():
    superdrug_coefficients = get_superdrug_coefficients()
    # Getting gene expression dictionary.
    exp_dct = file_operations.get_exp_dct()

    gene_coefficient_dct = {}
    for i, gene in enumerate(exp_dct.keys()):
        gene_coefficient_dct[gene] = float(superdrug_coefficients[i])

    sorted_dct = sorted(gene_coefficient_dct.items(), key=operator.itemgetter(1),
        reverse=True)[:250]

    out = open('./results/superdrug_top_250_genes.txt', 'w')
    for gene, coefficient in sorted_dct:
        out.write('%s\t%g\n' % (gene, coefficient))
    out.close()

if __name__ == '__main__':
    main()