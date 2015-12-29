### Author: Edward Huang

from collections import OrderedDict
from scipy.stats import fisher_exact
import operator
from sklearn.metrics import f1_score
import math
import sys

### Takes the drug-pathway scores from either PCA or linear regression L1
### and ranks by p-values. For either method, we arrive at a ranking of
### the most highly correlated pathways for each drug. We compare this to
### level 4 LINCS data by finding the size of the intersection of genes that
### have p-values below 0.0001, 0.1, and 0.04.

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print "Usage: " + sys.argv[0] + " pca/l1/exp"
        exit(1)
    method = sys.argv[1]


    out = open('./results/compare_lincs_and_%s.txt' % method, 'w')
    out.write('p-value threshold\tsize_intersection\tlincs\t%s' % method)
    out.write('\tneither\tfisher\'s p-value\n')
    # Results Processing
    for LOW_P_THRESHOLD in [1e-04, 0.01, 0.05]:
        # print 'Extracting data from %s output file...' % method
        results_dct = OrderedDict({})
        res_below_low_p = set([])
        if method == 'pca':
            results_file = open('./results/top_pathways_pca.txt', 'r')
        elif method == 'l1':
            results_file = open('./results/linear_regression_L1.txt', 'r')
        elif method == 'exp':
            results_file = open('./results/top_pathways_exp.txt', 'r')
        for i, line in enumerate(results_file):
            if method == 'pca' or method == 'exp':
                # Skip the first two line, as it contains summary of the results.
                if i < 2:
                    continue
                drug, path, score = line.strip().split('\t')
                # Remove differentiation between 1st and 2nd principal components.
                path = path[:-2]
            elif method == 'l1':
                drug, path, dont_use_this, score = line.split()
            if score == '[]':
                continue
            score = float(score)
            if score < LOW_P_THRESHOLD:
                res_below_low_p.add((drug, path))
            if (drug, path) not in results_dct:
                results_dct[(drug, path)] = score
            else:
                # For PCA, get the min value between the 1st and 2nd components.
                assert method == 'pca'
                results_dct[(drug, path)] = min(results_dct[(drug, path)], score)
        results_file.close()

        # print 'Extracting level 4 LINCS top pathways...'
        f = open('./results/top_pathways_lincs_lvl4.txt', 'r')
        lincs_below_low_p = set([])
        for i, line in enumerate(f):
            if i < 2:
                continue
            drug, cell_line, path, score, x, y, z = line.strip().split('\t')
            if float(score) < LOW_P_THRESHOLD:
                lincs_below_low_p.add((drug, path))
        f.close()

        # This just finds the overlap, but doesn't use hypergeometric.
        # print lincs_below_low_p
        # print res_below_low_p
        num_drug_pathway_pairs = len(results_dct) - len(lincs_below_low_p.union(res_below_low_p))
        lincs_and_res = len(lincs_below_low_p.intersection(res_below_low_p))
        lincs_not_res = len(lincs_below_low_p.difference(res_below_low_p))
        res_not_lincs = len(res_below_low_p.difference(lincs_below_low_p))
        # print method, 'p-value threshold: ' LOW_P_THRESHOLD
        # print lincs_and_res, lincs_not_res, res_not_lincs, num_drug_pathway_pairs
        o_r, p_val = fisher_exact([[lincs_and_res, lincs_not_res],
            [res_not_lincs, num_drug_pathway_pairs]])

        out.write('%f\t%d\t%d\t' % (LOW_P_THRESHOLD, lincs_and_res, lincs_not_res))
        out.write('%d\t%d\t%g\n' % (res_not_lincs, num_drug_pathway_pairs, p_val))

    out.close()