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
### level 4 LINCS data by using Fisher's exact test.

# RES_P_THRESHOLD = 0.05x
# LINCS_P_THRESHOLD = 0.1
# LOW_P_THRESHOLD = 0.05

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print "Usage: " + sys.argv[0] + " pca/l1/exp"
        exit(1)
    method = sys.argv[1]

    # Results Processing
    for LOW_P_THRESHOLD in [1e-04, 0.01, 0.05]:
        # print 'Extracting data from %s output file...' % method
        pathways = set([])
        results_dct = OrderedDict({})
        res_below_low_p = set([])
        if method == 'pca':
            results_file = open('./results/top_pathways_pca.txt', 'r')
        elif method == 'l1':
            results_file = open('./results/linear_regression_L1.txt', 'r')
        elif method == 'exp':
            results_file = open('./results/top_pathways_exp.txt', 'r')
        for i, line in enumerate(results_file):
            if method == 'pca':
                # Skip the first two line, as it contains summary of the results.
                if i < 2:
                    continue
                drug, path, score = line.strip().split('\t')
                # Remove differentiation between 1st and 2nd principal components.
                path = path[:-2]
            elif method == 'l1':
                drug, path, dont_use_this, score = line.split()
            elif method == 'exp':
                if i < 2:
                    continue
                line = line.strip().split('\t')
                drug, path, score = line[0], line[1], line[2]
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
            pathways.add(path)
        results_file.close()

        # print 'Sorting results by p values...'
        # sorted_res = sorted(results_dct.items(), key=operator.itemgetter(1))
        # res_dct = {}
        # for (drug, path), score in sorted_res:
        #     # Get only p-values lower than P_THRESHOLD to be in a drug's top paths.
        #     if score > RES_P_THRESHOLD:
        #         if drug not in res_dct:
        #             res_dct[drug] = []
        #         continue
        #     if drug not in res_dct:
        #         res_dct[drug] = [path]
        #     else:
        #         res_dct[drug] += [path]

        # print 'Extracting level 4 LINCS top pathways...'
        f = open('./results/top_pathways_lincs_lvl4.txt', 'r')
        lincs_dct = {}
        lincs_below_low_p = set([])
        for i, line in enumerate(f):
            if i < 2:
                continue
            drug, cell_line, path, score, x, y, z = line.strip().split('\t')
            size = len(lincs_dct)
            lincs_dct[(drug + '_' + cell_line, path)] = float(score)
            pathways.add(path)
            if float(score) < LOW_P_THRESHOLD:
                lincs_below_low_p.add((drug, path))
        f.close()

        # print 'Sorting level 4 LINCS top pathways...'
        # sorted_lincs = sorted(lincs_dct.items(), key=operator.itemgetter(1))
        # # lincs_top_dct = {}
        # for (drug_cell, path), score in sorted_lincs:
            # if score > LINCS_P_THRESHOLD:
            #     continue
            # if drug_cell not in lincs_top_dct:
            #     lincs_top_dct[drug_cell] = [path]
            # else:
            #     lincs_top_dct[drug_cell] += [path]

        # This just finds the overlap, but doesn't use hypergeometric.
        # print lincs_below_low_p
        # print res_below_low_p
        num_drug_pathway_pairs = len(results_dct) - len(lincs_below_low_p.union(res_below_low_p))
        lincs_and_res = len(lincs_below_low_p.intersection(res_below_low_p))
        lincs_not_res = len(lincs_below_low_p.difference(res_below_low_p))
        res_not_lincs = len(res_below_low_p.difference(lincs_below_low_p))
        print method, LOW_P_THRESHOLD
        print lincs_and_res, lincs_not_res, res_not_lincs, num_drug_pathway_pairs
        print fisher_exact([[lincs_and_res, lincs_not_res],
            [res_not_lincs, num_drug_pathway_pairs]])
        # exit()






    exit()






    # Compare each drug-cell line with pca pathways.
    fisher_dct = {}
    int_dct = {}
    for drug_cell in lincs_top_dct:
        drug, cell_line = drug_cell.split('_')
        lincs = set(lincs_top_dct[drug_cell])
        res = set(res_dct[drug])
        # if method == 'l1':
        #     res = set(res_dct[drug][:50])
        lincs_and_res = len(lincs.intersection(res))
        lincs_not_res = len(lincs.difference(res))
        res_not_lincs = len(res.difference(lincs))
        neither = len(pathways) - len(res.union(lincs))
        o_r, p_value = fisher_exact([[lincs_and_res, lincs_not_res], [res_not_lincs, neither]])
        fisher_dct[(drug, cell_line)] = p_value
        int_dct[(drug, cell_line)] = '%d\t%d\t%d' % (lincs_and_res, len(lincs), len(res))
    sorted_fisher = sorted(fisher_dct.items(), key=operator.itemgetter(1))

    print 'Writing results out to file...'
    out = open('./results/compare_lincs_and_' + method + '.txt', 'w')
    out.write('drug\tcell_line\tp_val\tintersection\tlincs\t%s\n' % method)
    for (drug, cell_line), p_val in sorted_fisher:
        out.write('%s\t%s\t%g\t%s\n' % (drug, cell_line, p_val, int_dct[drug, cell_line]))
    out.close()