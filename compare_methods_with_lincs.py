### Author: Edward Huang

# from scipy.stats import fisher_exact, combine_pvalues
import math
import sys
import fisher_test

### Takes the drug-pathway scores from either PCA or linear regression L1
### and ranks by p-values. For either method, we arrive at a ranking of
### the most highly correlated pathways for each drug. We compare this to
### level 4 LINCS data by finding the size of the intersection of genes that
### have p-values below 0.0001, 0.001, 0.01, and 0.05. Computes the Fisher's
### test on these intersections, differentiating between cell lines and drugs.

AFT_NUM = 3

# Different aggregation functions.
def sum_log_p_values(p_val_lst):
    return sum([math.log(p_val, 10) for p_val in p_val_lst])

# def fishers_method(p_val_lst):
#     return combine_pvalues(p_val_lst, method='fisher')[1]

def min_p_exp(p_val_lst):
    return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))

def compare_methods(AFT_NUM, method, in_filename, out_filename):
    out = open(out_filename, 'w')
    out.write('p-value\tdrug\tinter\tlincs\t%s' % method)
    out.write('\tneither\to_r\tfish-p\n')
    # Results Processing
    for p_thresh in [0.0001, 0.001, 0.01, 0.05]:
        pathways = set([])

        # If the method is embedding, get the number of pathways per drug that
        # meet the desired threshold from gene expression.
        if method in ['ppi', 'l1', 'genetic', 'literome', 'sequence']:
            f = open('./results/top_pathways_exp_hgnc.txt', 'r')
            # We need to find out for each threshold, the number of pathways for
            # each drug-pathway pair that are better than the threshld.
            exp_paths_per_drug = {}
            for i, line in enumerate(f):
                # Skip header lines.
                if i < 2:
                    continue
                drug, path, score = line.strip().split('\t')[:3]

                if score == '[]':
                    continue

                # If a Fisher's p-value is below the threshold, add it to the
                # dictionary.
                if float(score) <= p_thresh:
                    if drug not in exp_paths_per_drug:
                        exp_paths_per_drug[drug] = 1
                    else:
                        exp_paths_per_drug[drug] += 1
            f.close()

        # Get the top pathways for each drug in our target method.
        res_drug_dct = {}
        results_file = open(in_filename, 'r')
        for i, line in enumerate(results_file):
            # Skip header lines.
            if i < 2:
                continue
            drug, path, score = line.strip().split('\t')[:3]
            
            if score == '[]':
                continue
            
            # Remove differentiation between 1st and 2nd principal components.
            if method == 'pca':
                path = path[:-2]            

            pathways.add(path)

            # These methods are the ones that need expression to know how many
            # pathways to use for each drug.
            non_exp = ['ppi', 'l1', 'genetic', 'literome', 'sequence']
            if method not in non_exp and float(score) <= p_thresh:
                if drug not in res_drug_dct:
                    res_drug_dct[drug] = set([path])
                else:
                    res_drug_dct[drug].add(path)
            elif method in non_exp:
                # Skip drugs not in expression data.
                if drug not in exp_paths_per_drug:
                    continue
                if drug not in res_drug_dct:
                    res_drug_dct[drug] = set([path])
                elif drug in res_drug_dct:
                    # Skip the drug if we have reached the same number of paths
                    # as we have for the expression data.
                    if exp_paths_per_drug[drug] == len(res_drug_dct[drug]):
                        continue
                    else:
                        res_drug_dct[drug].add(path)
        results_file.close()

        # print 'Extracting level 4 LINCS top pathways...'
        f = open('./results/top_pathways_lincs_Aft_%s_hgnc.txt' % AFT_NUM, 'r')
        lincs_drug_path_dct = {}
        for i, line in enumerate(f):
            # Skip header lines.
            if i < 2:
                continue

            drug, cell_line, path, score = line.strip().split('\t')[:4]

            if drug not in res_drug_dct:
                continue

            pathways.add(path)
            score = float(score)

            if (drug, path) not in lincs_drug_path_dct:
                lincs_drug_path_dct[(drug, path)] = [score]
            else:
                lincs_drug_path_dct[(drug, path)] += [score]
        f.close()

        lincs_drug_dct = {}
        # Aggregate the p-values across cell lines for each drug-pathway pair.          
        for (drug, path) in lincs_drug_path_dct:
            p_val_lst = lincs_drug_path_dct[(drug, path)]
            # Change the function for other aggregation functions.
            drug_path_p_val = min_p_exp(p_val_lst)
            if drug_path_p_val <= p_thresh:
                if drug not in lincs_drug_dct:
                    lincs_drug_dct[drug] = set([path])
                else:
                    lincs_drug_dct[drug].add(path)

        for drug in lincs_drug_dct:
            lincs_below_p = lincs_drug_dct[drug]
            res_below_p = res_drug_dct[drug]
            lincs_and_res = len(lincs_below_p.intersection(res_below_p))
            lincs_not_res = len(lincs_below_p.difference(res_below_p))
            res_not_lincs = len(res_below_p.difference(lincs_below_p))
            neither = len(pathways) - len(lincs_below_p.union(res_below_p))

            f_table = [[lincs_and_res, lincs_not_res], [res_not_lincs, neither]]
            ft = fisher_test.FishersExactTest(f_table)
            p_val = ft.two_tail_p()

            out.write('%f\t%s\t%d\t%d\t' % (p_thresh, drug, lincs_and_res,
                lincs_not_res))
            out.write('%d\t%d\t%g\n' % (res_not_lincs, neither, p_val))
    out.close()

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print "Usage: " + sys.argv[0] + " pca/exp/ppi/l1/genetic/literome/sequence TOP_K"
        exit(1)
    method = sys.argv[1]
    if len(sys.argv) == 3:
        top_k = int(sys.argv[2])

    assert (method in ['pca', 'exp', 'ppi', 'l1', 'genetic',
        'literome', 'sequence'])

    if method in ['ppi', 'genetic', 'literome', 'sequence']:
        if method == 'ppi':
            dimensions = map(str, [50, 100, 500, 1000, 1500, 2000])
        else:
            dimensions = map(str, [50, 100, 500])
        for dim in dimensions:
            for suffix in ['U', 'US']:
                extension = '%s_0.8.%s' % (dim, suffix)
                in_filename = './results/embedding/%s_top_pathways_%s_top_%d.txt' % (method, extension, top_k)
                out_filename = './results/compare_lincs_Aft_%s_and_%s_%s_top_%d_hgnc.txt' % (AFT_NUM, method, extension, top_k)
                compare_methods(AFT_NUM, method, in_filename, out_filename)
    else:
        in_filename = './results/top_pathways_%s_hgnc.txt' % method
        out_filename = './results/compare_lincs_Aft_%s_and_%s_hgnc.txt' % (AFT_NUM, method)
        compare_methods(AFT_NUM, method, in_filename, out_filename)