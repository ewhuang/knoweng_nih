### Author: Edward Huang

from scipy.stats import fisher_exact, combine_pvalues
import math
import sys

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

def fishers_method(p_val_lst):
    return combine_pvalues(p_val_lst, method='fisher')[1]

def min_p_exp(p_val_lst):
    return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))

def compare_methods(AFT_NUM, method, in_filename, out_filename):
    out = open(out_filename, 'w')
    out.write('p-value\tdrug\tinter\tlincs\t%s' % method)
    out.write('\tneither\to_r\tfish-p\n')
    # Results Processing
    for p_threshold in [0.0001, 0.001, 0.01, 0.05]:
        pathways = set([])
        res_below_low_p = set([])
        # For if we want to see the intersections by drug.
        res_drug_dct = {}

        # If the method is embedding, then we must first open the top pathways
        # from the expression data.
        if method in ['ppi', 'l1', 'genetic', 'literome', 'sequence']:
            f = open('./results/top_pathways_exp_hgnc.txt', 'r')
            # We need to find out for each threshold, the number of pathways for
            # each drug-pathway pair that are better than the threshld.
            exp_num_paths_per_drug_dct = {}
            for i, line in enumerate(f):
                # Skip header lines.
                if i < 2:
                    continue
                raw_string = line.strip().split('\t')
                drug, path, score = raw_string[0], raw_string[1], raw_string[2]
                
                if score == '[]':
                    continue

                # If a Fisher's p-value is below the threshold, add it to the
                # dictionary.
                if float(score) <= p_threshold:
                    if drug not in exp_num_paths_per_drug_dct:
                        exp_num_paths_per_drug_dct[drug] = 1
                    else:
                        exp_num_paths_per_drug_dct[drug] += 1
            f.close()

        results_file = open(in_filename, 'r')
        for i, line in enumerate(results_file):
            # Skip the first two lines, as they contain summary of the results.
            if i < 2:
                continue
            raw_string = line.strip().split('\t')
            drug, path, score = raw_string[0], raw_string[1], raw_string[2]
            
            if score == '[]':
                continue
            
            # Remove differentiation between 1st and 2nd principal components.
            if method == 'pca':
                path = path[:-2]            

            pathways.add(path)

            if method not in ['ppi', 'l1', 'genetic', 'literome', 'sequence'] and float(score) <= p_threshold:
                if drug not in res_drug_dct:
                    res_drug_dct[drug] = set([path])
                else:
                    res_drug_dct[drug].add(path)
            elif method in ['ppi', 'l1', 'genetic', 'literome', 'sequence']:
                # Skip drugs not in expression data.
                if drug not in exp_num_paths_per_drug_dct:
                    continue
                if drug not in res_drug_dct:
                    res_drug_dct[drug] = set([path])
                elif drug in res_drug_dct:
                    # Skip the drug if we have reached the same number of paths
                    # as we have for the expression data.
                    if exp_num_paths_per_drug_dct[drug] == len(res_drug_dct[drug]):
                        continue
                    else:
                        res_drug_dct[drug].add(path)

        results_file.close()

        # print 'Extracting level 4 LINCS top pathways...'
        f = open('./results/top_pathways_lincs_Aft_%s_hgnc.txt' % AFT_NUM, 'r')
        lincs_drug_path_dct = {}
        for i, line in enumerate(f):
            # Skip first two lines, they're just results summaries.
            if i < 2:
                continue

            # Disregard x, y, z.
            drug, cell_line, path, score, x, y, z = line.strip().split('\t')

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
            # Change the sum_log_p_values() function for other aggregation
            # functions.
            drug_path_p_val = min_p_exp(p_val_lst)
            if drug_path_p_val <= p_threshold:
                if drug not in lincs_drug_dct:
                    lincs_drug_dct[drug] = set([path])
                else:
                    lincs_drug_dct[drug].add(path)

        # Loop through every single cell line and find Fisher's with the the
        # results from our methods (pca or exp).
        # for cell_line in lincs_cell_line_dct:
        #     lincs_below_low_p = lincs_cell_line_dct[cell_line]
        #     lincs_and_res = len(lincs_below_low_p.intersection(res_below_low_p))
        #     lincs_not_res = len(lincs_below_low_p.difference(res_below_low_p))
        #     res_not_lincs = len(res_below_low_p.difference(lincs_below_low_p))
        #     neither = len(drug_path_pairs) - len(lincs_below_low_p.union(res_below_low_p))

        #     o_r, p_val = fisher_exact([[lincs_and_res, lincs_not_res],
        #         [res_not_lincs, neither]])

        #     out.write('%f\t%s\t%d\t%d\t' % (p_threshold, cell_line, lincs_and_res, lincs_not_res))
        #     out.write('%d\t%d\t%g\t%g\n' % (res_not_lincs, neither, o_r, p_val))

        for drug in lincs_drug_dct:
            if drug not in res_drug_dct:
                continue
            lincs_below_low_p = lincs_drug_dct[drug]
            res_below_low_p = res_drug_dct[drug]
            lincs_and_res = len(lincs_below_low_p.intersection(res_below_low_p))
            lincs_not_res = len(lincs_below_low_p.difference(res_below_low_p))
            res_not_lincs = len(res_below_low_p.difference(lincs_below_low_p))
            neither = len(pathways) - len(lincs_below_low_p.union(res_below_low_p))

            o_r, p_val = fisher_exact([[lincs_and_res, lincs_not_res],
                [res_not_lincs, neither]])

            # if o_r > 1.0 and p_val < 0.05:
            #     print drug, o_r, p_val, p_threshold

            out.write('%f\t%s\t%d\t%d\t' % (p_threshold, drug, lincs_and_res, lincs_not_res))
            out.write('%d\t%d\t%g\t%g\n' % (res_not_lincs, neither, o_r, p_val))
    out.close()

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print "Usage: " + sys.argv[0] + " pca/exp/ppi/l1/genetic/literome/sequence TOP_K"
        exit(1)
    # AFT_NUM = sys.argv[1]
    method = sys.argv[1]
    top_k = int(sys.argv[2])

    assert (method in ['pca', 'exp', 'ppi', 'l1', 'genetic', 'literome', 'sequence'])

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