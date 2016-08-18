### Author: Edward Huang

import sys
import file_operations
from scipy.stats import fisher_exact
import time

### Takes the drug-pathway scores from either PCA or linear regression L1
### and ranks by p-values. For either method, we arrive at a ranking of
### the most highly correlated pathways for each drug. We compare this to
### level 4 LINCS data by finding the size of the intersection of genes that
### have p-values below 0.0001, 0.001, 0.01, and 0.05. Computes the Fisher's
### test on these intersections, differentiating between cell lines and drugs.
### Run time: 1 minute.

embedding_methods = ['ppi', 'genetic', 'literome', 'sequence']
p_thresh_range = [0.001, 0.005, 0.01, 0.05, 0.1]
corr_fname = './results/top_pathways_correlation_hgnc.txt'

# Superdrug pathway p-values.
# sppv = file_operations.get_superdrug_pathway_p_values()
lincs_drug_path_dct = file_operations.get_lincs_drug_path_dct()
pathways = set([])

def get_corr_path_dct(method_p_thresh):
    '''
    Returns a dictionary.
    Each key is a drug.
    Each value is the set of paths most associated with the drug via the
    correlation baseline method.
    '''
    corr_path_dct = {}

    corr_file = open(corr_fname, 'r')
    for i, line in enumerate(corr_file):
        if i < 2:
            continue
        drug, path, p_value = line.strip().split('\t')[:3]
        pathways.add(path)

        if p_value == '[]' or float(p_value) > method_p_thresh:
            continue

        if drug not in corr_path_dct:
            corr_path_dct[drug] = set([path])
        else:
            corr_path_dct[drug].add(path)
    corr_file.close()

    return corr_path_dct

def get_embedding_path_dct(in_filename, method_p_thresh, corr_path_dct):
    '''
    Given the correlation method path dictionary, find the top pathways for
    each drug as computed by embedding, where the number of top pathways for
    each drug is the same as the number as computed by correlation.
    '''
    method_drug_dct = {}
    for drug in corr_path_dct:
        method_drug_dct[drug] = set([])

    method_file = open(in_filename, 'r')
    for i, line in enumerate(method_file):
        if i < 2:
            continue
        drug, path, score = line.strip().split('\t')[:3]
        pathways.add(path)
        if score == '[]' or drug not in corr_path_dct:
            continue

        # if drug not in method_drug_dct:
        #     method_drug_dct[drug] = set([path])
        # elif drug in method_drug_dct:
        #     # Go to next path if we have reached the same number of paths as we
        #     # have for the co-expression data.
        if len(corr_path_dct[drug]) == len(method_drug_dct[drug]):
            continue
        method_drug_dct[drug].add(path)
    method_file.close()
    return method_drug_dct

def get_top_pathways(in_filename, method_p_thresh):
    '''
    Returns a dictionary.
    Each key is a drug. Each value is a set of paths associated with the drug.
    If the in_filename is not the correlation pathway file, then we first
    process that to get the number of pathways per drug. We keep only the top
    k pathways per drug, where k is the same number of top pathways for the drug
    via the correlation method.
    '''
    # Simply return the correlation dictionary if it's the desired method.
    corr_path_dct = get_corr_path_dct(method_p_thresh)
    if in_filename == corr_fname:
        return corr_path_dct
    method_drug_dct = get_embedding_path_dct(in_filename, method_p_thresh,
        corr_path_dct)
    return method_drug_dct

def compare_methods(method, in_filename, out_filename):
    out = open(out_filename, 'w')
    out.write('method_p\tlincs_p\tdrug\tinter\tlincs\t%s' % method)
    out.write('\tneither\to_r\tfish-p\n')
    # Results Processing
    for method_p_thresh in p_thresh_range:
        res_drug_dct = get_top_pathways(in_filename, method_p_thresh)

        for lincs_p_thresh in p_thresh_range:
            lincs_drug_dct = {}
            for (drug, path) in lincs_drug_path_dct:
                pathways.add(path)
                if drug not in res_drug_dct:
                    continue
                if lincs_drug_path_dct[(drug, path)] <= lincs_p_thresh:
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

                f_table = ([[lincs_and_res, lincs_not_res],
                    [res_not_lincs, neither]])
                o_r, p_val = fisher_exact(f_table)
                # ft = fisher_test.FishersExactTest(f_table)
                # p_val = ft.two_tail_p()

                out.write('%f\t%f\t%s\t%d\t%d\t' % (method_p_thresh,
                    lincs_p_thresh, drug, lincs_and_res, lincs_not_res))
                out.write('%d\t%d\t%g\n' % (res_not_lincs, neither, p_val))
    out.close()

def main():
    if (len(sys.argv) < 3):
        print "Usage: " + sys.argv[0] + " correlation/kw/ppi/genetic/literome/sequence TOP_K"
        exit(1)
    method = sys.argv[1]
    top_k = int(sys.argv[2])

    assert (method in ['correlation', 'kw'] + embedding_methods)

    dimensions = map(str, [50, 100, 500])
    if method in embedding_methods:
        if method == 'ppi':
            dimensions += map(str, [1000, 1500, 2000])
        for dim in dimensions:
            for suffix in ['U', 'US']:
                extension = '%s_0.8.%s' % (dim, suffix)
                in_filename = './results/embedding/%s_top_pathways_%s_top_%d.txt' % (method, extension, top_k)
                out_filename = './results/lincs_comparison_files/compare_lincs_Aft_3_and_%s_%s_top_%d_hgnc.txt' % (method, extension, top_k)
                compare_methods(method, in_filename, out_filename)
    elif method == 'correlation':
        out_filename = './results/lincs_comparison_files/compare_lincs_Aft_3_and_correlation_hgnc.txt'
        compare_methods(method, corr_fname, out_filename)
    # corrression top pathways computed with Krusakl-Wallis
    elif method == 'kw':
        in_filename = './results/top_pathways_kw.txt'
        out_filename = './results/lincs_comparison_files/compare_lincs_Aft_3_and_kw.txt'
        compare_methods(method, in_filename, out_filename)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))