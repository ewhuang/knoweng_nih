### Author: Edward Huang

from scipy.stats import fisher_exact
import sys

### Takes the drug-pathway scores from either PCA or linear regression L1
### and ranks by p-values. For either method, we arrive at a ranking of
### the most highly correlated pathways for each drug. We compare this to
### level 4 LINCS data by finding the size of the intersection of genes that
### have p-values below 0.0001, 0.001, 0.01, and 0.05. Computes the Fisher's
### test on these intersections, differentiating between cell lines and drugs.

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print "Usage: " + sys.argv[0] + "AFT_NUM pca/l1/exp"
        exit(1)
    aft_num = sys.argv[1]
    method = sys.argv[2]

    assert (method in ['pca', 'exp'])

    out = open('./results/compare_lincs_Aft_%s_and_%s.txt' % (aft_num, method), 'w')
    out.write('p-value\tdrug\tinter\tlincs\t%s' % method)
    out.write('\tneither\to_r\tfish-p\n')
    # Results Processing
    for p_threshold in [0.0001, 0.001, 0.01, 0.05]:
        # print 'Extracting data from %s output file...' % method
        pathways = set([])
        res_below_low_p = set([])
        # For if we want to see the intersections by drug.
        res_drug_dct = {}

        results_file = open('./results/top_pathways_%s.txt' % method, 'r')
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

            if float(score) <= p_threshold:
                res_below_low_p.add((drug, path))
                if drug not in res_drug_dct:
                    res_drug_dct[drug] = set([path])
                else:
                    res_drug_dct[drug].add(path)
        results_file.close()

        # print 'Extracting level 4 LINCS top pathways...'
        f = open('./results/top_pathways_lincs_Aft_%s.txt' % aft_num, 'r')
        lincs_cell_line_dct = {}
        lincs_drug_dct = {}
        for i, line in enumerate(f):
            # Skip first two lines, they're just results summaries.
            if i < 2:
                continue

            # Disregard x, y, z.
            drug, cell_line, path, score, x, y, z = line.strip().split('\t')
            pathways.add(path)
            if cell_line not in lincs_cell_line_dct:
                lincs_cell_line_dct[cell_line] = set([])
            
            if float(score) <= p_threshold:
                lincs_cell_line_dct[cell_line].add((drug, path))

                if drug not in lincs_drug_dct:
                    lincs_drug_dct[drug] = set([path])
                else:
                    lincs_drug_dct[drug].add(path)
        f.close()

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

            if o_r > 1.0 and p_val < 0.01:
                print drug

            out.write('%f\t%s\t%d\t%d\t' % (p_threshold, drug, lincs_and_res, lincs_not_res))
            out.write('%d\t%d\t%g\t%g\n' % (res_not_lincs, neither, o_r, p_val))
    out.close()