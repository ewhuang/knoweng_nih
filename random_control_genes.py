### Author: Edward Huang

import file_operations
import fisher_test
import math
import random
import sys

### This script randomly samples pathways for each drug, and performs a
### Fisher's test between those genes and LINCS. It does this 1,000 times.

def min_p_exp(p_val_lst):
    return 1 - math.pow((1 - min(p_val_lst)), len(p_val_lst))

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print "Usage: " + sys.argv[0] +" RUN_NUM"
        exit(1)
    run_num = int(sys.argv[1])

    p_thresh_range = [0.001, 0.005, 0.01, 0.05, 0.1]

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
        # If a Fisher's p-value is below 0.1, add it to the dictionary.
        score = float(score)
        if score <= 0.1:
            if drug in exp_paths_per_drug:
                exp_paths_per_drug[drug] += [score]
            else:
                exp_paths_per_drug[drug] = [score]
    f.close()

    # Extract the NCI pathway data.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()

    # For threshold 0.1, we sample pathways for each drug. The number of
    # pathways corresponds to the number of pathways below the threshold for
    # expression.
    pathways = nci_path_dct.keys()
    sampled_drug_paths = {}
    for drug in exp_paths_per_drug:
        exp_p_values = exp_paths_per_drug[drug]
        num_below_thresh = len([1 for p_val in exp_p_values if p_val <= 0.1])
        sampled_drug_paths[drug] = random.sample(pathways, num_below_thresh)

    # print 'Extracting level 4 LINCS top pathways...'
    f = open('./results/top_pathways_lincs_Aft_3_hgnc.txt', 'r')
    lincs_drug_path_dct = {}
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        drug, cell_line, path, score = line.strip().split('\t')[:4]
        score = float(score)
        if (drug, path) not in lincs_drug_path_dct:
            lincs_drug_path_dct[(drug, path)] = [score]
        else:
            lincs_drug_path_dct[(drug, path)] += [score]
    f.close()
    # Aggregate p-values by cell lines.
    for (drug, path) in lincs_drug_path_dct:
        p_val_lst = lincs_drug_path_dct[(drug, path)]
        # Change the function for other aggregation functions.
        lincs_drug_path_dct[(drug, path)] = min_p_exp(p_val_lst)

    table_dct = {}
    # Run script by 100 loops each.
    for run in range(100):
        # For each threshold, we get a subset of sampled pathways, and then compute
        # Fisher's with LINCS top pathways.
        for random_p_thresh in p_thresh_range:
            print random_p_thresh
            for drug in exp_paths_per_drug:
                num_below_thresh = len([1 for p_val in exp_p_values if p_val <=
                    random_p_thresh])
                # Get a subset of the sampled pathways with size equal to the
                # number of pathways below the current threshold for expression.
                sampled_paths = set(sampled_drug_paths[drug][:num_below_thresh])

                # Fisher's test between the sampled paths and LINCS top pathways,
                # chosen also by a p-value threshold.
                for lincs_p_thresh in p_thresh_range:
                    lincs_drug_dct = {}
                    for (drug, path) in lincs_drug_path_dct:
                        if drug not in exp_paths_per_drug:
                            continue
                        if lincs_drug_path_dct[(drug, path)] <= lincs_p_thresh:
                            if drug not in lincs_drug_dct:
                                lincs_drug_dct[drug] = set([path])
                            else:
                                lincs_drug_dct[drug].add(path)

                    for drug in lincs_drug_dct:
                        lincs_below_p = lincs_drug_dct[drug]
                        lincs_and_samp = len(lincs_below_p.intersection(sampled_paths))
                        lincs_not_samp = len(lincs_below_p.difference(sampled_paths))
                        samp_not_lincs = len(sampled_paths.difference(lincs_below_p))
                        neither = len(pathways) - len(lincs_below_p.union(sampled_paths))

                        f_table = [[lincs_and_samp, lincs_not_samp], [samp_not_lincs, neither]]
                        ft = fisher_test.FishersExactTest(f_table)
                        p_val = ft.two_tail_p()
                        if p_val <= 0.05:
                            if (random_p_thresh, lincs_p_thresh) in table_dct:
                                table_dct[(random_p_thresh, lincs_p_thresh)] += 1
                            else:
                                table_dct[(random_p_thresh, lincs_p_thresh)] = 1
    
    out = open('./results/random_control/random_table_%d.txt' % run_num, 'w')
    for random_p_thresh in p_thresh_range:
        for lincs_p_thresh in p_thresh_range:
            # If we recorded nothing for a pair of thresholds, then we just
            # write a 0.
            if (random_p_thresh, lincs_p_thresh) not in table_dct:
                out.write('0\t')
            else:
                out.write('%d\t' % table_dct[(random_p_thresh, lincs_p_thresh)])
        out.write('\n')
    out.close()