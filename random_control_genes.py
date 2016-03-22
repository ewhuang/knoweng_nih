### Author: Edward Huang

import file_operations
import fisher_test
import math
import random
from subtract_superdrug_from_pathways import get_superdrug_pathways
import sys

### This script randomly samples pathways for each drug, and performs a
### Fisher's test between those genes and LINCS. It does this 1,000 times.

# For a given set of superdrug pathways, finds the number of pathways for
# each drug-pathway pair that are better than the threshold.
def get_exp_paths_per_drug(random_p_thresh, superdrug_pathways):
    f = open('./results/top_pathways_exp_hgnc.txt', 'r')
    exp_paths_per_drug = {}
    for i, line in enumerate(f):
        # Skip header lines.
        if i < 2:
            continue
        drug, path, score = line.strip().split('\t')[:3]
        # We need to skip corresponding superdrug pathways for each
        # threshold.
        if score == '[]' or path in superdrug_pathways:
            continue
        # If a Fisher's p-value is below the threshold, add it to the
        # dictionary.
        score = float(score)
        if score <= random_p_thresh:
            if drug in exp_paths_per_drug:
                exp_paths_per_drug[drug] += [score]
            else:
                exp_paths_per_drug[drug] = [score]
    f.close()
    return exp_paths_per_drug

# Sample N pathways from the set of pathways not significant for the
# superdrug, where N is the number of pathways for expression below
# the threshold.
def sample_drug_paths(exp_paths_per_drug, pathways):
    drug_paths = {}
    for drug in exp_paths_per_drug:
        num_below_thresh = len(exp_paths_per_drug[drug])
        drug_paths[drug] = set(random.sample(pathways, num_below_thresh))
    return drug_paths

if __name__ == '__main__':
    if (len(sys.argv) != 2):
        print "Usage: " + sys.argv[0] +" NUM_RUNS"
        exit(1)
    NUM_RUNS = int(sys.argv[1])

    p_thresh_range = [0.001, 0.005, 0.01, 0.05, 0.1]

    # Extract the NCI pathway data.
    nci_path_dct, nci_genes = file_operations.get_nci_path_dct()
    pathways = set(nci_path_dct.keys())
    # Extract LINCS data.
    lincs_drug_path_dct = file_operations.get_lincs_drug_path_dct()

    table_dct = {}
    # Run script by NUM_RUNS loops each.
    for run in range(NUM_RUNS):
        print run
        for random_p_thresh in p_thresh_range:
            # Gets the pathways significant for the superdrug for the given
            # threshold.
            superdrug_pathways = set(get_superdrug_pathways(random_p_thresh))
            exp_paths_per_drug = get_exp_paths_per_drug(random_p_thresh,
                superdrug_pathways)
            # Get the set of pathways that are not in the superdrug pathways.
            pathways_minus_superdrug = pathways.difference(superdrug_pathways)
            sampled_drug_paths = sample_drug_paths(exp_paths_per_drug,
                pathways_minus_superdrug)

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
                    sampled_paths = sampled_drug_paths[drug]
                    lincs_below_p = lincs_drug_dct[drug]
                    lincs_and_samp = len(lincs_below_p.intersection(
                        sampled_paths))
                    lincs_not_samp = len(lincs_below_p.difference(
                        sampled_paths))
                    samp_not_lincs = len(sampled_paths.difference(
                        lincs_below_p))
                    neither = len(pathways) - len(lincs_below_p.union(
                        sampled_paths))

                    f_table = [[lincs_and_samp, lincs_not_samp],
                        [samp_not_lincs, neither]]
                    ft = fisher_test.FishersExactTest(f_table)
                    p_val = ft.two_tail_p()
                    if p_val <= 0.05:
                        if (random_p_thresh, lincs_p_thresh) in table_dct:
                            table_dct[(random_p_thresh, lincs_p_thresh)] += 1
                        else:
                            table_dct[(random_p_thresh, lincs_p_thresh)] = 1
    
    out = open('./results/random_table_minus_superdrug.txt', 'w')
    for random_p_thresh in p_thresh_range:
        for lincs_p_thresh in p_thresh_range:
            # If we recorded nothing for a pair of thresholds, then we just
            # write a 0.
            if (random_p_thresh, lincs_p_thresh) not in table_dct:
                out.write('0\t')
            else:
                out.write('%g\t' % (float(
                    table_dct[(random_p_thresh, lincs_p_thresh)]) / NUM_RUNS))
        out.write('\n')
    out.close()