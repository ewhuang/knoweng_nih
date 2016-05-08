### Author: Edward Huang

import file_operations
import time

### Various scripts to perform case studies in the paper.
### Run time: 6 seconds.

COMPARISON_P_THRESH = 0.05
main_folder = './results/lincs_comparison_files/'

def get_drugs_given_p_values(fname, target_method_p, target_lincs_p):
    '''
    Gets the drugs from a given filename. Takes in two p-values. One is the
    LINCS p-value and the other is the method p-value.
    '''
    significant_drugs = []
    f = open(fname, 'r')
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        method_p, lincs_p, comparison_p = (float(line[0]), float(line[1]),
            float(line[-1]))
        drug = line[2]
        # Count the number of p-values below the current threshold.
        if comparison_p < COMPARISON_P_THRESH and method_p == target_method_p\
            and lincs_p == target_lincs_p:
            significant_drugs += [drug]
    f.close()
    return set(significant_drugs)

def get_drug_differences(fname_a, fname_b, method_p, lincs_p):
    '''
    Gets the drug differences between the first file and the second file for a
    given slot in the summary table.
    '''
    a_drugs = get_drugs_given_p_values(fname_a, method_p, lincs_p)
    b_drugs = get_drugs_given_p_values(fname_b, method_p, lincs_p)
    return list(b_drugs.difference(a_drugs))

## For the drugs, get their p-values.
def get_p_values_given_drugs(fname, drugs):
    '''
    Given a top pathways filename, get the p-value for the list of drug-
    pathway pairs.
    '''
    name_to_brd_dct = get_name_to_brd_dct()
    drug_pathway_dict = {}
    f = open(fname, 'r')
    for i, line in enumerate(f):
        # Skip header line.
        if i < 2:
            continue
        line = line.strip().split('\t')
        drug, path, score = line[:3]

        # Convert back to BRD ID if 
        if drug in name_to_brd_dct:
            drug = name_to_brd_dct[drug]

        if drug in drugs:
            drug_pathway_dict[(drug, path)] = score
    f.close()
    return drug_pathway_dict

def get_name_to_brd_dct():
    name_to_brd_dct = {}
    f = open('./data/brd_to_name.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        name, cpd_id, brd_id = line.strip().split('\t')
        name_to_brd_dct[name] = brd_id
    f.close()
    return name_to_brd_dct

def get_brd_drug_to_name_dct():
    brd_drug_to_name_dct = {}
    f = open('./data/brd_to_name.txt', 'r')
    for i, line in enumerate(f):
        # Skip header line.
        if i == 0 or line == '\n':
            continue
        line = line.strip().split('\t')
        assert len(line) == 3
        drug_name = line[0]
        drug_id = line[2]
        assert drug_id not in brd_drug_to_name_dct
        brd_drug_to_name_dct[drug_id] = drug_name
    f.close()
    return brd_drug_to_name_dct

def main():
    brd_drug_to_name_dct = get_brd_drug_to_name_dct()

    drugs = get_drug_differences('%scompare_lincs_Aft_3_and_exp_hgnc.txt' % (
        main_folder), ('%scompare_lincs_Aft_3_and_ppi_50_0.8.U_top_250_hgnc.txt'
        % (main_folder)), 0.01, 0.005)

    netpath_drug_path_dct = get_p_values_given_drugs('./results/embedding/' + 
        'inverse_rankings/inverse_ppi_top_pathways_50_0.8.U_top_250.txt', drugs)

    lincs_drug_path_dct = file_operations.get_lincs_drug_path_dct()

    exp_drug_path_dct = get_p_values_given_drugs('./results/top_pathways_exp' +
        '_hgnc.txt', drugs)

    # Write out the relevant drug information to file.
    out = open('./results/exp_netpath_drug_difference.txt', 'w')
    out.write('drug\tpath\tlincs_p_value\tpearson_p_value\tnetpath_inverse_p_' +
        'value\n')
    for key in netpath_drug_path_dct:
        drug, path = key
        lincs_p_value = lincs_drug_path_dct[key]
        pearson_p_value = exp_drug_path_dct[key]
        netpath_p_value = netpath_drug_path_dct[key]
        out.write('%s\t%s\t%g\t%s\t%s\n' % (brd_drug_to_name_dct[drug], path,
            lincs_p_value, pearson_p_value, netpath_p_value))
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))