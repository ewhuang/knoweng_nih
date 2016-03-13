### Author: Edward Huang

from collections import OrderedDict

### This script summarizes the compare_lincs_Aft_NUM_and_method_hgnc.txt files
### with the following table format:
### method_p/fisher-p  0.001 0.005 0.01 0.05 0.1
### 0.001
### 0.005
### ...

# method_p_range = [0.0001, 0.001, 0.01, 0.05]
p_thresh_range = [0.001, 0.005, 0.01, 0.05, 0.1]
COMPARISON_P_THRESH = 0.01

def summarize_file_and_write(in_filename):
    f = open(in_filename, 'r')
    table_dct = OrderedDict({})
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        method_p, lincs_p, comparison_p = float(line[0]), float(line[1]), float(line[-1])
        if comparison_p < COMPARISON_P_THRESH:
            if (method_p, lincs_p) in table_dct:
                table_dct[(method_p, lincs_p)] += 1
            else:
                table_dct[(method_p, lincs_p)] = 1
    f.close()

    out = open(in_filename[:10] + 'summ_' + in_filename[10:], 'w')
    out.write('\t' + '\t'.join(map(str, p_thresh_range)) + '\n')
    for method_p in p_thresh_range:
        out.write(str(method_p))
        for lincs_p in p_thresh_range:
            if (method_p, lincs_p) not in table_dct:
                out.write('\t0')
            else:
                out.write('\t%d' % table_dct[(method_p, lincs_p)])
        out.write('\n')
    out.close()


if __name__ == '__main__':
    # Expression summaries.
    base = './results/compare_lincs_Aft_3_and_'
    exp_fname = base + 'exp_hgnc.txt'
    summarize_file_and_write(exp_fname)

    top_k = 250

    # Embedding summaries.
    for method in ['genetic', 'literome', 'sequence', 'ppi']: # SKIP PPI RIGHT NOW
        if method == 'ppi':
            dimensions = map(str, [50, 100, 500, 1000, 1500, 2000])
        else:
            dimensions = map(str, [50, 100, 500])
        for dim in dimensions:
            for suffix in ['U', 'US']:
                entity_vector_dct = OrderedDict({})

                extension = '%s_0.8.%s' % (dim, suffix)
                embedding_fname = base + method + '_' + extension +\
                    '_top_%d_hgnc.txt' % (top_k)
                summarize_file_and_write(embedding_fname)

    # L1 summaries.
    print 'Currently not comparing L1...'
    l1_fname = base + 'l1_hgnc.txt'
    # summarize_file_and_write(l1_fname)

    # Compare the expression summaries with embedding summaries.
    exp_file = open('./results/summ_compare_lincs_Aft_3_and_exp_hgnc.txt', 'r')
    exp_table = []
    for i, line in enumerate(exp_file):
        if i == 0:
            continue
        line = map(float, line.split())[1:]
        exp_table += line
    exp_file.close()

    # Get the embedding summaries.
    for method in ['genetic', 'literome', 'sequence','ppi']: # SKIP PPI RIGHT NOW
        if method == 'ppi':
            dimensions = map(str, [50, 100, 500, 1000, 1500, 2000])
        else:
            dimensions = map(str, [50, 100, 500])
        for num in dimensions:
            for suffix in ['U', 'US']:
                entity_vector_dct = OrderedDict({})

                extension = '%s_0.8.%s' % (num, suffix)
                embedding_fname = './results/summ_compare_lincs_Aft_3_and_' +\
                    method + '_' + extension + '_top_%d_hgnc.txt' % (top_k)
                embedding_table = []
                f = open(embedding_fname, 'r')
                for i, line in enumerate(f):
                    if i == 0:
                        continue
                    line = map(float, line.split())[1:]
                    embedding_table += line
                f.close()
                num_better_than_expression = 0
                for i, val in enumerate(embedding_table):
                    if val > exp_table[i]:
                        num_better_than_expression += 1
                if num_better_than_expression >= 8:
                    print embedding_fname[10:], num_better_than_expression