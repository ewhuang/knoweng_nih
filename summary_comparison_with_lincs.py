### Author: Edward Huang

from collections import OrderedDict
import sys

### This script summarizes the compare_lincs_Aft_NUM_and_method_hgnc.txt files
### with the following table format:
### method_p/fisher-p  0.001 0.005 0.01 0.05 0.1
### 0.001
### 0.005
### ...

p_thresh_range = [0.001, 0.005, 0.01, 0.05, 0.1]
main_folder = './results/'
subfolder = 'comparison_summary_'

def summarize_file_and_write(in_filename, comparison_p_thresh):
    f = open(main_folder + 'lincs_comparison_files/' + in_filename, 'r')
    table_dct = OrderedDict({})
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        method_p, lincs_p, comparison_p = (float(line[0]), float(line[1]),
            float(line[-1]))
        # Count the number of p-values below the current threshold.
        if comparison_p < comparison_p_thresh:
            if (method_p, lincs_p) in table_dct:
                table_dct[(method_p, lincs_p)] += 1
            else:
                table_dct[(method_p, lincs_p)] = 1
    f.close()

    # Write out the table.
    out = open('%s%s%g/summ_%s' % (main_folder, subfolder, comparison_p_thresh,
        in_filename), 'w')
    out.write('\t' + '\t'.join(map(str, p_thresh_range)) + '\n')
    for method_p in p_thresh_range:
        out.write(str(method_p))
        for lincs_p in p_thresh_range:
            # If we recorded nothing for a pair of thresholds, then we just
            # write a 0.
            if (method_p, lincs_p) not in table_dct:
                out.write('\t0')
            else:
                out.write('\t%d' % table_dct[(method_p, lincs_p)])
        out.write('\n')
    out.close()


if __name__ == '__main__':
    if (len(sys.argv) < 2):
        print "Usage: " + sys.argv[0] +" TOP_K"
        exit(1)
    top_k = int(sys.argv[1])

    # File name for all files.
    base_fname = 'compare_lincs_Aft_3_and_'
    exp_fname = base_fname + 'exp_hgnc.txt'
    
    print 'Currently not comparing L1...'

    for comparison_p_thresh in p_thresh_range:
        # File name for correlation.
        summarize_file_and_write(exp_fname, comparison_p_thresh)

        # File names for embedding networks.
        for method in ['genetic', 'literome', 'sequence', 'ppi']:
            if method == 'ppi':
                dimensions = map(str, [50, 100, 500, 1000, 1500, 2000])
            else:
                dimensions = map(str, [50, 100, 500])
            for dim in dimensions:
                for suffix in ['U', 'US']:
                    extension = '%s_0.8.%s' % (dim, suffix)
                    embedding_fname = base_fname + method + '_' + extension +\
                        '_top_%d_hgnc.txt' % (top_k)
                    summarize_file_and_write(embedding_fname,
                        comparison_p_thresh)

        # L1 summaries.
        # l1_fname = base_fname + 'l1_hgnc.txt'
        # # summarize_file_and_write(l1_fname)

        # Compare the expression summaries with embedding summaries.
        exp_file = open('%s%s%g/summ_%s' % (main_folder, subfolder,
            comparison_p_thresh, exp_fname), 'r')
        exp_table = []
        for i, line in enumerate(exp_file):
            if i == 0:
                continue
            line = map(float, line.split())[1:]
            exp_table += line
        exp_file.close()

        out = open('%s%s%g/best_files.txt' % (main_folder, subfolder,
            comparison_p_thresh), 'w')
        out.write('filename\tnum_better_than_correlation\n')
        # Get the embedding summaries.
        for method in ['genetic', 'literome', 'sequence','ppi']:
            if method == 'ppi':
                dimensions = map(str, [50, 100, 500, 1000, 1500, 2000])
            else:
                dimensions = map(str, [50, 100, 500])
            for num in dimensions:
                for suffix in ['U', 'US']:
                    entity_vector_dct = OrderedDict({})

                    extension = '%s_0.8.%s' % (num, suffix)
                    embedding_fname = '%s%s%g/summ_%s%s_%s_top_%d_hgnc.txt' % (main_folder,
                        subfolder, comparison_p_thresh, base_fname, method,
                        extension, top_k)
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
                        if val == exp_table[i]:
                            num_better_than_expression += 0.5
                    if num_better_than_expression >= 8:
                        out.write('%s\t%f\n' % (embedding_fname[embedding_fname.index('summ'):],
                            num_better_than_expression))