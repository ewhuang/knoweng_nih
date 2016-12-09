### Author: Edward Huang

from collections import OrderedDict
import os
import sys

### This script summarizes the compare_lincs_Aft_NUM_and_method_hgnc.txt files
### with the following table format:
### method_p/fisher-p  0.001 0.005 0.01 0.05 0.1
### 0.001
### 0.005
### ...

comparison_p_thresh = 0.05
p_thresh_range = [0.001, 0.005, 0.01, 0.05, 0.1]
results_folder = './results'

def count_below_comparison_p(in_filename, comparison_p_thresh):
    '''
    Reads a comparison file between LINCS and a particular method, dictated by
    in_filename. Counts the number of drugs for each pair of p-values that
    are better than comparison_p_thresh.
    Returns a dictionary.
    Key: (method p-value, lincs p-value) pair -> (float, float)
    Value: number of drugs better than comparison_p_thresh -> int
    '''
    table_dct = OrderedDict({})

    f = open('%s/lincs_z%s_max%s_comparison_files/%s' % (results_folder,
        lincs_z, lincs_max_num, in_filename), 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue

        line = line.split()
        method_p, lincs_p, comparison_p = map(float, (line[0], line[1], 
            line[-1]))
        drug = line[2]

        # Update the number of p-values below comparison_p_thresh.
        if comparison_p < comparison_p_thresh:
            if (method_p, lincs_p) in table_dct:
                table_dct[(method_p, lincs_p)] += 1
            else:
                table_dct[(method_p, lincs_p)] = 1
    f.close()

    return table_dct

def write_comparison_table(in_filename, table_dct):
    '''
    Write out the table created by count_below_comparison_p().
    '''
    out = open('%s/summ_%s' % (out_folder, in_filename), 'w')
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

def summarize_file_and_write(in_filename, comparison_p_thresh):
    table_dct = count_below_comparison_p(in_filename, comparison_p_thresh)
    write_comparison_table(in_filename, table_dct)

def read_correlation_table(corr_kw_fname):
    '''
    Reads the correlation with KW summary previously written by
    summarize_file_and_write.
    Return: list of table values -> list(float)
    '''
    correlation_table = []

    correlation_file = open('%s/summ_%s' % (out_folder, corr_kw_fname), 'r')
    for i, line in enumerate(correlation_file):
        if i == 0:
            continue
        line = map(float, line.split())[1:]
        correlation_table += line
    correlation_file.close()

    return correlation_table

def write_best_files(base_fname, embed_top_k, correlation_table):
    '''
    Determines what the best NETPATH files are, as compared to the correlation
    with KW files.
    '''
    out = open('%s/best_files.txt' % (out_folder), 'w')
    out.write('filename\tnum_better_than_kw\n')
    # Get the embedding summaries.
    for method in ['genetic', 'literome', 'sequence','ppi']:
        dimensions = map(str, [50, 100, 500])
        if method == 'ppi':
            dimensions += map(str, [1000, 1500, 2000])
        for num in dimensions:
            for suffix in ['U', 'US']:
                entity_vector_dct = OrderedDict({})

                extension = '%s_0.8.%s' % (num, suffix)
                embedding_fname = '%s/summ_%s%s_%s_embed%d.txt' % (
                    out_folder, base_fname, method, extension, embed_top_k)

                embedding_table = []
                f = open(embedding_fname, 'r')
                for i, line in enumerate(f):
                    if i == 0:
                        continue
                    line = map(float, line.split())[1:]
                    embedding_table += line
                f.close()

                num_better_than_correlation = 0
                for i, val in enumerate(embedding_table):
                    if val > correlation_table[i]:
                        num_better_than_correlation += 1
                    elif val == correlation_table[i]:
                        num_better_than_correlation += 0.5
                if num_better_than_correlation >= 8:
                    out.write('%s\t%d\n' % (embedding_fname[
                        embedding_fname.index('summ'):],
                        num_better_than_correlation))

if __name__ == '__main__':
    if (len(sys.argv) != 4):
        print "Usage: " + sys.argv[0] + " lincs_z lincs_max_num embed_top_k"
        exit(1)
    global lincs_z, lincs_max_num    
    lincs_z, lincs_max_num,  = sys.argv[1], sys.argv[2]
    assert lincs_max_num.isdigit()
    embed_top_k = int(sys.argv[3])

    # Create the outfile folder, if necessary.
    global out_folder
    out_folder = '%s/comparison_summary_%g_z%s_max%s' % (results_folder,
        comparison_p_thresh, lincs_z, lincs_max_num)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # File name for all files.
    base_fname = 'compare_lincs_Aft_3_and_'
    corr_kw_fname = 'compare_lincs_Aft_3_and_corr_kw.txt'
    
    # Summarizing correlation with KW file.
    summarize_file_and_write(corr_kw_fname, comparison_p_thresh)

    # Summarizing NETPATH files.
    for method in ['genetic', 'literome', 'sequence', 'ppi']:
        dimensions = map(str, [50, 100, 500])
        if method == 'ppi':
            dimensions += map(str, [1000, 1500, 2000])
        for dim in dimensions:
            for suffix in ['U', 'US']:
                extension = '%s_0.8.%s' % (dim, suffix)
                embedding_fname = '%s%s_%s_embed%d.txt' % (base_fname,
                    method, extension, embed_top_k)
                summarize_file_and_write(embedding_fname, comparison_p_thresh)

    # Comparing correlation with KW to NETPATH files.
    correlation_table = read_correlation_table(corr_kw_fname)

    write_best_files(base_fname, embed_top_k, correlation_table)