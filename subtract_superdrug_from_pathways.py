### Author: Edward Huang

from collections import OrderedDict
import sys

### Goes through the top pathways for any given method, and takes out the top
### pathways for the superdrug that have a p-value below the given threshold.
sub_dir = './results/'

def get_top_pathways(method):
    top_pathways = OrderedDict({})
    if method == 'exp':
        filename = 'top_pathways_exp_hgnc.txt'
    elif method == 'genetic':
        filename = 'embedding/genetic_top_pathways_500_0.8.US_top_250.txt'
    f = open(sub_dir + filename, 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        line = line.strip().split('\t')
        assert len(line) == 7 or len(line) == 3
        drug, path, score = line[:3]
        assert (drug, path) not in top_pathways
        top_pathways[(drug, path)] = score
    f.close()
    return top_pathways

def get_superdrug_pathways(superdrug_p_value):
    superdrug_pathways = []
    f = open('./results/superdrug_pathway_p_values.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        assert len(line) == 2
        path, p_val = line[0], float(line[1])
        if p_val <= superdrug_p_value:
            superdrug_pathways += [path]
    f.close()
    return superdrug_pathways

def main():
    if (len(sys.argv) != 3):
        print "Usage: " + sys.argv[0] +" method superdrug_p_value"
        exit(1)
    method = sys.argv[1]
    superdrug_p_value = float(sys.argv[2])
    assert method in ['exp', 'genetic']

    top_pathways = get_top_pathways(method)
    superdrug_pathways = get_superdrug_pathways(superdrug_p_value)

    out = open('./results/top_pathways_%s_subtract_superdrug.txt' %
        method, 'w')
    # Remove the superdrug pathways from the expression pathway data.
    for (drug, path) in top_pathways:
        if path in superdrug_pathways:
            continue
        out.write('%s\t%s\t%s\n' % (drug, path, top_pathways[(drug, path)]))
    out.close()


if __name__ == '__main__':
    main()