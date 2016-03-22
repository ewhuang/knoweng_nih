### Author: Edward Huang

from collections import OrderedDict
import sys

### Goes through the top pathways for any given method, and takes out the top
### pathways for the superdrug that have a p-value below the given threshold.

def get_exp_pathways(method):
    exp_pathways = OrderedDict({})
    if method == 'exp':
        filename = './results/top_pathways_exp_hgnc.txt'
    f = open(filename, 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        line = line.strip().split('\t')
        assert len(line) == 7
        drug, path, score = line[:3]
        assert (drug, path) not in exp_pathways
        exp_pathways[(drug, path)] = score
    f.close()
    return exp_pathways

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
    assert method in ['exp']

    exp_pathways = get_exp_pathways(method)
    superdrug_pathways = get_superdrug_pathways(superdrug_p_value)

    out = open('./results/top_pathways_exp_hgnc_subtract_superdrug.txt', 'w')
    # Remove the superdrug pathways from the expression pathway data.
    for (drug, path) in exp_pathways:
        if path in superdrug_pathways:
            continue
        out.write('%s\t%s\t%s\n' % (drug, path, exp_pathways[(drug, path)]))
    out.close()


if __name__ == '__main__':
    main()