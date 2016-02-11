### Author: Edward Huang

### This script does some basic preprocessing to prepare the l1 data for
### comparison with LINCS data.

if __name__ == '__main__':
    f = open('./results/nci_pathway_ranking_name.txt', 'r')
    out = open('./results/top_pathways_l1_hgnc.txt', 'w')
    # Write in the filler lines
    out.write('filler\ndrug\tpath\tscore\n')
    for line in f:
        drug, path, score, ignore = line.strip().split('\t')
        out.write('%s\t%s\t%s\n' % (drug.upper(), path, score))
    out.close()
    f.close()